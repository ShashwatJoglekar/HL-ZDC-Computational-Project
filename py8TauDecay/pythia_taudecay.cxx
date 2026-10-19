/*

Make sure the library paths are correct and $PYTHIADIR points to the install of desired Pythia8 version
Compilation command:

g++ pythia_taudecay.cxx -o pythia_taudecay \
-O2  -pedantic -W -Wall -Wshadow -fPIC                                 \
-I$HEPMCDIR/include -I$PYTHIADIR/include -I$PHOTOSDIR/include          \ 
-L$HEPMCDIR/lib  -Wl,-rpath,$HEPMCDIR/lib  -lHepMC                     \
-L$PYTHIADIR/lib -Wl,-rpath,$PYTHIADIR/lib -lpythia8 -ldl              \
-L$PHOTOSDIR/lib -Wl,-rpath,$PHOTOSDIR/lib -lPhotospp -lPhotosppHepMC         

To run on the example file "/path/to/my/lhe/file.lhe"

./pythia_taudecay /path/to/my/lhe/file

i.e. omit the ".lhe" extension. 

A new file, "/path/to/my/lhe/file.hepmc" will be created.

*/

// Standard headers
#include <fstream>

// pythia8 header files
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"

//PHOTOS header files
#include "Photos/Photos.h"
#include "Photos/PhotosHepMCParticle.h"
#include "Photos/PhotosHepMCEvent.h"

// namespaces
using namespace Pythia8; 
using namespace Photospp;

int main(int argc,char **argv) {

	// initilaize pythia object and tohepmc object
	Pythia pythia;
	HepMC::Pythia8ToHepMC ToHepMC;	

	// Read the input file name from command line argument
	// Also set the output .hepmc file name (inclusive of location, coming from command line)
	(void)argc;

	string fileName = argv[1];
	string inFileName = fileName; inFileName.append(".lhe");
	string outFileName = fileName; outFileName.append(".hepmc");

	// Tell pythia that input file is in lhe format, read input file. Setup Hepmc ascii dumper
	string pythiaReadStringInFile = "Beams:LHEF = "; 
	pythiaReadStringInFile.append(inFileName);
	pythia.readString("Beams:frameType = 4");
	pythia.readString(pythiaReadStringInFile);
	HepMC::IO_GenEvent ascii_io(outFileName, ios::out);

	// QED showering commands. Timelike --> FSR. Spacelike --> ISR (see [4])
	pythia.readString("TimeShower:QEDshowerByL = on");  
	pythia.readString("SpaceShower:QEDshowerByL = on");
	pythia.readString("TimeShower:QEDshowerByOther = off"); // see [3]

	// rseed
	pythia.readString("Random:setSeed = on");
	pythia.readString("Random:seed = 100");

	// force uncorrelated tau decays with external spins 
	pythia.readString("TauDecays:mode = 5");
	pythia.readString("TauDecays:externalMode = 0");

	// initialize pythia and internals
	pythia.init();

	// Initialize photos
	Photos::initialize();

	// Photos settings based on [1,2,3]
	Photos::setAlphaQED(0.00729735039);
	Photos::setInterference(true);
	Photos::setCorrectionWtForW(true);
	Photos::maxWtInterference(3.0);
	Photos::setMeCorrectionWtForW(false);
	Photos::setMeCorrectionWtForZ(false);
	Photos::setPairEmission(false);
	Photos::forceMassFrom4Vector(true);
	Photos::forceMassFromEventRecord(13);
	Photos::forceMassFromEventRecord(15);
	Photos::forceMass(11, 0.510998910);
	Photos::forceMassFromEventRecord(211);
	Photos::setTopProcessRadiation(false);
	Photos::createHistoryEntries(false, 3);
	Photos::setExponentiation(true);
	Photos::setInfraredCutOff(1.0E-7);
	Photos::setStopAtCriticalError(false);
	Photos::iniInfo();

	// Counters for number of photons added by PHOTOS
	long int photonAdded=0,twoAdded=0,moreAdded=0;

	// Parameters for charged particle filter
	double eta_range = 2.6;
	double pt_min = 3.0; // GeV

	// Counter for events that pass the filter
	long int filter_pass_ctr = 0;
	vector<long int> passed_event_numbers;

	// Begin event loop; generate until none left in input file.
	for (long int iEvent = 0; ; ++iEvent) {

		// Call pythia.next(). Exit if EOF reached
		if (!pythia.next()) {
			if (pythia.info.atEndOfFile()) { 
				break; 
			}
		}

		// initalize new hepmc event object, convert the pythia object to hepmc
		HepMC::GenEvent * HepMCEvt = new HepMC::GenEvent();
		ToHepMC.fill_next_event(pythia, HepMCEvt);

		// Counter for number of photons added by PHOTOS
		int buf = -HepMCEvt->particles_size();

		// Run PHOTOS on the event
		PhotosHepMCEvent evt(HepMCEvt);
		evt.process();

		// Set a counter for number of filter-passed charged particles. See references [8, 9]
		int n_filter_pass = 0;

		// Loop over all the particles in the HepMC event using the particle iterator
		for (HepMC::GenEvent::particle_const_iterator p = HepMCEvt->particles_begin(); p != HepMCEvt->particles_end(); ++p){

			// Skip the particle if it is not a stable particle
			if ((*p)->status() != 1) {continue;}

			// Get the pT and eta of the particle
			double temp_pt = (*p)->momentum().perp();
			double temp_eta = abs((*p)->momentum().eta());

			// figure out if the particle is charged based on the PID. This hardcoding is justified for pythia 
			// since these are the only particles that ever show up. If a rare decay pops up, check the warning message 
			// and modify the list of pids
			std::vector<int> l_charged_pids = {11, -11, 13, -13, 211, -211, 321, -321}; // elec, muon, pion+/-, kaon+/-
			std::vector<int> l_neutral_pids = {12, -12, 14, -14, 16, -16, 22, 130, 310, 311, 2112, -2112}; // nu, gamma, K0L/S,  K0, neutron, antineutron
			int pid = (*p)->pdg_id();
			bool is_charged, is_neutral;

			if (std::find(l_charged_pids.begin(), l_charged_pids.end(), pid) != l_charged_pids.end()){ is_charged = true; }
			else{ is_charged = false; }
			if (std::find(l_neutral_pids.begin(), l_neutral_pids.end(), pid) != l_neutral_pids.end()){ is_neutral = true; }
			else{ is_neutral = false; }

			// Print warning if particle is neither charged nor neutral
			if (!is_charged and !is_neutral){std::cout << "MYWARNING: rare particle with pdg_id " << pid << " unknown!!" << std::endl;  } 

			// Now filter the particle. If it passes, count the total
			if (is_charged and temp_pt >= pt_min and temp_eta <= eta_range){ n_filter_pass += 1; }
		}

		// If the number of filtered charged particles is 0 then skip the event. Else write the event record to output file
		if (n_filter_pass == 0){ delete HepMCEvt;  }
		else{

			// Count the number of photons added by PHOTOS
			buf += HepMCEvt->particles_size();
			if(buf==1){ photonAdded++;}
			else if(buf==2){ twoAdded++; }
			else if(buf>2){ moreAdded++; }

			// Increment the filter pass counter, append the event number to passed event number
			filter_pass_ctr += 1; 
			passed_event_numbers.push_back(iEvent);

			// Write the HEPMC file and delete it
			ascii_io << HepMCEvt; 
			delete HepMCEvt; 
		}
	}

	// Print how many photons added by photos, and the number of events that pass the filter
	std::cout << "Record of photons added to the event by PHOTOS: " << std::endl;
	std::cout << "Events with  1 photon  added: " << photonAdded << std::endl 
			  << "Events with  2 photons added: " << twoAdded << std::endl
			  << "Events with >2 photons added: " << moreAdded << std::endl;  
	std::cout << "Events which pass the filter: " << filter_pass_ctr << std::endl;

	// Write the passed event numbers to a file
	string passedEventsFileName = fileName; passedEventsFileName.append("_passed_filter.dat");
	std::ofstream f;
	f.open(passedEventsFileName);
	for (long unsigned int j = 0; j < passed_event_numbers.size(); j++) {f << passed_event_numbers[j] << std::endl;}
	f.close();

	return 0;   	
}