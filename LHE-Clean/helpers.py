"""
Helper functions and objects
"""

import numpy as np
import argparse
import os
import glob
import subprocess

import matrix_element_code


#-----------------------------------OBJECTS------------------------------------

class Particle():
	"""
	A class to store all the information associated to a particle
	"""	

	def __init__(self, in_line):
		"""
		Line information is: [PID, STATUS, MOTHER(S), COLOR(S), 4MOM(PX,PY,PZ,E), MASS, LIFETIME, SPIN]
		"""
		line = in_line.split()
		self.pid = int(line[0])
		self.status = int(line[1])
		self.mother1 = int(line[2])
		self.mother2 = int(line[3])
		self.color1 = int(line[4])
		self.color2 = int(line[5])
		self.px = float(line[6])
		self.py = float(line[7])
		self.pz = float(line[8])
		self.e = float(line[9])
		self.mass = float(line[10])
		self.lifetime = float(line[11])
		self.hel = float(line[12])

	def get_string(self):
		"""
		Return a formatted string with all the particle information for writing to output LHE file
		"""

		temp = [ str(self.pid), str(self.status), str(self.mother1), str(self.mother2), str(self.color1), str(self.color2),
		"{:.10e}".format(self.px), "{:.10e}".format(self.py), "{:.10e}".format(self.pz), "{:.10e}".format(self.e), "{:.10e}".format(self.mass),
		"{:.6e}".format(self.lifetime), "{:.6e}".format(self.hel), ]
		string = "\t" + "\t".join(temp) + "\n"
		return(string)

class Event():
	"""
	A container to store all the info of a yy->tautau event
	Line is formatted as [PID, STATUS, MOTHER(S), COLOR(S), 4MOM(PX,PY,PZ,E), MASS, LIFETIME, SPIN/HELICITY]
	"""
	def __init__(self, lines, weight):
		self.y1  = Particle(lines[0])
		self.y2  = Particle(lines[1])
		self.ta1 = Particle(lines[2])
		self.ta2 = Particle(lines[3])
		# Compute the scale -- avg pT of the taus
		self.scale = self.get_scale(self.ta1, self.ta2)
		# The event weight
		self.weight = float(weight)

	def get_scale(self, tau1, tau2):
		"""
		Get the average pT of the two taus in the event as the scale of the event
		"""
		px1, py1 = tau1.px, tau1.py
		px2, py2 = tau2.px, tau2.py
		avg_pt = (np.sqrt(px1**2+py1**2)+np.sqrt(px2**2+py2**2))/2.0
		return(avg_pt)

class InitBlock():
	"""
	A class that contains the integrated weight and error read from <init> block
	of an input LHE file
	"""

	def __init__(self, int_wt, int_err):
		self.int_wt  = float(int_wt)
		self.int_err = float(int_err)

#---------------------------------I/O HELPERS----------------------------------

def myprint(string, label, newline_above, newline_below):
	"""
	Print a string to terminal with in format "label:\tString"
	newline_above, newline_below do exactly what they say 
	"""
	if newline_above: print()
	print('%s:\t%s'%(label, string))
	if newline_below: print()

def get_file_names(args):
	"""
	Function to get all the file names from a given path to LHE files
	"""
	bad_input_flag = 0

	# Check if * is used in the path. if not, add it. This is needed for glob.glob()
	if args.input[-1] != '*': args.input += '*'

	# Sort the possible file names, select only args.nfiles files
	temp_files = sorted(glob.glob(args.input))
	if args.nfiles == -1: args.nfiles = len(temp_files)
	input_files = sorted(glob.glob(args.input))[args.skipfiles:args.skipfiles+args.nfiles]

	# Check that there actually are some files between skipfiles and skipfiles+nfiles
	if len(input_files) == 0: 
		bad_input_flag = 1 
		myprint('Number of files is 0. Please check that --skipfiles and --nfiles are consistent', 'ERROR', True, False)
		return(bad_input_flag, input_files)

	# Check if all the files are LHE files
	for i in range(len(input_files)):
		if input_files[i].split('.')[-1] != 'lhe': 
			bad_input_flag = 1
			myprint('Found a file without .lhe extension: %s'%(input_files[i]), 'ERROR', False, False)

	return(bad_input_flag, input_files)



def get_total_int_wt(init_blocks):
	"""
	Compute the average of a set of integrated weights from a given init_block list
	"""
	int_wts = [init_blocks[i].int_wt for i in range(len(init_blocks))]
	return(np.average(int_wts))

def get_total_int_err(init_blocks):
	"""
	Compute the total error of a set of integration errors from a given init_block list
	Total error is computed as 1/Nfiles*Sqrt(sum(errors^2))
	"""
	int_errs = np.array([init_blocks[i].int_err for i in range(len(init_blocks))])
	return( np.sqrt(np.sum(int_errs**2))/len(init_blocks) )

def write_header(outfile, init_block_merged):
	"""
	Use sed command to write the header to the output file
	"""

	# Get the merged xsec and error
	xsec = "{:.9e}".format(init_block_merged.int_wt)
	xsec_err = "{:.9e}".format(init_block_merged.int_err)

	# Define the LHEF output header lines
	output_header = [
		r"<LesHouchesEvents version=\"1.0\">\n" ,
		r"<!--\n" ,
		r"Gamma Gamma -> lepton+ lepton- events generated using StarLight \n",
		r"-->\n" ,
		r"<init>\n",
		r"13  -13  2.510000e+03  2.510000e+03  0  0  0  0  3  1\n" ,
		r"%s  %s  %s  9999\n"%(xsec, xsec_err, xsec),
		r"<\/init>\n"
	]

	# write the lines
	for line in output_header[::-1]: os.system(r"sed -i '1s/^/"+line+r"/' "+str(outfile))

def write_tail(outputFile):
	"""
	Write the tail lines to LHEF output file.. ust write a line to close the LHEF tag
	"""

	outputFile.write("</LesHouchesEvents>")

#---------------------------------EVENTLOOPER----------------------------------

def event_looper(args, file_names):
	"""
	Main event looping function for MG5 input files. Read input file event-by-event and do whatever needs to be done
	"""

	# Set a reading_event flag, a temp event holder
	reading_event = 0
	temp_event = [] 

	# list of all InitBlocks. This can be added to the beginning of the output file using bash magic
	init_blocks = []

	# Open the output file, if only conversion is requested. LHEF header is written at the end, using bash magic
	if args.output and not args.reweight: outfile = open(args.output, 'w')

	# if reweighting requested, open a new file to store the weights
	if args.reweight: rwt_output = open(args.output.split('.')[0]+'_weights.dat', 'w')
	new_weights = []
	rwt_factors = []

	event_counter = 0

	# Start looping over all input files
	for file_name in file_names:

		# open the file
		f = open(file_name, 'r')
		myprint('Reading events from %s'%(file_name), 'INFO', False, False)

		# holders for the init block
		reading_init  = 0
		temp_init  = []

		# Start reading the current input file
		for line in f:
			
			# Switch reading_event on/off based on which line in the input LHE file we are at
			# When reading_event is switched off, append only the "Event" object with lines with actual particle info
			if "<event>"  in line: reading_event = 1
			if "</event>" in line: 

				reading_event = 0
				# Get the event weight (i.e. cross section), add a new event object to the file
				weight = temp_event[1].split()[2]
				
				# Now we have the event to do all the work with!
				event = Event(temp_event[2:6], weight)

				# if no reweighting has been requested, then we can directly go to writing output LHE file
				if not args.reweight:

					# Define the event open and close tags
					event_open, event_close = "<event>\n", "</event>\n"

					# Switch the incoming particle PID to electrons if later decaying by Tauola (1) or Pythia (2) [no change]
					if args.format == 1:
						event.y1.pid, event.y2.pid = "-11", "11"

					# Define the event header
					event_header = "\t".join(["4",								# 4 particles
											  "1",								# process id
											  "{:.6e}".format(event.weight),	# event weight
											  "{:.6e}".format(event.scale),		# event scale
											  "7.297352e-3",					# alpha_EM value (hardcoded!)
											  "2.569093e-01"]) + "\n"			# alpha_S value (hardcoded but not used anywhere)

					# Define all the particle info lines
					line_photon1 = event.y1.get_string()
					line_photon2 = event.y2.get_string()
					line_lepton1 = event.ta1.get_string()
					line_lepton2 = event.ta2.get_string()

					# Write the event open tag and event header
					outfile.write(event_open)
					outfile.write(event_header)

					# Write the photon strings to the output file
					outfile.write(line_photon1)
					outfile.write(line_photon2)

					# Write the lepton strings to the output file
					outfile.write(line_lepton1)
					outfile.write(line_lepton2)

					# Write the event close line
					outfile.write(event_close)

				# if reweighting has been requested, we only write the event weights as a separate file (no LHE file)
				if args.reweight:

					# get the [new_weight, reweight_factor], write both to the rwt_output file
					new_weight, rwt_factor = matrix_element_code.reweight_event(event, args)
					new_weights.append(new_weight)
					rwt_factors.append(rwt_factor)
					rwt_output.write("%s %s\n"%("{:.8e}".format(new_weight), "{:.8e}".format(rwt_factor)))

				# clear the event holder, proceed
				temp_event = []

				# Status update
				event_counter += 1
				if event_counter%10000 == 0: print('INFO:\tAt event number %s'%(event_counter))

			# While reading is true, add a particle info line to "temp_event"
			if reading_event == 1:
				temp_event.append(line)

			# Same switch on/off logic to read the init block
			if "<init>" in line: reading_init = 1
			if "</init>" in line:
				reading_init = 0
				# Get the integrated weight (xsec) and the associated error
				integrated_weight = temp_init[2].split()[0]
				integration_error = temp_init[2].split()[1]

				# Store the weight and error as InitBlock object, to be used when writing LHE file header
				init_object = InitBlock(integrated_weight, integration_error)
				init_blocks.append(init_object)

			if reading_init == 1:
				temp_init.append(line)


		# close, finished processing the current file
		f.close()
		myprint('Finished reading events from %s'%(file_name), 'INFO', False, False)

	# merge the init blocks, get the total xsec and error
	total_int_wt  = get_total_int_wt(init_blocks)
	total_int_err = get_total_int_err(init_blocks)
	init_block_merged = InitBlock(total_int_wt, total_int_err)	

	new_weights = np.array(new_weights)
	rwt_factors = np.array(rwt_factors)

	# if reweighting was requested, no output file exists, so just break
	if not args.output: return()

	# if output file exists, do the remaining work
	if args.output and not args.reweight:

		# write the tail and close
		write_tail(outfile)
		outfile.close()

		# do bash magic and write the header
		write_header(args.output, init_block_merged)

	# if reweighting was requested, just close the file, print the new xsec and error
	if args.reweight: 
		
		rwt_output.close()
		
		# Define the old xsec and error
		old_xsec, old_xsec_err = init_block_merged.int_wt, init_block_merged.int_err

		# compute the new xsec and error
		# Error of the reweighting factor is defined as the sqrt(variance/sqrt(N))
		error_rwfct = np.round(np.sqrt(np.var(rwt_factors)/np.sqrt(event_counter)), 5)
		mean_rwfct  = np.round(np.average(rwt_factors), 5)

		new_xsec = np.average(new_weights)
		new_xsec_err = old_xsec_err*mean_rwfct + old_xsec*error_rwfct

		myprint('Old xsec and error (in nb): %s +- %s'%(np.round(old_xsec/1.0e6, 3), np.round(old_xsec_err/1.0e6,3)), 'INFO', True, False)
		myprint('New xsec and error (in nb): %s +- %s'%(np.round(new_xsec/1.0e6, 3), np.round(new_xsec_err/1.0e6,3)), 'INFO', False, True)

		del new_weights, rwt_factors
