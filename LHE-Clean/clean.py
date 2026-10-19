"""

Python script to do various things with LHE files.
Main uses:
	-- Converting MG5 LHE files so that they are suitable for UPC Pythia/Tauola decays
	-- Reweighting parton level events for given atau and dtau

"""

import argparse
import numpy as np
from time import ctime

import helpers

def main():
	"""
	Main!
	"""

	# Argument parser
	parser = argparse.ArgumentParser()

	# input/output files/directories
	parser.add_argument("-i", "--input", required=True, type=str, help="Inpute LHE file / path to input LHE file(s)")
	parser.add_argument("-o", "--output", required=False, type=str, help="Output LHE file / path to output LHE file" )

	# Number of files to consider/skip
	parser.add_argument("-n", "--nfiles", default=-1, type=int, required=False, help="Number of input LHE files to process. Default -1 = all files")
	parser.add_argument("-s", "--skipfiles", required=False, default=0, type=int, help="Number of input files to skip" )

	# Format for LHE conversion, tauola/pythia decays
	parser.add_argument("-f", "--format", type=int, required=False, default=2, help="Taudecay format for output LHE files. 1=Tauola, 2=Pythia (default)")

	# Arguments for reweighting
	parser.add_argument("-r", "--reweight", action='store_true', help="Reweight events? Default = False (no reweighting). No output written")
	parser.add_argument("-vf2o", "--vertexf2old", type=float, required=False, default=0.0, help="Value of old F2(0) from the gamma-f-fbar vertex Default = 0.0")
	parser.add_argument("-vf3o", "--vertexf3old", type=float, required=False, default=0.0, help="Value of old F3(0) from the gamma-f-fbar vertex. Default = 0.0")
	parser.add_argument("-vf2n", "--vertexf2new", type=float, required=False, default=0.0, help="Value of new F2(0) from the gamma-f-fbar vertex. Default = 0.0")
	parser.add_argument("-vf3n", "--vertexf3new", type=float, required=False, default=0.0, help="Value of new F3(0) from the gamma-f-fbar vertex. Default = 0.0")

	# Define the arguments
	args = parser.parse_args()

	# Print the banner
	print()
	print('========================================================')
	print('INFO:\tWelcome to LHE file converter + reweighting tool')
	print('========================================================')
	print()

	# print the starting time
	start_time = ctime()
	helpers.myprint('Starting the script at: %s'%(start_time), 'INFO', False, True)

	# Print the input arguments
	helpers.myprint('The following arguments have been provided:', 'INFO', False, False)
	for arg in vars(args): helpers.myprint('  --%s = %s'%(arg, getattr(args, arg)), 'INFO', False, False)
	print()

	# Some sanity checks on the input arguments
	if args.nfiles == 0:
		helpers.myprint('You set the number of input files to 0! Not doing anything...', 'ERROR', True, False)
		helpers.myprint('Aborting...', 'ERROR', False, True)
		return()	
	if args.format not in [1,2]:
		helpers.myprint('Format for output LHE files not recognized! Please check --help!', 'ERROR', True, False)
		helpers.myprint('Aborting...', 'ERROR', False, True)
		return()
	if (not args.reweight) and (not args.output):
		helpers.myprint('Only converting file formats, but no output given', 'ERROR', True, False)
		helpers.myprint('Aborting...', 'ERROR', False, True)
		return()
	if (not args.output) and (args.reweight):
		helpers.myprint('Reweighting requested but no output name given', 'ERROR', True, False)
		helpers.myprint('Aborting...', 'ERROR', False, True)
		return()

	# Get all the input LHE filenames
	bad_input_flag, file_names = helpers.get_file_names(args)
	if bad_input_flag == 1:
		helpers.myprint('Not processing anything. Please check the path to the files!', 'ERROR', True, False)
		helpers.myprint('Aborting...', 'ERROR', False, True)
		return()

	# Print the names of the input files just to be sure
	helpers.myprint('The following LHE files will be read for conversion + reweighting:', 'INFO', True, False)
	for f_ in file_names: print(f_)
	print()

	# Start main event processing loop
	helpers.myprint('Starting main event processing loop!', 'INFO', True, True)
	helpers.event_looper(args, file_names)
	helpers.myprint('Finished processing all events! Have a nice day!', 'INFO', True, False)

	# print the end time
	end_time = ctime()
	helpers.myprint('Time: %s'%(end_time), 'INFO', False, True)

if __name__ == "__main__": 
	main()