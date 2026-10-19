Usage if you have exactly one file:

python3 clean.py --input /path/to/lhe/file.lhe --output /path/to/cleaned/lhe/file.lhe 

Usage if you have multiple LHE files, then the script will clean and merge then:

python3 clean.py --input /path/to/lhe/files/* --output /path/to/cleaned/and/merged/file.lhe 
