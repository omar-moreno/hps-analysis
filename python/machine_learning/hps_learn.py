#!/usr/bin/python

###############
#   Imports   #
###############

import argparse
import sys
import ROOT as root
from ConfigReader import ConfigReader

############
#   Main   #
############

def main():

    # Parse all command line arguments using the argparse module
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-c", "--config", help="Configuration file")
    args = parser.parse_args()

    if args.config is None: 
        print '[ hps-learn ]: A configuration file needs to be specified!'
        sys.exit(2)

    config_reader = ConfigReader(args.config)

    signal_files = config_reader.get_signal_files()
    print "Signal files: " + str(signal_files)
    

if __name__ == "__main__":
    main()
