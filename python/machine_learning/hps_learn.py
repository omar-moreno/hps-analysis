#!/usr/bin/python

###############
#   Imports   #
###############

import argparse
import ROOT as root
import ConfigReader

############
#   Main   #
############

def main():

    # Parse all command line arguments using the argparse module
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-c", "--config", help="Configuration file")
    args = parser.parse_args()

    config_reader = ConfigReader(args.config)

if __name__ == "__main__":
    main()
