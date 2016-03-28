#!/usr/bin/python

import argparse
import sys
import yaml
import BumpHuntToyAnalysis

def parse_config(config_file) :

    print "Loading configuration from " + str(config_file)
    config = open(config_file, 'r')
    return yaml.load(config)

def main() : 
   
    # Parse all command line arguments using the argparse module
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-c", "--config", 
                        help="Configuration file describing analysis paramters.")
    args = parser.parse_args()

    if args.config is None: 
        print 'A configuration file needs to be specified!'
        sys.exit(2)
    
    config = parse_config(args.config)
   
    bh = BumpHuntToyAnalysis.BumpHuntToyAnalysis()

    for input_file in config["Files"] : 
        print str(input_file)
        bh.process(input_file)

    bh.make_plots()

if __name__ == "__main__":
    main()
