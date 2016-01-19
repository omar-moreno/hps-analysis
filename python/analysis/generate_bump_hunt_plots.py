#!/usr/bin/python

import argparse
import sys
from rootpy.io import root_open
from rootpy.tree import Tree
import rootpy.plotting.root2matplotlib as rplt

def main() : 
    
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input", help="ROOT file")
    args = parser.parse_args()

    if args.input is None : 
        print 'A ROOT file needs to be specified'
        sys.exit(2)

    root_file = root_open(args.input)
    
    tree = root_file.results

    tree.Draw("yield")

if __name__ == "__main__" : 
    main()

