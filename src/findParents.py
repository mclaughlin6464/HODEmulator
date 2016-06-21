#!/bin/bash
#@Author: Sean McLaughlin
#This module runs rockstar's findParents on a collection of outlists and writes them to a new location.

import argparse
from subprocess import call
from os.path import isdir


desc = "Run ROCKSTAR's find_parents on a collcetion of outlists and write them to a new location. ."
parser = argparse.ArgumentParser(description = desc)
parser.add_argument('dirname', metavar='dirname', type=str,
                    help='The name of the directory to collect the outlists from')
parser.add_argument('outputname', metavar='outputname', type=str,
                    help='The name of the directory to store the finished outlists')
parser.add_argument('--box_size', type = float, help = 'Size of the Box. Optional')
args = parser.parse_args()

for input in (args.dirname, args.outputname):
    assert isdir(input)
