GeneralLHCbAnalysis
===================

A set of scripts for generically producing an LHCb analysis

For example you can use an input file (or a set of input files) in the ROOT TTree format
to generate a list of branches with ./dumpTreeBranches.py (run with -h or --help to see options)

Then you can use these list of branches to create a new package with a few default classes
which will setup an analysis from these tree with ./setupNewAnalysis.py (run with -h or --help
to see options)

The basic package structure will have src and interface directories for your classes.
You can make an analysis class to select, reject or shuffle events which inherits from BaseAnalyser.cc

Their is a makefile to compile it. It will then by default be run by a script called runAnalysis.py
