#!/usr/bin/env python
# -----------------------------------------------------------------------------
#  File:        train.py
#  Description: Example of Random Grid Search to find the results of an
#               ensemble cuts. 
#  Created:     10-Jan-2015 Harrison B. Prosper and Sezen Sekmen
# -----------------------------------------------------------------------------
import os, sys, re
from rgsutil import *
from string import *
from ROOT import *
# -----------------------------------------------------------------------------
def main():
    print "="*80
    print "\t=== Example 1 - Box Cuts ==="
    print "="*80

    # ---------------------------------------------------------------------
    # Load the RGS shared library and check that the various input files
    # exist.
    # ---------------------------------------------------------------------
    if gSystem.Load("libRGS") < 0:
        error("unable to load libRGS")

    # Name of file containing cut definitions
    # Format of file:
    #   variable-name  cut-type (>, <, <>, |>, |<, ==)
    varfilename = "example1.cuts"
    if not os.path.exists(varfilename):
        error("unable to open variables file %s" % varfilename)

    # Name of signal file
    sigfilename = "../data/vbf13TeV_train.root"
    if not os.path.exists(sigfilename):
        error("unable to open signal file %s" % sigfilename)

    # Name of background file        
    bkgfilename = "../data/ggf13TeV_train.root"
    if not os.path.exists(bkgfilename):
        error("unable to open background file %s" % bkgfilename)

    # ---------------------------------------------------------------------
    #  Create RGS object
    #  
    #   The file (cutdatafilename) of cut-points is usually a signal file,
    #   which ideally differs from the signal file on which the RGS
    #   algorithm is run.
    # ---------------------------------------------------------------------
    cutdatafilename = sigfilename
    start      = 0           # start row 
    maxcuts    = 1000        # maximum number of cut-points to consider
    treename   = "Analysis"  # name of Root tree 
    weightname = "weight"    # name of event weight variable
    
    rgs = RGS(cutdatafilename, start, maxcuts, treename, weightname)

    # ---------------------------------------------------------------------
    #  Add signal and background data to RGS object.
    #  Weight each event using the value in the field weightname, if
    #  present.
    #  NB: We asssume all files are of the same format.
    # ---------------------------------------------------------------------
    start    = 0 #  start row
    numrows  =-1 #  scan all the data from the files
    # The last (optional) argument is a string, which, if given, will be
    # appended to the "count" and "fraction" variables. The "count" variable
    # contains the number of events that pass per cut-point, while "fraction"
    # is count / total, where total is the total number of events per file.
    # If no string is given, the default is to append an integer to the
    # "count" and "fraction" variables, starting at 0, in the order in which
    # the files are added to the RGS object.
    rgs.add(bkgfilename, start, numrows, "_b")
    rgs.add(sigfilename, start, numrows, "_s")

    # ---------------------------------------------------------------------	
    #  Run RGS and write out results
    # ---------------------------------------------------------------------	    
    rgs.run(varfilename)

    # Write to a root file
    rgsfilename = "%s.root" % nameonly(varfilename)
    rgs.save(rgsfilename)

    # Write to a text file
    rgsfilename = "%s.txt" % nameonly(varfilename)
    rgs.save(rgsfilename)    
# -----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print "\tciao!\n"



