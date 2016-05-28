#!/usr/bin/env python
# -----------------------------------------------------------------------------
#  File:        train.py
#  Description: Example of Random Grid Search to find cuts that best separates
#               T2tt from TTJets using ladder cuts
#  Created:     10-Jan-2015 Harrison B. Prosper
# -----------------------------------------------------------------------------
import os, sys, re
from string import *
from ROOT import *
# -----------------------------------------------------------------------------
def error(message):
    print "** %s" % message
    exit(0)
def nameonly(s):
    import posixpath
    return posixpath.splitext(posixpath.split(s)[1])[0]    
# -----------------------------------------------------------------------------
def main():
    print "="*80
    print "\t=== Example 2 - Ladder Cuts ==="
    print "="*80
    
    # ---------------------------------------------------------------------
    # Load the RGS shared library and check that the various input files
    # exist.
    # ---------------------------------------------------------------------
    if gSystem.Load("libRGS") < 0:
        error("unable to load libRGS")

    # Name of file containing cut definitions
    # Format of file for ladder cuts:
    #   \ladder number of cut-points
    #      variable-name cut-type (>, <, |>, |<, ==)
    #           :           :
    #   \end        
    varfilename = "example2.cuts"
    if not os.path.exists(varfilename):
        error("unable to open variables file %s" % varfilename)

    sigfilename = "../data/T2tt_mStop_850_mLSP_100.root"
    if not os.path.exists(sigfilename):
        error("unable to open signal file %s" % sigfilename)

    bkgfilename = "../data/TTJets.root"
    if not os.path.exists(bkgfilename):
        error("unable to open background file %s" % bkgfilename)

    # ---------------------------------------------------------------------
    #  Create RGS object
    #  Need:
    #   A file of cut-points - usually a signal file, which ideally is
    #   not the same as the signal file on which the RGS algorithm is run.
    # ---------------------------------------------------------------------
    print "==> create RGS object"
    cutfilename = sigfilename
    start   = 0    
    maxcuts = 2000 #  maximum number of cut-points to consider
    treename= "Analysis"
    rgs = RGS(cutfilename, start, maxcuts, treename)

    # ---------------------------------------------------------------------
    #  Add signal and background data to RGS object
    #  Weight each event using the value in the field weightname, if it
    #  exists.
    #  NB: We asssume all files are of the same format
    # ---------------------------------------------------------------------
    start   = 0
    numrows =-1 #  Load all the data from the files
    rgs.add(bkgfilename, start, numrows, '_b')
    rgs.add(sigfilename, start, numrows, '_s')

    # ---------------------------------------------------------------------	
    #  Run RGS and write out results
    # ---------------------------------------------------------------------	    
    rgs.run(varfilename)
    
    rgsfilename = "%s.root" % nameonly(varfilename)
    rgs.save(rgsfilename)
    
    rgsfilename = "%s.txt" % nameonly(varfilename)
    rgs.save(rgsfilename)    
# -----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print "\tciao!\n"



