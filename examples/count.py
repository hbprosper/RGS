#!/usr/bin/env python
# ----------------------------------------------------------------------------
#  File:        count.py
#  Description: sum weights of ntuples
#  Created: ??-Oct-2016 HBP
# ----------------------------------------------------------------------------
import os, sys, re
from ROOT import *
from string import *
from rgsutil import *
from time import sleep
# ----------------------------------------------------------------------------
VARS = '''
    double weight
    int njet
    int nb
    int nW
    double j1pT
    double MR
    double R2
'''
VARS = map(split, split(strip(VARS), '\n'))
TREENAME = 'RGSinput'
# ----------------------------------------------------------------------------
def main():
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        sys.exit('''
    Usage:
        python count.py root-file

    Example:
        python count.py data/SUSY_gluino_m1354.root
        ''')    
    filenames = [filename]
        
    # do some PyRoot magic to read from tree
    struct = 'struct Event{'
    for vtype, vname in VARS:
        struct += '%s %s;' % (vtype, vname)
    struct += '};'    
    gROOT.ProcessLine(struct)
    from ROOT import Event
    event = Event()
    
    # loop over files    
    for filename in filenames:
        print "\n%s" % filename
        
        # open root file
        ntuple = TFile(filename)
        # get tree
        tree = ntuple.Get(TREENAME)

        # set addresses of fields
        for vtype, vname in VARS:
            tree.SetBranchAddress(vname, AddressOf(event, vname))
        
        # loop over events in current file
        nevents = tree.GetEntries()
        print "  entries:  %10d" % nevents
        w1 = 0.0
        w2 = 0.0
        for index in xrange(nevents):
            tree.GetEntry(index)

            if not (event.njet >= 3):  continue
            if not (event.nb   >= 1):  continue
            if not (event.nW   >= 1):  continue
            if not (event.j1pT > 200): continue
            if not (event.MR > 800):   continue
            if not (event.R2 > 0.08):  continue                
                
            w = event.weight
            w1 += w
            w2 += w*w
    w2 = sqrt(w2)
    print '  count:\t%10.3f +/- %-10.3f\n' % (w1, w2)
# -------------------------------------------------------------------------
try:     
    main()
except KeyboardInterrupt:
    print '\nciao!\n'
