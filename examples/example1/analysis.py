#!/usr/bin/env python
# ---------------------------------------------------------------------
#  File:        analysis.py
#  Description: Analyze the results of RGS and find the best cuts.
#               Definitions:
#                 1. A cut is a threshold on a single variable.
#                    e.g., x > xcut
#                 2. A cut-point is the AND of a sequence of cuts. This
#                    can be visualized as a point in the space of cuts.
#                 3. A box cut is a two-sided threshold.
#                    e.g., (x > xlow) and (x < xhigh)
#                 4. A ladder cut is the OR of cut-pooints.
# ---------------------------------------------------------------------
#  Created:     10-Jan-2015 Harrison B. Prosper and Sezen Sekmen
# ---------------------------------------------------------------------
import os, sys, re
from string import *
from histutil import *
from time import sleep
from array import array
from ROOT import *
# ---------------------------------------------------------------------
def error(message):
    print "** %s" % message
    exit(0)
# ---------------------------------------------------------------------
def plotData():

    msize = 0.30 # marker size
    
    xbins =   50
    xmin  =  0.0
    xmax  = 10.0

    ybins =    50
    ymin  =   0.0
    ymax  =5000.0    

    cmass = TCanvas("fig_example1_VBF_ggF", "VBF/ggF",
                    10, 10, 500, 500)    
    
    # -- background
    hb = mkhist2("hb",
                 "#Delta#font[12]{#eta_{jj}}",
                 "#font[12]{m_{jj}} (GeV)",
                 xbins, xmin, xmax,
                 ybins, ymin, ymax,
                 color=kMagenta+1)
    hb.Sumw2()
    hb.SetMarkerSize(msize)
    hb.GetYaxis().SetTitleOffset(1.80)
    
    bntuple = Ntuple('../data/ggf13TeV_test.root', 'Analysis')
    btotal  = 0.0
    total  = 0
    for ii, event in enumerate(bntuple):
        btotal += event.weight
        total  += 1
        hb.Fill(event.deltaetajj, event.massjj, event.weight)
        if total % 100 == 0:
            cmass.cd()
            hb.Draw('p')
            cmass.Update()
    
    # -- signal
    hs = mkhist2("hs",
                 "#Delta#font[12]{#eta_{jj}}",
                 "#font[12]{m_{jj}} (GeV)",
                 xbins, xmin, xmax,
                 ybins, ymin, ymax,
                 color=kCyan+1)
    hs.Sumw2()
    hs.SetMarkerSize(msize)
    hs.GetYaxis().SetTitleOffset(1.80)
    
    sntuple = Ntuple('../data/vbf13TeV_test.root', 'Analysis')
    stotal  = 0.0
    total   = 0
    for event in sntuple:
        stotal += event.weight
        total  += 1
        hs.Fill(event.deltaetajj, event.massjj, event.weight)
        if total % 100 == 0:
            cmass.cd()
            hs.Draw('p')
            cmass.Update()

    cmass.cd()
    hs.Draw('p')
    hb.Draw('p same')
    cmass.Update()
    return (cmass, hs, hb)
# ---------------------------------------------------------------------
def main():
    print "="*80
    print "\t=== Example 1 - Analyze results of Box Cuts ==="
    print "="*80

    resultsfilename = "example1.root"
    treename = "RGS"
    print "\n\topen RGS file: %s"  % resultsfilename
    ntuple = Ntuple(resultsfilename, treename)
    
    variables = ntuple.variables()
    for name, count in variables:
        print "\t\t%-30s\t%5d" % (name, count)        
    print "\tnumber of cut-points: ", ntuple.size()

    # -------------------------------------------------------------
    # Plot results of RGS, that is, the fraction of events that
    # pass a given cut-point.
    #  1. Loop over cut points and compute a significance measure
    #     for each cut-point.
    #  2. Find cut-point with highest significance.
    # -------------------------------------------------------------
    # Set up a standard Root graphics style (see histutil.py in the
    # python directory).
    setStyle()

    # Create a 2-D histogram for ROC plot
    msize = 0.30  # marker size for points in ROC plot
    
    xbins =   50  # number of bins in x (background)
    xmin  =  0.0  # lower bound of x
    xmax  =  1.0  # upper bound of y

    ybins =   50
    ymin  =  0.0
    ymax  =  1.0

    color = kBlue+1
    hist  = mkhist2("hroc",
                    "#font[12]{#epsilon_{B}}",
                    "#font[12]{#epsilon_{S}}",
                    xbins, xmin, xmax,
                    ybins, ymin, ymax,
                    color=color)
    hist.SetMinimum(0)
    hist.SetMarkerSize(msize)


    # loop over all cut-points, compute a significance measure Z
    # for each cut-point, and find the cut-point with the highest
    # significance and the associated cuts.
    print "\tfilling ROC plot..."	
    bestZ = -1      # best Z value
    bestRow = -1    # row with best cut-point

    for row, cuts in enumerate(ntuple):
        fb = cuts.fraction_b  #  background fraction
        fs = cuts.fraction_s  #  signal fraction
        b  = cuts.count_b     #  background count
        s  = cuts.count_s     #  signal count
                
        #  Plot fs vs fb
        hist.Fill(fb, fs)
        	
        # Compute measure of significance
        #   Z  = sign(LR) * sqrt(2*|LR|)
        # where LR = log(Poisson(s+b|s+b)/Poisson(s+b|b))
        Z = 0.0
        if b > 1:
            Z = 2*((s+b)*log((s+b)/b)-s)
            absZ = abs(Z)
            if absZ != 0:
                Z = Z*sqrt(absZ)/absZ                    
        if Z > bestZ:
            bestZ = Z
            bestrow = row

    # -------------------------------------------------------------            
    # Write out best cut
    # -------------------------------------------------------------
    ntuple.read(bestrow)
    print "\nBest cuts"
    bestcuts = {}
    for name, count in variables:    
        if name[0:5] in ['count', 'fract', 'cutpo']: continue
        var = ntuple.get(name)
        bestcuts[name] = var
        print "\t%s" % name
        for ii in xrange(len(var)):
            print "\t\t%10.2f" % var[ii]
    print
    
    print "Yields and relative efficiencies"
    for name, count in variables:        
        if not (name[0:5] in ['count', 'fract']): continue
        var = ntuple.get(name)
        print "\t%-30s %10.3f" % (name, var)
        if name[0:5] == "fract":
            print
    
    # -------------------------------------------------------------
    # Save plots
    # -------------------------------------------------------------
    print "\t== plot ROC ==="	
    croc = TCanvas("fig_%s_ROC" % nameonly(resultsfilename),
                   "ROC", 520, 10, 500, 500)
    croc.cd()
    hist.Draw()
    croc.Update()
    croc.SaveAs(".png")


    print "\t=== plot cuts ==="

    cmass, hs, hb = plotData()
    
    xbins= hs.GetNbinsX()
    xmin = hs.GetXaxis().GetBinLowEdge(1)
    xmax = hs.GetXaxis().GetBinUpEdge(xbins)

    ybins= hs.GetNbinsY()
    ymin = hs.GetYaxis().GetBinLowEdge(1)
    ymax = hs.GetYaxis().GetBinUpEdge(ybins)
    
    hcut = TH2Poly('hcut', '', xmin, xmax, ymin, ymax)
    hcut.AddBin(bestcuts['deltaetajj'][0], bestcuts['massjj'][0],
                bestcuts['deltaetajj'][1], bestcuts['massjj'][1])
    cmass.cd()
    hcut.Draw('same')
    cmass.SaveAs('.png')
    
    sleep(5)
# ---------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print "bye!"


