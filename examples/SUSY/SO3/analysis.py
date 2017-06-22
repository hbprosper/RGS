#!/usr/bin/env python
# ---------------------------------------------------------------------
#  File:        analysis.py
#  Description: SO1: Analyze the results of RGS staircase cuts and find
#               the best cuts.
# ---------------------------------------------------------------------
#  Created:     10-Jan-2015 Harrison B. Prosper and Sezen Sekmen
#               15-Oct-2016 HBP now refer to staircase cuts
#               17-Jun-2017 HBP adapt to latest version of OuterHull
# ---------------------------------------------------------------------
import os, sys, re
from string import *
from rgsutil import *
from time import sleep
from ROOT import *

sys.path.append('../../python')
from rgsexamples import *
# ---------------------------------------------------------------------
NAME = 'SO3'
def cut(event):
    skip = \
      (event.njet <     3)
    return skip
# ---------------------------------------------------------------------
def main():
    global cut
    
    print "="*80
    print "\t=== %s: find best cuts ===" % NAME
    print "="*80
    
    treename = "RGS"
    varfilename  = "%s.cuts" % NAME
    resultsfilename= "%s.root" % NAME

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

    cmass, hs, hb = fill_MR_R2(cut, 'SO3')
    
    # Create a 2-D histogram for ROC plot
    msize = 0.30  # marker size for points in ROC plot
    
    xbins =  10000   # number of bins in x (background)
    xmin  =  0.0    # lower bound of x
    xmax  =  1.0    # upper bound of y

    ybins =  50
    ymin  =  0.0
    ymax  =  1.0

    color = kBlue+1
    hroc  = mkhist2("hroc",
                    "#font[12]{#epsilon_{B}}",
                    "#font[12]{#epsilon_{S}}",
                    xbins, xmin, xmax,
                    ybins, ymin, ymax,
                    color=color)
    hroc.SetMinimum(0)
    hroc.SetMarkerSize(msize)

    hs.Scale(1.0/hs.Integral())
    hb.Scale(1.0/hb.Integral())
    
    cmass.cd()
    hb.Draw('box')
    hs.Draw('box same')
    cmass.Update()
    gSystem.ProcessEvents()
    
    # initialize an object to determine the outer hull of the
    # staircase cuts, that is, the minimum set of cuts that comprise
    # the staircase cut
    cutdirs = []
    for t in getCutDirections(varfilename):
        token = t[0]
        if token == '\\staircase':
            continue
        elif token == '\\end':
            continue
        else:
            cutdirs.append(t)
            
    xbins = hs.GetNbinsX()
    xmin  = hs.GetXaxis().GetBinLowEdge(1)
    xmax  = hs.GetXaxis().GetBinLowEdge(xbins)+hs.GetXaxis().GetBinWidth(xbins)

    ybins = hs.GetNbinsY()
    ymin  = hs.GetYaxis().GetBinLowEdge(1)
    ymax  = hs.GetYaxis().GetBinLowEdge(ybins)+hs.GetYaxis().GetBinWidth(ybins)

    outerHull = OuterHull(xmin, xmax, ymin, ymax, cutdirs)

    # loop over all cut-points, compute a significance measure Z
    # for each cut-point, and find the cut-point with the highest
    # significance and the associated cuts.

    bestZ   = -1    # best Z value
    bestRow = -1    # row with best cut-point

    for row, cuts in enumerate(ntuple):
        b  = cuts.count_b #  background count        
        s  = cuts.count_s #  signal count

        fb = cuts.fraction_b
        fs = cuts.fraction_s        
        
        hroc.Fill(fb, fs)

        Z = signalSignificance(s, b)
        if Z > bestZ:
            bestZ = Z
            bestRow = row
            
        # add staircase cut to outer hull object
        outerHull.add(Z, cuts.MR, cuts.R2)
        
    Z, outerhull, cutpoints = outerHull(0)
    # -------------------------------------------------------------            
    # Write out results
    # -------------------------------------------------------------
    bestcuts = writeSUSYResults('r_%s.txt' % NAME,
                                '%s.cuts' % NAME,
                                ntuple, variables,
                                bestRow, bestZ,
                                outerhull)

    # -------------------------------------------------------------
    # plot
    # -------------------------------------------------------------
  
    print "\n\t=== plot outer hull"
    color = kGreen+3
    cut = outerHull(0)
        
    cmass.cd()
    hs.Draw('box')
    hb.Draw('boxsame')    
    outerHull.draw(cut,
                   hullcolor=kGreen+3,
                   lstyle=1,
                   lwidth=3)
    
    drawSUSYLegend('with preselection', left=False)
    cmass.Update()
    gSystem.ProcessEvents()    
    cmass.SaveAs('.pdf')
    
    print "\t== plot ROC ==="	
    croc = TCanvas("h_%s_ROC" % NAME, "ROC", 520, 10, 500, 500)
    croc.cd()
    croc.SetLogx()
    hroc.Draw()
    croc.Update()
    gSystem.ProcessEvents()    
    croc.SaveAs(".pdf")  

    sleep(5)
# ---------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print '\nciao!'


