#!/usr/bin/env python
# ---------------------------------------------------------------------
#  File:        analysis.py
#  Description: HO1: Analyze the results of RGS and find the best cuts.
#               Definitions:
#                 1. A one-sided cut is a threshold on a single
#                    variable.
#                    e.g., x > xcut
#                 2. A cut-point is the AND of a sequence of cuts. This
#                    can be visualized as a point in the space of cuts.
#                 3. A two-sided cut is a two-sided threshold.
#                    e.g., (x > xlow) and (x < xhigh)
#                 4. A staircase cut is the OR of cut-points.
# ---------------------------------------------------------------------
#  Created:     10-Jan-2015 Harrison B. Prosper and Sezen Sekmen
# ---------------------------------------------------------------------
import os, sys, re
from string import *
from rgsutil import *
from time import sleep
from ROOT import *

sys.path.append('../../python')
from rgsexamples import *

# ---------------------------------------------------------------------
def cut(event):
    return False
# ---------------------------------------------------------------------
def main():
    print "="*80
    print "\t=== HO1: obtain best uni-directional cuts ==="
    print "="*80

    NAME = 'HO1'
    
    resultsfilename = "%s.root" % NAME

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

    cmass, hs, hb = fillZ1massZ2mass(cut)
    
    # Create a 2-D histogram for ROC plot
    msize = 0.30  # marker size for points in ROC plot
    
    xbins =   50  # number of bins in x (background)
    xmin  =  0.0  # lower bound of x
    xmax  =  0.1  # upper bound of y

    ybins =   50
    ymin  =  0.0
    ymax  =  0.3

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

    totals = ntuple.totals()

    t_VBF, et1 = totals[0]
    t_ggF, et2 = totals[1]
    t_ZZ,  et3 = totals[2]

    ts = t_VBF + t_ggF
    tb = t_ZZ

    for row, cuts in enumerate(ntuple):
        c_VBF = cuts.count_VBF
        c_ggF = cuts.count_ggF
        c_ZZ  = cuts.count_ZZ 
        s  = c_VBF + c_ggF
        b  = c_ZZ
        
        fs = s / ts
        fb = b / tb
                
        #  Plot fs vs fb
        hist.Fill(fb, fs)
        	
        # Compute measure of significance
        #   Z  = sign(LR) * sqrt(2*|LR|)
        # where LR = log(Poisson(s+b|s+b)/Poisson(s+b|b))
        Z = signalSignificance(s, b)
        if Z > bestZ:
            bestZ = Z
            bestrow = row

    # -------------------------------------------------------------            
    # Write out best cut
    # -------------------------------------------------------------
    bestcuts = writeHZZResults('r_%s.txt' % NAME,
                               '%s.cuts'  % NAME,
                               ntuple, variables,
                               bestrow, bestZ,
                               totals)
    
    # -------------------------------------------------------------
    # Save plots
    # -------------------------------------------------------------
    print "\t== plot ROC ==="	
    croc = TCanvas("h_%s_ROC" % NAME,
                   "ROC", 520, 10, 400, 400)
    croc.cd()
    hist.Draw()
    croc.Update()
    croc.SaveAs(".pdf")        
    gSystem.ProcessEvents()    


    print "\t=== plot uni-directional cuts ==="
    
    xbins= hs.GetNbinsX()
    xmin = hs.GetXaxis().GetBinLowEdge(1)
    xmax = hs.GetXaxis().GetBinUpEdge(xbins)

    ybins= hs.GetNbinsY()
    ymin = hs.GetYaxis().GetBinLowEdge(1)
    ymax = hs.GetYaxis().GetBinUpEdge(ybins)

    xcut = array('d')
    xcut.append(xmin)
    xcut.append(bestcuts['Z1mass'])
    xcut.append(bestcuts['Z1mass'])

    ycut = array('d')
    ycut.append(bestcuts['Z2mass'])
    ycut.append(bestcuts['Z2mass'])
    ycut.append(ymin)    
    hcut = TGraph(3, xcut, ycut)

    cmass.cd()
    hs.Draw('box')
    hb.Draw('boxsame')
    hcut.Draw('same')    
    drawHZZLegend()    
    cmass.Update()
    gSystem.ProcessEvents()    
    cmass.Print("h_%s.pdf" % NAME)
    
    sleep(5)
# ---------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print "bye!"


