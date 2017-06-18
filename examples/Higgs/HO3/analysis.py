#!/usr/bin/env python
# ---------------------------------------------------------------------
#  File:        analysis.py
#  Description: Analysis of the results of RGS and find the best cuts.
#               Definitions:
#                 1. A one-sided cut is a threshold on a single variable.
#                    e.g., x > xcut
#                 2. A cut-point is the AND of a sequence of one=sided cuts.
#                    A cut-point can be visualized as a point in the
#                    space of cuts.
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
NAME = 'HO3'
# ---------------------------------------------------------------------
def cut(event):
    skip = \
      (event.massjj <=     0) or \
      (event.mass4l <=   100) or \
      (event.mass4l >=   150) or \
      (event.Z1mass <= 58.35) or \
      (event.Z1mass >= 94.58) or \
      (event.Z2mass <= 17.00) or \
      (event.Z2mass >= 52.11)
    return skip
# ---------------------------------------------------------------------
def fill():

    xbins =   50
    xmin  =  0.0
    xmax  = 10.0

    ybins =    50
    ymin  =   0.0
    ymax  =2000.0    
    
    cmass = TCanvas("h_%s" % NAME, "VBF/ggF/ZZ",
                    10, 10, 425, 400)    
    
    # -- background
    hb = mkhist2("hb",
                 "|#Delta#font[12]{#eta_{jj}}|",
                 "#font[12]{m_{jj}} (GeV)",
                 xbins, xmin, xmax,
                 ybins, ymin, ymax,
                 color=kMagenta+1)
    hb.GetYaxis().SetTitleOffset(1.80)
    
    bntuple = Ntuple(['../../data/ntuple_HZZ4L.root',
                      '../../data/ntuple_ZZ4L.root'],
                      'Analysis')
    for ii, event in enumerate(bntuple):
        if cut(event): continue
        hb.Fill(event.detajj, event.massjj, event.weight)
        if ii % 1000 == 0:
            cmass.cd()
            hb.Draw('box')
            cmass.Update()
            gSystem.ProcessEvents()
            
    # -- signal
    hs = mkhist2("hs",
                 "|#Delta#font[12]{#eta_{jj}}|",
                 "#font[12]{m_{jj}} (GeV)",
                 xbins, xmin, xmax,
                 ybins, ymin, ymax,
                 color=kAzure+1)
    hs.GetYaxis().SetTitleOffset(1.80)
    
    sntuple = Ntuple('../../data/ntuple_HZZ4L_VBF.root', 'Analysis')
    for event in sntuple:
        if cut(event): continue
        hs.Fill(event.detajj, event.massjj, event.weight)
        if ii % 1000 == 0:
            cmass.cd()
            hs.Draw('box')
            cmass.Update()
            gSystem.ProcessEvents()
            
    cmass.cd()
    hs.Scale(1.0/hs.Integral())
    hb.Scale(1.0/hb.Integral())
    hs.Draw('box')
    hb.Draw('box same')
    cmass.Update()
    gSystem.ProcessEvents()    
    return (cmass, hs, hb)
# ---------------------------------------------------------------------
def main():
    print "="*80
    print "\t=== HO3 find best cut ==="
    print "="*80

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

    cmass, hs, hb = fill()
    
    # Create a 2-D histogram for ROC plot
    msize = 0.30  # marker size for points in ROC plot
    
    xbins =   50   # number of bins in x (background)
    xmin  =  0.0   # lower bound of x
    xmax  =  0.02  # upper bound of y

    ybins =   50
    ymin  =  0.0
    ymax  =  0.2

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
    bestZ   = -1    # best Z value
    bestRow = -1    # row with best cut-point

    totals = ntuple.totals()
    
    ts, es  = totals[0]
    t1, et1 = totals[1]
    t2, et2 = totals[2]
    tb = t1 + t2
    for row, cuts in enumerate(ntuple):
        s  = cuts.count_VBF #  signal count
        b1 = cuts.count_ggF #  background 1
        b2 = cuts.count_ZZ  #  background 2
        b  = b1 + b2
        
        fs = s / ts
        fb = b / tb
                                    
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
                                   '%s.cuts' % NAME,
                                ntuple, variables,
                                bestrow, bestZ,
                                totals)
        
    # -------------------------------------------------------------
    # Make plots
    # -------------------------------------------------------------
    print "\n\t== plot ROC ==="	
    croc = TCanvas("h_%s_ROC" % NAME,
                   "ROC", 520, 10, 425, 400)
    croc.cd()
    hist.Draw()
    croc.Update()
    gSystem.ProcessEvents()    
    croc.SaveAs(".pdf")    

    print "\n\t=== plot one-sided cuts ==="
    xbins= hs.GetNbinsX()
    xmin = hs.GetXaxis().GetBinLowEdge(1)
    xmax = hs.GetXaxis().GetBinUpEdge(xbins)

    ybins= hs.GetNbinsY()
    ymin = hs.GetYaxis().GetBinLowEdge(1)
    ymax = hs.GetYaxis().GetBinUpEdge(ybins)

    xcut = array('d')
    xcut.append(bestcuts['detajj'])
    xcut.append(bestcuts['detajj'])
    xcut.append(xmax)

    ycut = array('d')
    ycut.append(ymax)
    ycut.append(bestcuts['massjj'])
    ycut.append(bestcuts['massjj'])
    hcut = TGraph(3, xcut, ycut)

    cmass.cd()
    hs.Draw('box')
    hb.Draw('boxsame')
    hcut.Draw('same')
    drawHZZLegend()
    cmass.Update()
    gSystem.ProcessEvents()    
    cmass.SaveAs('.pdf')

    print "\n\t=== plot approximation to optimal discriminant ==="
    # approximate D = p(x|S)/[p(x|S)+p(x|B)]
    # using histograms
    hmvd = hs.Clone('hmvd')    
    hsum = hs.Clone('hsum')
    hsum.Add(hb)
    hmvd.Divide(hsum)

    cmvd = TCanvas("h_%s_D" % NAME, "D(x, y)", 1040, 10, 425, 400)
    cmvd.cd()
    cmvd.SetRightMargin(0.14)
    hmvd.Draw('cont4z')
    drawHZZLegend(left=False, postfix='2')
    cmvd.Update()
    gSystem.ProcessEvents()
    cmvd.SaveAs('.pdf')
    
    sleep(5)
# ---------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print "bye!"


