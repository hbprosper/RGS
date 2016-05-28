#!/usr/bin/env python
# ---------------------------------------------------------------------
#  File:        analysisFromTxt.py
#  Description: Analyze the results of RGS ladder cuts and find the
#               best cuts.
# ---------------------------------------------------------------------
#  Created:     10-Jan-2015 Harrison B. Prosper and Sezen Sekmen
# ---------------------------------------------------------------------
import os, sys, re
from string import *
from histutil import *
from time import sleep
from LadderCut import LadderCut
from ROOT import *
# ---------------------------------------------------------------------
def error(message):
    print "** %s" % message
    exit(0)
# ---------------------------------------------------------------------
def main():
    print "="*80
    print "\t=== Example 2 - Analyze results of Ladder Cuts ==="    
    print "="*80

    setStyle()

    msize = 0.30
    xbins = 40
    xmin  = 0.0
    xmax  =4000.0

    ybins = 40
    ymin  = 0.0
    ymax  = 1.0    

    cmass = TCanvas("fig_T2tt_TTJets", "", 10, 10, 700, 350)    
    # divide canvas canvas along x-axis
    cmass.Divide(2,1)
    
    # -- background

    hb = mkhist2("hb",
                 "M_{R} (GeV)",
                 "R^{2}",
                 xbins, xmin, xmax,
                 ybins, ymin, ymax,
                 color=kMagenta+1)
    hb.Sumw2()
    hb.SetMarkerSize(msize)
    bntuple = Ntuple('../data/TTJets.root', 'Analysis')
    for ii, event in enumerate(bntuple):
        hb.Fill(event.MR, event.R2)
        if ii % 100 == 0:
            cmass.cd(2)
            hb.Draw('p')
            cmass.Update()
    
    # -- signal
    hs = mkhist2("hs",
                 "M_{R} (GeV)",
                 "R^{2}",
                 xbins, xmin, xmax,
                 ybins, ymin, ymax,
                 color=kCyan+1)
    hs.Sumw2()
    hs.SetMarkerSize(msize)
    sntuple = Ntuple('../data/T2tt_mStop_850_mLSP_100.root',
                     'Analysis')
    for ii, event in enumerate(sntuple):
        hs.Fill(event.MR, event.R2)        
        if ii % 100 == 0:
            cmass.cd(2)
            hs.Draw('p')
            cmass.Update()

    # compute D = p(x|S)/[p(x|S)+p(x|B)]
    hD = hs.Clone('hD'); hD.Scale(1.0/hD.Integral())
    hB = hb.Clone('hB'); hB.Scale(1.0/hB.Integral())
    
    hSum = hD.Clone('hSum')
    hSum.Add(hB)
    hD.Divide(hSum)

    cmass.cd(1)
    hD.Draw('cont')
    
    cmass.cd(2)
    hs.Draw('p')
    hb.Draw('p same')
    cmass.Update()
    
    # -------------------------------------------------------------
    #  Plot results of RGS
    # -------------------------------------------------------------

    resultsfilename = "example2.txt"
    print "\n\topen RGS file: %s"  % resultsfilename
    table = Table(resultsfilename)
    variables = table.variables()
    for name, count in variables:
        print "\t\t%-30s\t%d" % (name, count)
    totalcutpoints = len(table)
    print "\tnumber of cut-points: %d" % totalcutpoints
    
    bmax = 0.30
    smax = 1.00
    color= kBlue+1
    hist = mkhist2("hroc",
                   "#font[12]{#epsilon_{t#bar{t}}}",
                   "#font[12]{#epsilon_{T2tt}}",
                   xbins, 0.0, 0.30,
                   ybins, 0.7, 1.00,
                   color=color)
    hist.SetMinimum(0)

    laddercut = LadderCut(xmin, xmax, ymin, ymax)
    print "\tfilling ROC plot..."	
    for row, cuts in enumerate(table):
        fb = cuts("fraction_b")  #  background fraction
        fs = cuts("fraction_s")  #  signal fraction 
        b  = cuts("count_b")     #  background count
        s  = cuts("count_s")     #  signal count
   
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
 
        # add ladder cut to ladder plot
        R2 = cuts('R2')
        MR = cuts('MR')
        laddercut.add(Z, R2, MR)

    print "\n\t=== best ladder cut"
    signif, outerhull, cutpoints = laddercut(0)
    delim = ''
    for ii, cutpoint in enumerate(outerhull):
        print "\t%4s\t(R2 > %8.4f) AND (MR > %8.1f)" % (delim,
                                                        cutpoint[0],
                                                        cutpoint[1])
        delim = 'OR'
        
    print "\n\t=== plot ladder cut"
    cmass.cd(2)    	
    laddercut.draw()
    cmass.Update()
    cmass.SaveAs('.png')

    # -------------------------------------------------------------
    # roc plot
    # -------------------------------------------------------------
    print "\n\t=== plot ROC ==="	
    croc = TCanvas("fig_example2_ROC", "RGS", 516, 10, 500, 500)
    croc.cd()
    hist.Draw()
    croc.Update()
    croc.SaveAs(".png")    
    sleep(5)
# ---------------------------------------------------------------------
main()



