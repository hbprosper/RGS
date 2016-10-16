#!/usr/bin/env python
# ---------------------------------------------------------------------
#  File:        analysis.py
#  Description: Analyze the results of RGS staircase cuts and find the
#               best cuts.
# ---------------------------------------------------------------------
#  Created:     10-Jan-2015 Harrison B. Prosper and Sezen Sekmen
#               15-Oct-2016 HBP now refer to staircase cuts
# ---------------------------------------------------------------------
import os, sys, re
from string import *
from rgsutil import *
from time import sleep
from ROOT import *
# ---------------------------------------------------------------------
def main():
    print "="*80
    print "\t=== Example 2 - Analyze results of Staircase Cuts ==="    
    print "="*80

    setStyle()

    msize = 0.15
    xbins = 50
    xmin  = 0.0
    xmax  =2000.0

    ybins = 50
    ymin  = 0.0
    ymax  = 0.5    

    cmass = TCanvas("fig_example2", "", 10, 10, 1000, 500)    
    # divide canvas canvas along x-axis
    cmass.Divide(2, 1)
    
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
        hb.Fill(event('MR'), event('R2'))
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
        hs.Fill(event('MR'), event('R2'))        
        if ii % 100 == 0:
            cmass.cd(2)
            hs.Draw('p')
            cmass.Update()

    # approximate D = p(x|S)/[p(x|S)+p(x|B)]
    # using histograms
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

    # initialize an object to determine the outer hull of the
    # ladder cuts, that is, the minimum set of cuts that comprise
    # the ladder cut
    varfilename = 'example2.cuts'
    cutdirs     = getCutDirections(varfilename)[1:-1]
    outerHull   = OuterHull(xmin, xmax, ymin, ymax, cutdirs)
    
    # -------------------------------------------------------------
    #  Plot results of RGS
    # -------------------------------------------------------------
    resultsfilename = "example2.root"
    print "\n\topen RGS file: %s"  % resultsfilename    
    ntuple = Ntuple(resultsfilename, 'RGS')
    print "number of cut-points: ", ntuple.size()
    
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

    for row, cuts in enumerate(ntuple):
        fb = cuts('fraction_b')   #  background fraction
        fs = cuts('fraction_s')   #  signal fraction
        b  = cuts('count_b')      #  background count
        s  = cuts('count_s')      #  signal count
               
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
            else:
                Z = 0.0
 
        # add ladder cut to outer hull object
        outerHull.add(Z, cuts('MR'), cuts('R2'))
        
    print "\n\t=== plot outer hulls"
    xname, xdir = cutdirs[0]
    yname, ydir = cutdirs[1]

    for ii, color, lstyle in [(0,    kGreen+3, 1),
                              (1000, kBlue,    7),
                              (1500, kBlack,   9)]:

        # draw outer full of specified ladder cut
        cut = outerHull(ii)
        
        cmass.cd(1)
        outerHull.draw(cut, hullcolor=color, lstyle=lstyle,
                       lwidth=3,
                       plotall=False)

        cmass.cd(2)
        outerHull.draw(cut, hullcolor=color, lstyle=lstyle,
                       lwidth=3,
                       plotall=ii==0)        
        cmass.Update()
        
        # print out the cuts defining the outer hull        
        print '\ncut number %d' % ii
        Z, outerhull, cutpoints = cut
        OR = ''
        for ii, cutpoint in enumerate(outerhull):
            xcut = cutpoint[0]
            ycut = cutpoint[1]
            print "\t%4s\t(%s %s %8.3f)\tAND\t(%s %s %8.3f)" % (OR,
                                                                xname,
                                                                xdir, xcut,
                                                                yname,
                                                                ydir, ycut)
            OR = 'OR'
        
    # -------------------------------------------------------------
    # roc plot
    # -------------------------------------------------------------
    print "\n\t=== plot ROC ==="	
    croc = TCanvas("fig_example2_ROC", "RGS", 700, 500, 500, 500)
    croc.cd()
    hist.Draw()
    croc.Update()

    # saveplots
    cmass.SaveAs('.png')    
    croc.SaveAs(".png")    
    sleep(5)
# ---------------------------------------------------------------------
main()



