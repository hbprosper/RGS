import os, sys, re
from array import array
from rgsutil import *
from string import *
from math import *
from ROOT import *
# -----------------------------------------------------------------------
# used in examples
# -----------------------------------------------------------------------
def fillZ1massZ2mass(cut):
    xbins =  50
    xmin  =   0.0
    xmax  = 150.0

    ybins =  50
    ymin  =   0.0
    ymax  = 150.0    

    # Make a canvas and set its margins
    cmass = TCanvas("cmass", "cmass", 10,10, 425, 400)
    cmass.SetBottomMargin(0.15)
    cmass.SetLeftMargin(0.15)
    
    # -- background
    hb = mkhist2("hb",
                 "#font[12]{m}_{Z1} (GeV)",
                 "#font[12]{m}_{Z2} (GeV)",
                 xbins, xmin, xmax,
                 ybins, ymin, ymax,
                 color=kMagenta+1)
    hb.GetYaxis().SetTitleOffset(1.25)
    
    bntuple = Ntuple('../../data/ntuple_ZZ4L.root', 'Analysis')
    total = 0
    for ii, event in enumerate(bntuple):
        if cut(event): continue

        hb.Fill(event.Z1mass, event.Z2mass)
        if total % 1000 == 0:
            cmass.cd()
            hb.Draw('box')
            cmass.Update()
            gSystem.ProcessEvents()
        total += 1

    cmass.cd()
    hb.Draw('box')
    cmass.Update()
    gSystem.ProcessEvents()
    
    # -- signal
    hs = mkhist2("hs",
                 "#font[12]{m}_{Z1} (GeV)",
                 "#font[12]{m}_{Z2} (GeV)",
                 xbins, xmin, xmax,
                 ybins, ymin, ymax,
                 color=kAzure+1)
    hs.GetYaxis().SetTitleOffset(1.25)
    
    sntuple = Ntuple(['../../data/ntuple_HZZ4L_VBF.root',
                      '../../data/ntuple_HZZ4L.root'], 'Analysis')
    total   = 0
    for event in sntuple:
        if cut(event): continue

        hs.Fill(event.Z1mass, event.Z2mass, event.weight)
        if total % 1000 == 0:
            cmass.cd()
            hs.Draw('box')
            cmass.Update()
            gSystem.ProcessEvents()    
        total  += 1
        
    cmass.cd()
    hs.Draw('box')
    cmass.Update()
    gSystem.ProcessEvents()
    hs.Scale(1.0/hs.Integral())
    hb.Scale(1.0/hb.Integral())
    return (cmass, hs, hb)

def writeHZZResults(filename, varfilename, ntuple, variables,
                    bestrow, bestZ, totals, outerhull=None):
    
    cutdir  = {}
    cutdirs = []
    for t in getCutDirections(varfilename):
        token = t[0]
        if token == '\\':
            continue
        else:
            cutdirs.append(t)
            cutdir[token] = t[1]
            
    ntuple.read(bestrow)
    
    out = open(filename, 'w')

    print 
    record = "Yields before optimization"
    out.write('%s\n' % record); print record
    
    record = "\tVBF:    %10.3f +/- %-10.1f" % (totals[0][0],
                                                   totals[0][1])
    out.write('%s\n' % record); print record

    record = "\tggF:    %10.3f +/- %-10.1f" % (totals[1][0],
                                                   totals[1][1])
    out.write('%s\n' % record); print record    

    record = "\tZZ:     %10.3f +/- %-10.1f" % (totals[2][0],
                                                   totals[2][1])
    out.write('%s\n' % record); print record

    s = totals[0][0]
    b = totals[1][0] + totals[2][0]
    Z = signalSignificance(s, b)

    record = "\nZ values"
    out.write('%s\n' % record); print record
    
    record = "  before optimization:  %10.3f" % Z
    out.write('%s\n' % record); print record

    record = "  after optimization:   %10.3f" % bestZ
    out.write('%s\n' % record); print record    


    record = "Best cuts"
    out.write('\n%s\n' % record); print; print record

    bestcuts = {}
    
    if outerhull:
        OR = ''
        for ii, cutpoint in enumerate(outerhull):
            xname, xdir = cutdirs[0]
            yname, ydir = cutdirs[1]
            xcut = cutpoint[0]
            ycut = cutpoint[1]
            record = "\t%4s\t(%s %s %8.3f)\tAND\t(%s %s %8.3f)" % (OR,
                                                                xname,
                                                                xdir, xcut,
                                                                yname,
                                                                ydir, ycut)
            OR = 'OR'
            out.write('%s\n' % record); print record

        out.write('\n' % record); print
        
    for name, count in variables:    
        if name[0:5] in ['count', 'fract', 'cutpo']: continue
        var = ntuple(name)
        bestcuts[name] = var
        if type(var) == type(0.0):
            rec = '%3s %6.2f' % (cutdir[name], var)
            record = "\t%-10s\t%10s" % (name, rec)
            out.write('%s\n' % record); print record
        else:
            record = "\t%-10s\t%10.2f\t%10.2f" % (name,
                                                    min(var[0], var[1]),
                                                    max(var[0], var[1]))
            out.write('%s\n' % record); print record            
    print
    out.write('\n')

    record = "Yields after optimization (and relative efficiencies)"
    out.write('%s\n' % record); print record    

    for name, count in variables:        
        if not (name[0:5] in ['count', 'fract']): continue
        var = ntuple(name)
        record = "\t%-15s %10.3f" % (name, var)
        out.write('%s\n' % record); print record                    
        if name[0:5] == "fract":
            print
            out.write('\n')
    out.close()
    return bestcuts
# ---------------------------------------------------------------------
def drawHZZLegend(stitle="H #rightarrow 4#font[12]{l}",
                    btitle="ZZ #rightarrow 4#font[12]{l}",
                left=True,
                      postfix=''):                       
    lgoffset1 = 0.595
    lgoffset2 = 0.583
    if left:
        lgoffset1 = 0.23
        lgoffset2 = 0.21
        
    # Title
    t = TLatex(lgoffset1, 0.85, "njets #geq 2")
    t.SetTextSize(0.045)
    t.SetTextFont(42)
    t.SetNDC()
    t.Draw("same")
    SetOwnership(t, 0)    

    # Make the histogram legend    
    hs_d = TH1D('hs_d%s' % postfix, '', 4,0,4)
    hs_d.SetMarkerStyle(20)
    hs_d.SetMarkerColor(kMagenta+1)
    hs_d.SetMarkerSize(0.8)
    SetOwnership(hs_d, 0)
    
    hb_d = TH1D('hb_d%s' % postfix, '', 4,0,4)
    hb_d.SetMarkerStyle(20)    
    hb_d.SetMarkerColor(kAzure+1)
    hb_d.SetMarkerSize(0.8)    
    SetOwnership(hb_d, 0)
    
    l = TLegend(lgoffset2, 0.67, 0.90, 0.82)
    l.SetBorderSize(0)
    l.SetFillStyle(0000)
    l.SetTextSize(0.045)
    l.SetTextFont(42)
    l.AddEntry(hs_d, stitle, 'p')
    l.AddEntry(hb_d, btitle, 'p')
    l.Draw("same")
    SetOwnership(l, 0)
# -----------------------------------------------------------------------
# SUSY
# -----------------------------------------------------------------------
def fill_MR_R2(Cut, name='SUSY'):

    xbins =   50
    xmin  =    0.0
    xmax  = 3000.0

    ybins =   50
    ymin  =    0.0
    ymax  =    0.5    
    
    cmass = TCanvas("h_%s" % name, name, 10, 10, 500, 500)    
    
    # -- background
    hb = mkhist2("hb",
                 "M_{R} (GeV)",
                 "R^{2}",
                 xbins, xmin, xmax,
                 ybins, ymin, ymax,
                 color=kMagenta+1)
    hb.GetYaxis().SetTitleOffset(1.80)
    
    bntuple = Ntuple('../../data/ttbar.root', 'RGSinput')
    count = 0
    for event in bntuple:
        if Cut(event): continue
            
        hb.Fill(event.MR, event.R2, event.weight)
        if count % 1000 == 0:
            cmass.cd()
            hb.Draw('box')
            cmass.Update()
            gSystem.ProcessEvents()
        count += 1
        if count >= 10000: break
            
    # -- signal
    hs = mkhist2("hs",
                 "M_{R} (GeV)",
                 "R^{2}",
                 xbins, xmin, xmax,
                 ybins, ymin, ymax,
                 color=kAzure+1)
    hs.GetYaxis().SetTitleOffset(1.80)
    
    sntuple = Ntuple('../../data/SUSY_gluino_m1354.root', 'RGSinput')
    count = 0
    for event in sntuple:
        if Cut(event): continue
            
        hs.Fill(event.MR, event.R2, event.weight)
        if count % 1000 == 0:
            cmass.cd()
            hs.Draw('box')
            cmass.Update()
            gSystem.ProcessEvents()
        count += 1
        if count >= 10000: break
            
    cmass.cd(2)
    hs.Draw('box')
    hb.Draw('box same')
    cmass.Update()
    gSystem.ProcessEvents()    
    return (cmass, hs, hb)

def writeSUSYResults(filename, varfilename, ntuple, variables,
                     bestrow, bestZ, outerhull):
    
    cutdir  = {}
    cutdirs = []
    for t in getCutDirections(varfilename):
        token = t[0]
        if token[0] == '\\':
            continue
        else:
            cutdirs.append(t)
            cutdir[token] = t[1]

    totals = ntuple.totals()
    ntuple.read(bestrow)
    
    out = open(filename, 'w')

    print 
    record = "Yields before optimization"
    out.write('%s\n' % record); print record
    
    record = "\tttbar:   %12.1f +/- %-10.1f" % (totals[0][0],
                                                   totals[0][1])
    out.write('%s\n' % record); print record

    record = "\tSUSY:    %12.1f +/- %-10.1f" % (totals[1][0],
                                                   totals[1][1])
    out.write('%s\n' % record); print record    

    b = totals[0][0]
    s = totals[1][0]
    Z = signalSignificance(s, b)

    record = "\nZ values"
    out.write('%s\n' % record); print record
    
    record = "  before optimization:  %10.3f" % Z
    out.write('%s\n' % record); print record

    record = "  after optimization:   %10.3f" % bestZ
    out.write('%s\n' % record); print record    

    record = "Best cuts"
    out.write('\n%s\n' % record); print; print record

    bestcuts = {}
    
    OR = ''
    for ii, cutpoint in enumerate(outerhull):
        xname, xdir = cutdirs[0]
        yname, ydir = cutdirs[1]
        xcut = cutpoint[0]
        ycut = cutpoint[1]
        record = "\t%4s\t(%s %s %8.3f)\tAND\t(%s %s %8.3f)" % (OR,
                                                            xname,
                                                            xdir, xcut,
                                                            yname,
                                                            ydir, ycut)
        OR = 'OR'
        out.write('%s\n' % record); print record

    out.write('%s\n' % record); print

    for name, cdir in cutdirs[2:]:    
        var = ntuple(name)
        bestcuts[name] = var
        if type(var) == type(0.0):
            rec = '%3s %6.2f' % (cutdir[name], var)
            record = "\t%-10s\t%10s" % (name, rec)
            out.write('%s\n' % record); print record
        else:
            record = "\t%-10s\t%10.2f\t%10.2f" % (name,
                                                    min(var[0], var[1]),
                                                    max(var[0], var[1]))
            out.write('%s\n' % record); print record            
    print
    out.write('\n')

    record = "Yields after optimization (and relative efficiencies)"
    out.write('%s\n' % record); print record    

    for name, count in variables:        
        if not (name[0:5] in ['count', 'fract']): continue
        var = ntuple(name)
        record = "\t%-15s %10.3f" % (name, var)
        out.write('%s\n' % record); print record                    
        if name[0:5] == "fract":
            print
            out.write('\n')
    out.close()
    return bestcuts
# ---------------------------------------------------------------------
def drawSUSYLegend(lgtitle='',
                   stitle="pp #rightarrow #tilde{g}#tilde{g}",
                   btitle="pp #rightarrow t#bar{t}",
                   left=True,
                      postfix=''):                       
    lgoffset1 = 0.595
    lgoffset2 = 0.583
    if left:
        lgoffset1 = 0.23
        lgoffset2 = 0.21
        
    # Title
    t = TLatex(lgoffset1, 0.85, lgtitle)
    t.SetTextSize(0.045)
    t.SetTextFont(42)
    t.SetNDC()
    t.Draw("same")
    SetOwnership(t, 0)    

    # Make the histogram legend    
    hs_d = TH1D('hs_d%s' % postfix, '', 4,0,4)
    hs_d.SetMarkerStyle(20)
    hs_d.SetMarkerColor(kMagenta+1)
    hs_d.SetMarkerSize(0.8)
    SetOwnership(hs_d, 0)
    
    hb_d = TH1D('hb_d%s' % postfix, '', 4,0,4)
    hb_d.SetMarkerStyle(20)    
    hb_d.SetMarkerColor(kAzure+1)
    hb_d.SetMarkerSize(0.8)    
    SetOwnership(hb_d, 0)
    
    l = TLegend(lgoffset2, 0.67, 0.90, 0.82)
    l.SetBorderSize(0)
    l.SetFillStyle(0000)
    l.SetTextSize(0.045)
    l.SetTextFont(42)
    l.AddEntry(hs_d, stitle, 'p')
    l.AddEntry(hb_d, btitle, 'p')
    l.Draw("same")
    SetOwnership(l, 0)
