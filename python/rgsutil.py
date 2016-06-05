#-----------------------------------------------------------------------------
# File: rgsutil
# Description: A collection of RGS utilities.
#
# Created: 18-Jan-2015 Harrison B. Prosper and Sezen Sekmen
# Udpated: 21-Oct-2015 HBP & SS - allow flexibility of plotting the outer
#                      hull of any ladder cut
#          28-May-2016 HBP - rename LadderPlot OuterHull
#-----------------------------------------------------------------------------
import os, sys, re
from array import array
from string import split, strip, atoi, atof, replace, joinfields
from math import *
from ROOT import *
#-----------------------------------------------------------------------------
def getCWD():
    return split(os.environ['PWD'],'/')[-1]

def error(message):
    print "** %s" % message
    exit(0)

# Return the name of a file with the extension and path stripped away.
def nameonly(s):
    import posixpath
    return posixpath.splitext(posixpath.split(s)[1])[0]

def getEntries(filename, treename):
    tfile = TFile(filename)
    if not tfile.IsOpen():
        sys.exit('** cannot open file %s' % filename)
    t = tfile.Get(treename)
    if not t:
        sys.exit('** cannot get tree %s' % treename)
    n = t.GetEntries()
    return n    
# ----------------------------------------------------------------------------
def getCutDirections(varfilename):
    if not os.path.exists(varfilename):
        error("unable to open variables file %s" % varfilename)
    records = map(split,
                  filter(lambda x: x[0] != '#',
                        filter(lambda x: x != '',
                               map(strip, open(varfilename)))))       
    return records
# ----------------------------------------------------------------------------    
class OuterHull:
    def __init__(self, xmin, xmax, ymin, ymax,
                 cutdirs=[('x', '>'), ('y', '>')],
                 color=kRed,
                 hullcolor=kBlack,
                 hullwidth=2):
        self.cuts = []
        self.cutdirs = cutdirs
        self.xmin, self.xmax, self.ymin, self.ymax = xmin, xmax, ymin, ymax
        self.color=color
        self.hullcolor=hullcolor
        self.hullwidth=hullwidth
        self.plots = []
        
    def add(self, significance, x, y):
        cutdirs = self.cutdirs
        # sort all cut-points in this ladder cut
        # in increasing or decreasing order of y value
        # depending on the cut directions
        cutpoints = [None]*len(y)
        for ii in xrange(len(y)):
            cutpoints[ii] = [y[ii], x[ii]]
        cutpoints.sort()
        
        xname, xdir = cutdirs[0]
        yname, ydir = cutdirs[1]
        if (xdir == '>' and ydir == '<') or (xdir == '<' and ydir == '<'):
            cutpoints.reverse()
        
        # find outer hull by picking cut-points such that
        # x value decreases or increases monotonically, depending on
        # the cut directions.
        outerhull = [cutpoints[0]]
        for ii in xrange(1, len(cutpoints)):
            y0, x0 = outerhull[-1]  # previous point
            y1, x1 = cutpoints[ii]  # current point
            if   (xdir == '>' and ydir == '>') or (xdir == '>' and ydir == '<'):
                if x1 < x0: outerhull.append(cutpoints[ii])
            else:
                if x1 > x0: outerhull.append(cutpoints[ii])
                                        
        # place significance in position 1, so that outer hulls can be
        # sorted according to the user-defined significance measure
        self.cuts.append((significance, outerhull, cutpoints))
        self.sort = True
        
    # return outer hull of specified ladder cut
    def __call__(self, point=0):
        cuts = self.cuts
        # check for sensible point index
        if point < 0: return None
        if point > len(cuts)-1: return None
        if self.sort:
            self.sort = False
            # order outer hulls according to user-defined significances        
            self.cuts.sort()
            self.cuts.reverse()
        Z, outerhull, cutpoints = cuts[point]    
        for i, (y, x) in enumerate(outerhull):
            outerhull[i] = (x, y)
        return (Z, outerhull, cutpoints)

    def plot(self, cutpoint, xmin, xmax, ymin, ymax, color, plotit=True):
        xname, xdir = self.cutdirs[0]
        yname, ydir = self.cutdirs[1]
        x = array('d')
        y = array('d')
        yy, xx = cutpoint
        xx = min(xx, xmax); xx = max(xx, xmin)
        yy = min(yy, ymax); yy = max(yy, ymin)
        
        if xdir == '>':
            x.append(xmax); x.append(xx); x.append(xx)
            y.append(yy);   y.append(yy)
            if ydir == '>':
                y.append(ymax)
            else:
                y.append(ymin)
        else:
            x.append(xmin); x.append(xx); x.append(xx)
            y.append(yy);   y.append(yy)
            if ydir == '>':
                y.append(ymax)
            else:
                y.append(ymin)

        if plotit:
            poly = TPolyLine(len(x), x, y)
            poly.SetLineWidth(1)
            poly.SetLineColor(color)
            return poly
        else:
            return (x, y)
        
    def draw(self, point=0, option='l same', hullcolor=kBlack, plotall=False):
        cuts = self.__call__(point)
        if cuts == None: return
        significance, outerhull, cutpoints = cuts
            
        xmin, xmax, ymin, ymax = self.xmin, self.xmax, self.ymin, self.ymax
        color = self.color
        plot  = self.plot
        plots = self.plots
        if hullcolor != None:
            self.hullcolor = hullcolor
        hullcolor = self.hullcolor

        # plot cut-points of current ladder, if requested
        if plotall:
            for cutpoint in cutpoints:
                plots.append(plot(cutpoint,
                                  self.xmin, self.xmax,
                                  self.ymin, self.ymax,
                                  color))
                plots[-1].Draw(option)

        # plot outer hull of current ladder
        x = array('d')
        y = array('d')
        xx, yy = plot(outerhull[0],
                      self.xmin, self.xmax,
                      self.ymin, self.ymax,
                      color,
                      False)
        x.append(xx[0]); y.append(yy[0])
        x.append(xx[1]); y.append(yy[1])
        x.append(xx[2]); y.append(yy[2])
        for cutpoint in outerhull[1:]:
            xx, yy = plot(cutpoint,
                          self.xmin, self.xmax,
                          self.ymin, self.ymax,
                          color,
                          False)
            # the y-value of the previous point should be
            # equal to the y-value of the current point.
            # note: we skip the first point of the triplet
            # because it is a doppelganger of the previous point.
            y[-1] = yy[0]
            x.append(xx[1]); y.append(yy[1])
            x.append(xx[2]); y.append(yy[2]) 
            
        hull = TPolyLine(len(x), x, y)
        hull.SetLineWidth(2)
        hull.SetLineColor(hullcolor)
        hull.Draw(option)           
        plots.append(hull)
#---------------------------------------------------------------------------
