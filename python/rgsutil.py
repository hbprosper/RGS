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
from ROOT import TPolyLine, kRed, kBlack
#-----------------------------------------------------------------------------
def error(message):
    print "** %s" % message
    exit(0)

# Return the name of a file with the extension and path stripped away.
def nameonly(s):
    import posixpath
    return posixpath.splitext(posixpath.split(s)[1])[0]

def count(filename, treename):
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
    records = map(split, filter(lambda x: x[0]!='#', open(varfilename)))
    return records
# ----------------------------------------------------------------------------    
class OuterHull:
    def __init__(self, xmin, xmax, ymin, ymax,
                 color=kRed,
                 hullcolor=kBlack,
                 hullwidth=2):
        self.cuts = []
        self.xmin, self.xmax, self.ymin, self.ymax = xmin, xmax, ymin, ymax
        self.color=color
        self.hullcolor=hullcolor
        self.hullwidth=hullwidth
        self.plots = []
        
    def add(self, significance, x, y, xdir=1, ydir=1):
        # sort all cut-points in this ladder cut
        # in increasing y value)
        cutpoints = [None]*len(y)
        for ii in xrange(len(y)):
            cutpoints[ii] = (y[ii], x[ii])
        cutpoints.sort()

        # find outer hull by picking cut-points such that
        # x value decreases monotonically
        outerhull = [cutpoints[0]]  # start with cut-point with lowest y value 
        for ii in xrange(1, len(cutpoints)):
            y0, x0 = outerhull[-1]
            y1, x1 = cutpoints[ii]
            # if the x-cut direction is > and the current cut-point
            # has a lower x value than that of the previous cut-point
            # then this point is on the outer hull.
            #
            # if the x-cut direction is < the condition is reversed
            if xdir > 0:
                if x1 < x0:
                    outerhull.append(cutpoints[ii])
            else:
                if x1 > x0:
                    outerhull.append(cutpoints[ii])
                                        
        # place significance in position 1, so that outer hulls can be
        # sorted according to the user-defined significance measure
        self.cuts.append((significance, outerhull, cutpoints))    

    # return outer hull of specified ladder cut
    def __call__(self, point=0):
        cuts = self.cuts
        # check for sensible point index
        if point < 0: return None
        if point > len(cuts)-1: return None
        # order outer hulls according to user-defined significances
        cuts.sort()
        cuts.reverse()
        return cuts[point]

    def plot(self, cutpoint, xmax, ymax, color):
        x = array('d')
        y = array('d')
        yy, xx = cutpoint
        y.append(ymax); x.append(xx)
        y.append(yy);   x.append(xx)
        y.append(yy);   x.append(xmax)
        poly = TPolyLine(len(x), x, y)
        poly.SetLineWidth(1)
        poly.SetLineColor(color)
        return poly
        
    def draw(self, point=0, hullcolor=kBlack, plotall=False):
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
                                  xmax, ymax,
                                  color))
                plots[-1].Draw('l same')

        # plot outer hull of current ladder
        for ii, cutpoint in enumerate(outerhull):
            if ii == 0:
                ymax, xx = outerhull[ii+1]
            elif ii < len(outerhull)-1:
                yy, xmax = outerhull[ii-1]                    
                ymax, xx = outerhull[ii+1]
            else:
                yy, xmax = outerhull[ii-1]
                ymax = self.ymax

            plots.append(plot(cutpoint,
                              xmax, ymax,
                              self.hullcolor))
            plots[-1].SetLineWidth(self.hullwidth)
            plots[-1].Draw('l same')                
