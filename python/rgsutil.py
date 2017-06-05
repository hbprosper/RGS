#-----------------------------------------------------------------------------
# File: rgsutil
# Description: A collection of RGS utilities.
#
# Created: 18-Jan-2015 Harrison B. Prosper and Sezen Sekmen
# Udpated: 21-Oct-2015 HBP & SS - allow flexibility of plotting the outer
#                      hull of any ladder cut
#          28-May-2016 HBP - rename LadderPlot OuterHull
#          02-Sep-2016 HBP - generalize OuterHull to permit drawing of
#                      multiple outer hulls.
#          25-Sep-2016 HBP - finally, bug in OuterHull fixed!
#-----------------------------------------------------------------------------
import os, sys, re
from array import array
from string import split, strip, atoi, atof, replace, joinfields
from math import *
from ROOT import *
#-----------------------------------------------------------------------------
TEXTFONT=42
NDIVX=510
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

def signalSignificance(s, b, which=1):
    # Compute various measures of significance
    Z = 0.0
    if which == 1:
        #   Z  = sign(LR) * sqrt(2*|LR|)
        # where LR = log(Poisson(s+b|s+b)/Poisson(s+b|b))
        if b > 1:
            Z = 2*((s+b)*log((s+b)/b)-s)
            absZ = abs(Z)
            if absZ != 0:
                Z = Z*sqrt(absZ)/absZ
            else:
                Z = 0.0
    return Z

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
                 cutdirs=[('x', '>'), ('y', '>')]):
        self.xmin, self.xmax, self.ymin, self.ymax = xmin, xmax, ymin, ymax
        self.cutdirs = cutdirs                
        self.cuts  = []
        self.plots = []
        self.sort  = True
        
    def add(self, significance, x, y):
        cutdirs = self.cutdirs
        # sort all cut-points in this ladder-cut
        # in increasing value of y
        cutpoints = [None]*len(y)
        for ii in xrange(len(y)):
            cutpoints[ii] = [y[ii], x[ii]]
        cutpoints.sort()
        
        # now switch order of x and y
        for ii, (yy, xx) in enumerate(cutpoints):
            cutpoints[ii] = [xx, yy]
                
        # find outer hull by picking cut-points such that
        # x value decreases or increases monotonically,
        # depending on the cut directions.
        xdir = cutdirs[0][1]
        ydir = cutdirs[1][1]
        
        outerhull = [cutpoints[0]]
        for ii in xrange(1, len(cutpoints)):
            x0, y0 = outerhull[-1]  # previous point
            x1, y1 = cutpoints[ii]  # current point
            if   xdir == '>'  and ydir == '>':
                if x1 < x0: outerhull.append(cutpoints[ii])
            elif xdir == '>'  and ydir == '<':
                if x1 < x0: outerhull.append(cutpoints[ii])
            elif xdir == '<'  and ydir == '>':
                if x1 > x0: outerhull.append(cutpoints[ii])
            elif xdir == '<'  and ydir == '<':
                if x1 > x0: outerhull.append(cutpoints[ii])
            else:
                pass            
                                        
        # place significance in position 1, so that outer hulls can be
        # sorted according to the user-defined significance measure
        self.cuts.append((significance, outerhull, cutpoints))

        # cutpoints will be sorted at the first call to retrieve an
        # outer hull of cuts-points
        self.sort = True
        
    # return outer hull of specified ladder cut
    def __call__(self, point=0):
        # check for sensible point index
        if point < 0: return None
        if point > len(self.cuts)-1: return None
        # first time around, sort cuts
        if self.sort:
            self.sort = False
            # order outer hulls according to user-defined significances        
            self.cuts.sort()
            self.cuts.reverse()
        Z, outerhull, cutpoints = self.cuts[point]    
        return (Z, outerhull, cutpoints)

    def plot(self, cutpoint, xmin, xmax, ymin, ymax, color,
             return_plot=True):
        xname, xdir = self.cutdirs[0]
        yname, ydir = self.cutdirs[1]
        x = array('d')
        y = array('d')
        xx, yy = cutpoint
        xx = min(xx, xmax); xx = max(xx, xmin)
        yy = min(yy, ymax); yy = max(yy, ymin)
        
        if xdir == '>':
            x.append(xmax); x.append(xx)
            y.append(yy);   y.append(yy)
            
            x.append(xx)
            if ydir == '>':
                y.append(ymax)
            else:
                y.append(ymin)
        else:
            x.append(xmin); x.append(xx)
            y.append(yy);   y.append(yy)
            
            x.append(xx)
            if ydir == '>':
                y.append(ymax)
            else:
                y.append(ymin)
                
        if return_plot:
            poly = TPolyLine(len(x), x, y)
            poly.SetLineWidth(1)
            poly.SetLineColor(color)
            return poly
        else:
            return (x, y)
        
    def draw(self, laddercuts,
             option='l same',
             color=kRed,
             hullcolor=kBlack,
             lstyle=1,
             lwidth=2,
             plotall=False):
        if laddercuts == None: return
        
        significance, outerhull, cutpoints = laddercuts
        
        # plot cut-points of current ladder, if requested
        if plotall:
            for cutpoint in cutpoints:
                self.plots.append(self.plot(cutpoint,
                                            self.xmin, self.xmax,
                                            self.ymin, self.ymax,
                                            color))
                self.plots[-1].Draw(option)
                
        # plot outer hull of current ladder
        x = array('d')
        y = array('d')
        xx, yy = self.plot(outerhull[0],
                           self.xmin, self.xmax,
                           self.ymin, self.ymax,
                           color,
                           return_plot=False)
        x.append(xx[0]); y.append(yy[0])
        x.append(xx[1]); y.append(yy[1])
        x.append(xx[2]); y.append(yy[2])

        for cutpoint in outerhull[1:]:
            xx, yy = self.plot(cutpoint,
                               self.xmin, self.xmax,
                               self.ymin, self.ymax,
                               color,
                               return_plot=False)
            # the y-value of the previous point should be
            # equal to the y-value of the current point.
            # note: we skip the first point of the triplet
            # because it is a duplicate of the previous point.
            y[-1] = yy[0]
            x.append(xx[1]); y.append(yy[1])
            x.append(xx[2]); y.append(yy[2]) 

        hull = TPolyLine(len(x), x, y)
        hull.SetLineWidth(lwidth)
        hull.SetLineStyle(lstyle)
        hull.SetLineColor(hullcolor)
        hull.Draw(option)
        self.plots.append(hull)
#---------------------------------------------------------------------------
def setStyle():
    style = TStyle("rgsutil", "rgsutil")
    style.SetPalette(1)
    
    # For the canvases
    style.SetCanvasBorderMode(0)
    style.SetCanvasColor(kWhite)
    style.SetCanvasDefH(500) #Height of canvas
    style.SetCanvasDefW(500) #Width of canvas
    style.SetCanvasDefX(0)   #Position on screen
    style.SetCanvasDefY(0)

    # For the pads
    style.SetPadBorderMode(0)
    style.SetPadColor(kWhite)
    style.SetPadGridX(kFALSE)
    style.SetPadGridY(kFALSE)
    style.SetGridColor(kGreen)
    style.SetGridStyle(3)
    style.SetGridWidth(1)

    # For the frames
    style.SetFrameBorderMode(0)
    style.SetFrameBorderSize(1)
    style.SetFrameFillColor(0)
    style.SetFrameFillStyle(0)
    style.SetFrameLineColor(1)
    style.SetFrameLineStyle(1)
    style.SetFrameLineWidth(1)

    # For the histograms
    style.SetHistLineColor(kBlack)
    style.SetHistLineStyle(0)
    style.SetHistLineWidth(2)

    style.SetEndErrorSize(2)
    #style.SetErrorX(0.)

    style.SetMarkerSize(0.4)
    style.SetMarkerStyle(20)

    # For the fit/function:
    style.SetOptFit(1)
    style.SetFitFormat("5.4g")
    style.SetFuncColor(2)
    style.SetFuncStyle(1)
    style.SetFuncWidth(1)

    # For the date:
    style.SetOptDate(0)

    # For the statistics box:
    style.SetOptFile(0)
    style.SetOptStat("")
    # To display the mean and RMS:
    # style.SetOptStat("mr") 
    style.SetStatColor(kWhite)
    style.SetStatFont(TEXTFONT)
    style.SetStatFontSize(0.03)
    style.SetStatTextColor(1)
    style.SetStatFormat("6.4g")
    style.SetStatBorderSize(1)
    style.SetStatH(0.2)
    style.SetStatW(0.3)

    # Margins:
    style.SetPadTopMargin(0.05)
    style.SetPadBottomMargin(0.16)
    style.SetPadLeftMargin(0.20)
    style.SetPadRightMargin(0.10)

    # For the Global title:
    style.SetOptTitle(0) 
    style.SetTitleFont(TEXTFONT)
    style.SetTitleColor(1)
    style.SetTitleTextColor(1)
    style.SetTitleFillColor(10)
    style.SetTitleFontSize(0.05)

    # For the axis titles:
    style.SetTitleColor(1, "XYZ")
    style.SetTitleFont(TEXTFONT, "XYZ")
    style.SetTitleSize(0.05, "XYZ")
    style.SetTitleXOffset(1.25)
    style.SetTitleYOffset(1.40)

    # For the axis labels:
    style.SetLabelColor(1, "XYZ")
    style.SetLabelFont(TEXTFONT, "XYZ")
    style.SetLabelOffset(0.020, "XYZ")
    style.SetLabelSize(0.05, "XYZ")

    # For the axis:
    style.SetAxisColor(1, "XYZ")
    style.SetStripDecimals(kTRUE)
    style.SetTickLength(0.03, "XYZ")
    style.SetNdivisions(NDIVX, "XYZ")
    # To get tick marks on the opposite side of the frame
    style.SetPadTickX(1)  
    style.SetPadTickY(1)

    # Change for log plots:
    style.SetOptLogx(0)
    style.SetOptLogy(0)
    style.SetOptLogz(0)

    # Postscript options:
    style.SetPaperSize(20.,20.)
    style.cd()
#------------------------------------------------------------------------------
def getarg(args, key, d):
    if args.has_key(key):
        return args[key]
    else:
        return d
#------------------------------------------------------------------------------
def mkhist1(hname, xtitle, ytitle, nbins, xmin, xmax, **args):
    ymin   = getarg(args, 'ymin', None)
    ymax   = getarg(args, 'ymax', None)
    color  = getarg(args, 'color',   kBlack)
    lstyle = getarg(args, 'lstyle',  1)
    lwidth = getarg(args, 'lwidth',  1)
    ndivx  = getarg(args, 'ndivx',   505)
    ndivy  = getarg(args, 'ndivy',   510)

    h = TH1F(hname, "", nbins, xmin, xmax)		
    h.SetLineColor(color)
    h.SetLineStyle(lstyle)
    h.SetLineWidth(lwidth)

    h.SetMarkerSize(0.8)
    h.SetMarkerColor(color)
    h.SetMarkerStyle(20)

    h.GetXaxis().SetTitle(xtitle)
    h.GetXaxis().SetTitleOffset(1.2);
    h.GetXaxis().SetLimits(xmin, xmax)
    h.SetNdivisions(ndivx, "X")

    h.GetYaxis().SetTitle(ytitle)
    h.GetYaxis().SetTitleOffset(1.6)
    if ymin != None: h.SetMinimum(ymin)
    if ymax != None: h.SetMaximum(ymax)
    h.SetNdivisions(ndivy, "Y")
    return h
#------------------------------------------------------------------------------
def mkhist2(hname, xtitle, ytitle,
        nbinx, xmin, xmax,	
        nbiny, ymin, ymax, **args):
    color  = getarg(args, 'color',   kBlack)
    mstyle = getarg(args, 'mstyle',  20)
    msize  = getarg(args, 'msize',   0.5)
    ndivx  = getarg(args, 'ndivx',   505)
    ndivy  = getarg(args, 'ndivy',   505)

    h = TH2F(hname, "", nbinx, xmin, xmax, nbiny, ymin, ymax)
    h.SetLineColor(color)
    h.SetMarkerColor(color)
    h.SetMarkerSize(msize)
    h.SetMarkerStyle(mstyle)

    #h.GetXaxis().CenterTitle()
    h.GetXaxis().SetTitle(xtitle)
    h.GetXaxis().SetTitleOffset(1.3)
    h.SetNdivisions(ndivx, "X")

    #h.GetYaxis().CenterTitle()
    h.GetYaxis().SetTitle(ytitle)
    h.GetYaxis().SetTitleOffset(1.3)
    h.SetNdivisions(ndivy, "Y")
    return h
#------------------------------------------------------------------------------
class Buffer:

    def __init__(self, buffer, buffermap, variable):
        self.buffer = buffer
        self.buffermap = buffermap
        self.variable = variable

    def __getattr__(self, variable):
        if self.buffermap.has_key(variable):
            jj = self.buffermap[variable]
            return self.buffer[jj].__getattribute__(variable)
        else:
            raise AttributeError(variable)

    def __call__(self, variable):
        return self.__getattr__(variable)
    
    def __str__(self):
        s = 'Event variables:\n'
        for tname, name, maxcount in self.variable:
            s += "  %-12s %-24s:\t" % (tname, name)
            s += "%s\n" % self.__getattr__(name)
        return s

class Ntuple:
    # "self" is Python's equivalent of the "this" pointer in C++
    # self points to the memory allocated for the object

    def __init__(self, filename, treename, firstrow=0, nrows=None):

        # cache inputs
        self.status = 0

        from random import randint
        self.postfix  = randint(1, 1000000)
        if type(filename) == type(""):
            self.filename = [filename]
        else:
            self.filename = filename

        self.currentTreeNumber = -1
        self.treename = treename
        self.nrows = nrows

        # make sure all files exist
        fnames = []
        for fname in self.filename:
            if not os.path.exists(fname):
                print "** Ntuple *** "\
                      "root file %s not found" % fname
                sys.exit(0)

            # check file size
            size = os.path.getsize(fname)
            if size == 0:
                print "== ZERO length == %s" % fname
                continue
            else:
                fnames.append(fname)

        self.filename = fnames
        if len(self.filename) == 0:
            self.status = -3
            return

        # create a chain of files
        self.chain = TChain(treename)
        if not self.chain:
            print "*** Ntuple *** can't create chain %s" % treename
            sys.exit(0)

        # ------------------------------------
        # add files to chain
        # ------------------------------------
        for fname in self.filename:
            try:
                self.chain.Add(fname)
            except:
                print "*** problem with file %s" % fname

        self.entries = self.chain.GetEntries()
        self.tree = self.chain
        tree = self.tree
            
        if tree == None:
            print "** problem accessing Tree - "\
              "perhaps the name %s is wrong?" % treename
            sys.exit(0)

        # get names of variables from root file
        branches  = tree.GetListOfBranches()

        # get number of variables
        try:
            nbranches = branches.GetEntries()
        except:
            print "** ====>  problem accessing branches\n"
            self.status = -1
            return

        bnamemap = {}        
        self.vars = []
        for i in xrange(nbranches):
            # get the ith branch (aka variable)
            bname = branches[i].GetName()
            			
            # just in case, check for duplicates!
            if bnamemap.has_key(bname):
                print '** duplicate branch name: %s' % bname
                continue
            else:
                bnamemap[bname] = 1
            
            # assume one leaf/branch
            # Get all leaves associated with this branch
            leaves = branches[i].GetListOfLeaves()
            if leaves == None:
                print "No leaves found!"
                sys.exit(0)

            leaf = leaves[0]
            if leaf == None:
                print "No leaf found"
                sys.exit(0)

            leafname = leaf.GetName()

            # get leaf type (int, float, double, etc.)
            tname = leaf.GetTypeName()

            #check for leaf counter
            flag = Long(0)
            leafcounter = leaf.GetLeafCounter(flag)
            if leafcounter:
                maxcount = leafcounter.GetMaximum()
            else:
                maxcount = leaf.GetLen()

            # store type and variable name
            self.vars.append( (tname, bname, maxcount) )
                
        nlen = len(self.vars)

        # create a map of variable name to column number
        self.varmap = {}
        for ind, var in enumerate(self.vars):
            self.varmap[var] = ind

        nentries = self.entries
        if self.nrows != None:
            self.entries = min(self.nrows, nentries) 
        else:			
            self.entries = nentries
        # ------------------------------------
        # set up branches as a struct
        # Root has a limit on how long a
        # string it can cope with, so split
        # into multiple strings
        # ------------------------------------

        bufferCount = 0
        newBuffer = True
        rec = ""
        bufferName = ""
        maxlength  = 2000
        self.buffermap  = {}
        self.buffer = []

        self.varname= []
        for count, (tname, name, maxcount) in enumerate(self.vars):
            self.varname.append((name, maxcount))
            
            # keep track of map from variable name to buffer count
            self.buffermap[name] = bufferCount

            if newBuffer:
                newBuffer = False
                bufferName = "S%d_%d" % (self.postfix, bufferCount)
                rec = "struct %s {" % bufferName

            if maxcount == 1:
                rec += "%s %s;" % (tname, name)
            else:				
                rec += "%s %s[%d];" % (tname, name, maxcount)

            if (len(rec) > maxlength) or \
                   (count >= len(self.vars)-1):
                rec += "};"
                newBuffer = True

                # evaluate it

                gROOT.ProcessLine(rec)

                # now import struct
                exec("from ROOT import %s" % bufferName)

                # add to list of buffers
                self.buffer.append(eval("%s()" % bufferName))

                # remember to update buffer count
                bufferCount += 1

        # create a generic event object
        self.event = Buffer(self.buffer, self.buffermap, self.vars)

        # Now that addresses are stable, give address of each variable
        for tname, name, maxcount in self.vars:
            jj = self.buffermap[name]
            tree.SetBranchAddress(name, AddressOf(self.buffer[jj], name))

        self.status = 0
        # initialize row number
        self.row = firstrow

    # destructor
    def __del__(self):
            pass

    def size(self):
        return int(self.entries)

    def numEntries(self):
        return int(self.entries)

    def totals(self):
        self.file = self.tree.GetCurrentFile()
        htotal = self.file.Get('total')
        if htotal:
            t = []
            for ii in xrange(htotal.GetNbinsX()):
                t.append((htotal.GetBinContent(ii+1),
                          htotal.GetBinError(ii+1)))
            return t
        else:
            return []
    
    def __len__(self):
        return int(self.entries)
    
    def read(self, row):
        localentry = self.chain.LoadTree(row)
        if self.chain.GetTreeNumber() != self.currentTreeNumber:
            self.currentTreeNumber = self.chain.GetTreeNumber()
            # Update branch addresses
            self.tree  = self.chain.GetTree()
            for tname, name, maxcount in self.vars:
                jj = self.buffermap[name]
                self.tree.SetBranchAddress(name,
                                           AddressOf(self.buffer[jj], name))

        self.tree.GetEntry(localentry)

    def treeNumber(self):
        return (self.currentTreeNumber,
                self.filename[self.currentTreeNumber])

    def good(self):
        return self.status == 0

    def get(self, variable):
        if self.buffermap.has_key(variable):
            jj = self.buffermap[variable]
            return self.buffer[jj].__getattribute__(variable)
        else:
            return None

    def __call__(self, variable):
        return self.get(variable)
    
    def __str__(self):
        rec = ''
        for ii, (tname, bname, maxcount) in enumerate(self.vars):
            rec += "\t%4d\t%-12s\t%-32s\t%d\n" % (ii, tname, bname, maxcount)
        return rec

    def ls(self):
        print self.__str__()

    def variables(self):
        return self.varname
    
    # Implement Python iterator protocol
    def __iter__(self):
        return self

    def next(self):
        if self.row > self.entries-1:
            self.row = 0
            raise StopIteration
        else:
            self.read(self.row)
            self.row += 1
            return self.event
        
#------------------------------------------------------------------------------
class Row:
    def __init__(self, rownumber, varmap, data):
        self.row = rownumber
        self.varmap= varmap
        self.data  = data
        self.items = map(lambda x: (x[1],x[0]), self.varmap.items())
        self.items.sort()

        # Initialize row counter
        self.col = 0
        self.maxcol = len(self.items)-1

    def __del__(self):
        pass

    def __call__(self, variable):
        if not self.varmap.has_key(variable): return None
        index = self.varmap[variable]
        if len(index) == 1:
            return self.data[index[0]]
        else:
            return self.data[index[0]:index[1]+1]
        
    def __str__(self):
        strrep = "row: %d\n" % self.row
        lname = ''
        for index, name in self.items:
            for ii in index:
                v = self.data[ii]
                strvalue = ''                
                if   type(v) == type(1.0):
                    strvalue += "%12.3f" % v
                elif type(v) == type(1):
                    strvalue += "%12d" % v                        
                elif type(v) == type(""):
                    strvalue += "%12s" % v
                strrep += "%4d %-16s %s\n" % (ii, name, strvalue)
        return strip(strrep)

    # Implement Python iterator protocol	
    def __iter__(self):
        return self

    def next(self):
        if self.col > self.maxcol:
            self.col = 0
            raise StopIteration
        else:
            index, name = self.items[self.col]
            if len(index) == 1: 
                value = self.data[index[0]]
            else:
                value = self.data[index[0]:index[1]+1]
            self.col += 1
            return (name, value)

    def __len__(self):
        return len(self.varmap)

    def __getitem__(self, key):
        if type(key) != type(1): return None
        if key < -len(self.varmap): return None
        if key > len(self.varmap)-1: return None
        return self.data[key]

def tonumber(x):
    try:
        y = atof(x)
    except:
        y = x
    return y

class Table:

    def __init__(self, filename, nrows=-1):
        try:
            myfile = open(filename, 'r')
        except:
            print "*** can't read file %s" % filename
            sys.exit(0)

        # Read header and check for array variables.
        # An array variable is identified by the syntax "name number"
        # Expand header by replicating the name appending to it an number
        
        records = myfile.readlines()
        t = split(records[0])
        
        self.varname= []
        self.header = []
        self.varmap = {}
        i = 0
        while i < len(t):
            name = t[i] # should be a name
            try:
                size = atoi(name)
                sys.exit('** wrong header syntax\n'\
                         '** found an integer in column %d where a string was '\
                         'expected\n' % i)
            except:
                pass
            
            # check if this is an array variable by looking ahead
            size = 1
            array_type = False
            if i < len(t)-1:
                try:
                    size = atoi(t[i+1])
                    array_type = True
                except:
                    pass

            if array_type:
                # this is an array variable so expand header names
                # (if size > 1),
                # but first cache starting position of array
                # and its size
                self.varmap[name] = []
                if size > 1:
                    self.varmap[name] = [len(self.header),
                                         len(self.header)+size-1]
                    for j in xrange(size):
                        newname = '%s[%d]' % (name, j)
                        self.header.append(newname)
                else:
                    self.varmap[name] = [len(self.header)]
                    self.header.append(name)
                i += 2
            else:
                # this is scalar variable so just add it to the
                # header as-is
                self.varmap[name] = [len(self.header)]
                self.header.append(name)
                i += 1
            # cache original names and size
            self.varname.append((name, size))
            
        ## # print some stuff
        ## print 
        ## lname = ''
        ## for ii, name in enumerate(self.header):
        ##     basename = split(name, '[')[0]
        ##     if basename != lname:
        ##         extra = self.varmap[basename]
        ##     else:
        ##         extra = ''
        ##     lname = basename
        ##     print '%5d\t%-16s\t%s' % (ii, self.header[ii], extra) 

        # load data into memory
        rownumber = 0
        self.data = []
        index = 1
        for record in records[1:]:
            # Convert to numbers
            record = map(tonumber, split(record))
            self.data.append(record)
            rownumber += 1
            if nrows > 0:
                if rownumber >= nrows:
                    break
        myfile.close()

        # Initialize row counter
        self.rownumber = 0
        self.maxrow = len(self.data)-1
        self.maxcol = len(self.header)-1
        
        # Create a name to index map for rows
        self.rowmap = {} # empty map
        for index in xrange(len(self.data)):
            self.rowmap['%d' % index] = index            

    def __del__(self):
        pass

    def __call__(self, rownumber, variable=None):
        if rownumber < 0: return None
        if rownumber > self.maxrow: return None
        if variable == None:
            return Row(rownumber, self.varmap, self.data[rownumber])
        else:
            if not self.varmap.has_key(variable): return None
            index = self.varmap[variable]
            if len(index) == 1:
                return self.data[rownumber][index[0]]
            else:
                return self.data[rownumber][index[0]:index[1]+1]

    # Implement Python iterator protocol
    def __iter__(self):
        return self

    def next(self):
        if self.rownumber > self.maxrow:
            self.rownumber = 0
            raise StopIteration
        else:
            data = Row(self.rownumber,
                       self.varmap,
                       self.data[self.rownumber])
            self.rownumber += 1
            return data

    def row(self, rownumber):
        if rownumber > self.maxrow: return None
        data = Row(rownumber, self.varmap, self.data[rownumber])
        return data

    def variables(self):
        return self.varname

    def numRows(self):
        return self.maxrow+1

    def numColumns(self):
        return len(self.header)

    def __len__(self):
        return self.numRows()

    def __getitem__(self, key):
        if type(key) != type(1): return None
        if key < -(self.maxrow+1): return None
        if key > self.maxrow: return None
        return Row(key, self.varmap, self.data[key])


PDGCode = '''
//-----------------------------------------------------------------------------
// pdgid to name map
// HBP
//-----------------------------------------------------------------------------
#include <map>
//-----------------------------------------------------------------------------
struct PDG
{
  std::map<int, std::string> namemap;

  PDG() : namemap(std::map<int, std::string>())
  {
      namemap[1]	= "d";
      namemap[-1]	= "d~";
      namemap[2]	= "u";
      namemap[-2]	= "u~";
      namemap[3]	= "s";
      namemap[-3]	= "s~";
      namemap[4]	= "c";
      namemap[-4]	= "c~";
      namemap[5]	= "b";
      namemap[-5]	= "b~";
      namemap[6]	= "t";
      namemap[-6]	= "t~";
      namemap[7]	= "b'";
      namemap[-7]	= "b'~";
      namemap[8]	= "t'";
      namemap[-8]	= "t'~";
      namemap[11]	= "e^-";
      namemap[-11]	= "e^+";
      namemap[12]	= "nu_e";
      namemap[-12]	= "nu_e~";
      namemap[13]	= "mu^-";
      namemap[-13]	= "mu^+";
      namemap[14]	= "nu_mu";
      namemap[-14]	= "nu_mu~";
      namemap[15]	= "tau^-";
      namemap[-15]	= "tau^+";
      namemap[16]	= "nu_tau";
      namemap[-16]	= "nu_tau~";
      namemap[17]	= "tau'^-";
      namemap[-17]	= "tau'^+";
      namemap[18]	= "nu_tau'";
      namemap[-18]	= "nu_tau'~";
      namemap[21]	= "g";
      namemap[22]	= "gamma";
      namemap[23]	= "Z^0";
      namemap[24]	= "W^+";
      namemap[-24]	= "W^-";
      namemap[25]	= "H_1^0";
      namemap[32]	= "Z_2^0";
      namemap[33]	= "Z_3^0";
      namemap[34]	= "W_2^+";
      namemap[-34]	= "W_2^-";
      namemap[35]	= "H_2^0";
      namemap[36]	= "H_3^0";
      namemap[37]	= "H^+";
      namemap[-37]	= "H^-";
      namemap[39]	= "G";
      namemap[41]	= "R^0";
      namemap[-41]	= "R~^0";
      namemap[42]	= "LQ_c";
      namemap[-42]	= "LQ_c~";
      namemap[51]	= "H_L^0";
      namemap[52]	= "H_1^++";
      namemap[-52]	= "H_1^--";
      namemap[53]	= "H_2^+";
      namemap[-53]	= "H_2^-";
      namemap[54]	= "H_2^++";
      namemap[-54]	= "H_2^--";
      namemap[55]	= "H_4^0";
      namemap[-55]	= "H_4~^0";
      namemap[81]	= "generator-specific+81";
      namemap[-81]	= "generator-specific-81";
      namemap[82]	= "generator-specific+82";
      namemap[-82]	= "generator-specific-82";
      namemap[83]	= "generator-specific+83";
      namemap[-83]	= "generator-specific-83";
      namemap[84]	= "generator-specific+84";
      namemap[-84]	= "generator-specific-84";
      namemap[85]	= "generator-specific+85";
      namemap[-85]	= "generator-specific-85";
      namemap[86]	= "generator-specific+86";
      namemap[-86]	= "generator-specific-86";
      namemap[87]	= "generator-specific+87";
      namemap[-87]	= "generator-specific-87";
      namemap[88]	= "generator-specific+88";
      namemap[-88]	= "generator-specific-88";
      namemap[89]	= "generator-specific+89";
      namemap[-89]	= "generator-specific-89";
      namemap[90]	= "generator-specific+90";
      namemap[-90]	= "generator-specific-90";
      namemap[91]	= "generator-specific+91";
      namemap[-91]	= "generator-specific-91";
      namemap[92]	= "generator-specific+92";
      namemap[-92]	= "generator-specific-92";
      namemap[93]	= "generator-specific+93";
      namemap[-93]	= "generator-specific-93";
      namemap[94]	= "generator-specific+94";
      namemap[-94]	= "generator-specific-94";
      namemap[95]	= "generator-specific+95";
      namemap[-95]	= "generator-specific-95";
      namemap[96]	= "generator-specific+96";
      namemap[-96]	= "generator-specific-96";
      namemap[97]	= "generator-specific+97";
      namemap[-97]	= "generator-specific-97";
      namemap[98]	= "generator-specific+98";
      namemap[-98]	= "generator-specific-98";
      namemap[99]	= "generator-specific+99";
      namemap[-99]	= "generator-specific-99";
      namemap[100]	= "generator-specific+100";
      namemap[-100]	= "generator-specific-100";
      namemap[101]	= "geantino";
      namemap[102]	= "charged-geantino";
      namemap[110]	= "reggeon";
      namemap[111]	= "pi^0";
      namemap[113]	= "rho(770)^0";
      namemap[115]	= "a_2(1320)^0";
      namemap[117]	= "rho_3(1690)^0";
      namemap[119]	= "a_4(2040)^0";
      namemap[130]	= "K_L^0";
      namemap[211]	= "pi^+";
      namemap[-211]	= "pi^-";
      namemap[213]	= "rho(770)^+";
      namemap[-213]	= "rho(770)^-";
      namemap[215]	= "a_2(1320)^+";
      namemap[-215]	= "a_2(1320)^-";
      namemap[217]	= "rho_3(1690)^+";
      namemap[-217]	= "rho_3(1690)^-";
      namemap[219]	= "a_4(2040)^+";
      namemap[-219]	= "a_4(2040)^-";
      namemap[221]	= "eta";
      namemap[223]	= "omega(782)";
      namemap[225]	= "f_2(1270)";
      namemap[227]	= "omega_3(1670)";
      namemap[229]	= "f_4(2050)";
      namemap[310]	= "K_S^0";
      namemap[311]	= "K^0";
      namemap[-311]	= "K~^0";
      namemap[313]	= "K*(892)^0";
      namemap[-313]	= "K*(892)~^0";
      namemap[315]	= "K*_2(1430)^0";
      namemap[-315]	= "K*_2(1430)~^0";
      namemap[317]	= "K*_3(1780)^0";
      namemap[-317]	= "K*_3(1780)~^0";
      namemap[319]	= "K*_4(2045)^0";
      namemap[-319]	= "K*_4(2045)~^0";
      namemap[321]	= "K^+";
      namemap[-321]	= "K^-";
      namemap[323]	= "K*(892)^+";
      namemap[-323]	= "K*(892)^-";
      namemap[325]	= "K*_2(1430)^+";
      namemap[-325]	= "K*_2(1430)^-";
      namemap[327]	= "K*_3(1780)^+";
      namemap[-327]	= "K*_3(1780)^-";
      namemap[329]	= "K*_4(2045)^+";
      namemap[-329]	= "K*_4(2045)^-";
      namemap[331]	= "eta'(958)";
      namemap[333]	= "phi(1020)";
      namemap[335]	= "f'_2(1525)";
      namemap[337]	= "phi_3(1850)";
      namemap[411]	= "D^+";
      namemap[-411]	= "D^-";
      namemap[413]	= "D*(2010)^+";
      namemap[-413]	= "D*(2010)^-";
      namemap[415]	= "D*_2(2460)^+";
      namemap[-415]	= "D*_2(2460)^-";
      namemap[421]	= "D^0";
      namemap[-421]	= "D~^0";
      namemap[423]	= "D*(2007)^0";
      namemap[-423]	= "D*(2007)~^0";
      namemap[425]	= "D*_2(2460)^0";
      namemap[-425]	= "D*_2(2460)~^0";
      namemap[431]	= "D_s^+";
      namemap[-431]	= "D_s^-";
      namemap[433]	= "D*_s^+";
      namemap[-433]	= "D*_s^-";
      namemap[435]	= "D*_s2(2573)^+";
      namemap[-435]	= "D*_s2(2573)^-";
      namemap[441]	= "eta_c(1S)";
      namemap[443]	= "J/psi(1S)";
      namemap[445]	= "chi_c2(1P)";
      namemap[511]	= "B^0";
      namemap[-511]	= "B~^0";
      namemap[513]	= "B*^0";
      namemap[-513]	= "B*~^0";
      namemap[515]	= "B*_2^0";
      namemap[-515]	= "B*_2~^0";
      namemap[521]	= "B^+";
      namemap[-521]	= "B^-";
      namemap[523]	= "B*^+";
      namemap[-523]	= "B*^-";
      namemap[525]	= "B*_2^+";
      namemap[-525]	= "B*_2^-";
      namemap[531]	= "B_s^0";
      namemap[-531]	= "B_s~^0";
      namemap[533]	= "B*_s^0";
      namemap[-533]	= "B*_s~^0";
      namemap[535]	= "B*_s2^0";
      namemap[-535]	= "B*_s2~^0";
      namemap[541]	= "B_c^+";
      namemap[-541]	= "B_c^-";
      namemap[543]	= "B*_c^+";
      namemap[-543]	= "B*_c^-";
      namemap[545]	= "B*_c2^+";
      namemap[-545]	= "B*_c2^-";
      namemap[551]	= "eta_b(1S)";
      namemap[553]	= "Upsilon(1S)";
      namemap[555]	= "chi_b2(1P)";
      namemap[557]	= "Upsilon_3(1D)";
      namemap[611]	= "T^+";
      namemap[-611]	= "T^-";
      namemap[613]	= "T*^+";
      namemap[-613]	= "T*^-";
      namemap[621]	= "T^0";
      namemap[-621]	= "T~^0";
      namemap[623]	= "T*^0";
      namemap[-623]	= "T*~^0";
      namemap[631]	= "T_s^+";
      namemap[-631]	= "T_s^-";
      namemap[633]	= "T*_s^+";
      namemap[-633]	= "T*_s^-";
      namemap[641]	= "T_c^0";
      namemap[-641]	= "T_c~^0";
      namemap[643]	= "T*_c^0";
      namemap[-643]	= "T*_c~^0";
      namemap[651]	= "T_b^+";
      namemap[-651]	= "T_b^-";
      namemap[653]	= "T*_b^+";
      namemap[-653]	= "T*_b^-";
      namemap[661]	= "eta_t";
      namemap[663]	= "theta";
      namemap[711]	= "L^0";
      namemap[-711]	= "L~^0";
      namemap[713]	= "L*^0";
      namemap[-713]	= "L*~^0";
      namemap[721]	= "L^-";
      namemap[-721]	= "L^+";
      namemap[723]	= "L*^-";
      namemap[-723]	= "L*^+";
      namemap[731]	= "L_s^0";
      namemap[-731]	= "L_s~^0";
      namemap[733]	= "L*_s^0";
      namemap[-733]	= "L*_s~^0";
      namemap[741]	= "L_c^-";
      namemap[-741]	= "L_c^+";
      namemap[743]	= "L*_c^-";
      namemap[-743]	= "L*_c^+";
      namemap[751]	= "L_b^0";
      namemap[-751]	= "L_b~^0";
      namemap[753]	= "L*_b^0";
      namemap[-753]	= "L*_b~^0";
      namemap[761]	= "L_t^-";
      namemap[-761]	= "L_t^+";
      namemap[763]	= "L*_t^-";
      namemap[-763]	= "L*_t^+";
      namemap[771]	= "eta_l";
      namemap[773]	= "theta_l";
      namemap[811]	= "H^+";
      namemap[-811]	= "H^-";
      namemap[813]	= "H*^+";
      namemap[-813]	= "H*^-";
      namemap[821]	= "H^0";
      namemap[-821]	= "H~^0";
      namemap[823]	= "H*^0";
      namemap[-823]	= "H*~^0";
      namemap[831]	= "H_s^+";
      namemap[-831]	= "H_s^-";
      namemap[833]	= "H*_s^+";
      namemap[-833]	= "H*_s^-";
      namemap[841]	= "H_c^0";
      namemap[-841]	= "H_c~^0";
      namemap[843]	= "H*_c^0";
      namemap[-843]	= "H*_c~^0";
      namemap[851]	= "H_b^+";
      namemap[-851]	= "H_b^-";
      namemap[853]	= "H*_b^+";
      namemap[-853]	= "H*_b^-";
      namemap[861]	= "H_t^0";
      namemap[-861]	= "H_t~^0";
      namemap[863]	= "H*_t^0";
      namemap[-863]	= "H*_t~^0";
      namemap[871]	= "H_l^+";
      namemap[-871]	= "H_l^-";
      namemap[873]	= "H*_l^+";
      namemap[-873]	= "H*_l^-";
      namemap[881]	= "eta_h";
      namemap[883]	= "theta_H";
      namemap[990]	= "pomeron";
      namemap[1103]	= "dd_1";
      namemap[-1103]	= "dd_1~";
      namemap[1112]	= "Delta(1620)^-";
      namemap[1114]	= "Delta^-";
      namemap[-1114]	= "Delta~^+";
      namemap[1116]	= "Delta(1905)^-";
      namemap[1118]	= "Delta(1950)^-";
      namemap[1212]	= "Delta(1620)^0";
      namemap[1214]	= "N(1520)^0";
      namemap[1216]	= "Delta(1905)^0";
      namemap[1218]	= "N(2190)^0";
      namemap[2101]	= "ud_0";
      namemap[-2101]	= "ud_0~";
      namemap[2103]	= "ud_1";
      namemap[-2103]	= "ud_1~";
      namemap[2112]	= "n^0";
      namemap[-2112]	= "n~^0";
      namemap[2114]	= "Delta^0";
      namemap[-2114]	= "Delta~^0";
      namemap[2116]	= "N(1675)^0";
      namemap[2118]	= "Delta(1950)^0";
      namemap[2122]	= "Delta(1620)^+";
      namemap[2124]	= "N(1520)^+";
      namemap[2126]	= "Delta(1905)^+";
      namemap[2128]	= "N(2190)^+";
      namemap[2203]	= "uu_1";
      namemap[-2203]	= "uu_1~";
      namemap[2212]	= "p^+";
      namemap[-2212]	= "p~^-";
      namemap[2214]	= "Delta^+";
      namemap[-2214]	= "Delta~^-";
      namemap[2216]	= "N(1675)^+";
      namemap[2218]	= "Delta(1950)^+";
      namemap[2222]	= "Delta(1620)^++";
      namemap[2224]	= "Delta^++";
      namemap[-2224]	= "Delta~^--";
      namemap[2226]	= "Delta(1905)^++";
      namemap[2228]	= "Delta(1950)^++";
      namemap[3101]	= "sd_0";
      namemap[-3101]	= "sd_0~";
      namemap[3103]	= "sd_1";
      namemap[-3103]	= "sd_1~";
      namemap[3112]	= "Sigma^-";
      namemap[-3112]	= "Sigma~^+";
      namemap[3114]	= "Sigma*^-";
      namemap[-3114]	= "Sigma*~^+";
      namemap[3116]	= "Sigma(1775)^-";
      namemap[-3116]	= "Sigma~(1775)^-";
      namemap[3118]	= "Sigma(2030)^-";
      namemap[-3118]	= "Sigma~(2030)^-";
      namemap[3122]	= "Lambda^0";
      namemap[-3122]	= "Lambda~^0";
      namemap[3124]	= "Lambda(1520)^0";
      namemap[-3124]	= "Lambda~(1520)^0";
      namemap[3126]	= "Lambda(1820)^0";
      namemap[-3126]	= "Lambda~(1820)^0";
      namemap[3128]	= "Lambda(2100)^0";
      namemap[-3128]	= "Lambda~(2100)^0";
      namemap[3201]	= "su_0";
      namemap[-3201]	= "su_0~";
      namemap[3203]	= "su_1";
      namemap[-3203]	= "su_1~";
      namemap[3212]	= "Sigma^0";
      namemap[-3212]	= "Sigma~^0";
      namemap[3214]	= "Sigma*^0";
      namemap[-3214]	= "Sigma*~^0";
      namemap[3216]	= "Sigma(1775)^0";
      namemap[-3216]	= "Sigma~(1775)^0";
      namemap[3218]	= "Sigma(2030)^0";
      namemap[-3218]	= "Sigma~(2030)^0";
      namemap[3222]	= "Sigma^+";
      namemap[-3222]	= "Sigma~^-";
      namemap[3224]	= "Sigma*^+";
      namemap[-3224]	= "Sigma*~^-";
      namemap[3226]	= "Sigma(1775)^+";
      namemap[-3226]	= "Sigma~(1775)^+";
      namemap[3228]	= "Sigma(2030)^+";
      namemap[-3228]	= "Sigma~(2030)^+";
      namemap[3303]	= "ss_1";
      namemap[-3303]	= "ss_1~";
      namemap[3312]	= "Xi^-";
      namemap[-3312]	= "Xi~^+";
      namemap[3314]	= "Xi*^-";
      namemap[-3314]	= "Xi*~^+";
      namemap[3322]	= "Xi^0";
      namemap[-3322]	= "Xi~^0";
      namemap[3324]	= "Xi*^0";
      namemap[-3324]	= "Xi*~^0";
      namemap[3334]	= "Omega^-";
      namemap[-3334]	= "Omega~^+";
      namemap[4101]	= "cd_0";
      namemap[-4101]	= "cd_0~";
      namemap[4103]	= "cd_1";
      namemap[-4103]	= "cd_1~";
      namemap[4112]	= "Sigma_c^0";
      namemap[-4112]	= "Sigma_c~^0";
      namemap[4114]	= "Sigma*_c^0";
      namemap[-4114]	= "Sigma*_c~^0";
      namemap[4122]	= "Lambda_c^+";
      namemap[-4122]	= "Lambda_c~^-";
      namemap[4132]	= "Xi_c^0";
      namemap[-4132]	= "Xi_c~^0";
      namemap[4201]	= "cu_0";
      namemap[-4201]	= "cu_0~";
      namemap[4203]	= "cu_1";
      namemap[-4203]	= "cu_1~";
      namemap[4212]	= "Sigma_c^+";
      namemap[-4212]	= "Sigma_c~^-";
      namemap[4214]	= "Sigma*_c^+";
      namemap[-4214]	= "Sigma*_c~^-";
      namemap[4222]	= "Sigma_c^++";
      namemap[-4222]	= "Sigma_c~^--";
      namemap[4224]	= "Sigma*_c^++";
      namemap[-4224]	= "Sigma*_c~^--";
      namemap[4232]	= "Xi_c^+";
      namemap[-4232]	= "Xi_c~^-";
      namemap[4301]	= "cs_0";
      namemap[-4301]	= "cs_0~";
      namemap[4303]	= "cs_1";
      namemap[-4303]	= "cs_1~";
      namemap[4312]	= "Xi'_c^0";
      namemap[-4312]	= "Xi'_c~^0";
      namemap[4314]	= "Xi*_c^0";
      namemap[-4314]	= "Xi*_c~^0";
      namemap[4322]	= "Xi'_c^+";
      namemap[-4322]	= "Xi'_c~^-";
      namemap[4324]	= "Xi*_c^+";
      namemap[-4324]	= "Xi*_c~^-";
      namemap[4332]	= "Omega_c^0";
      namemap[-4332]	= "Omega_c~^0";
      namemap[4334]	= "Omega*_c^0";
      namemap[-4334]	= "Omega*_c~^0";
      namemap[4403]	= "cc_1";
      namemap[-4403]	= "cc_1~";
      namemap[4412]	= "Xi_cc^+";
      namemap[-4412]	= "Xi_cc~^-";
      namemap[4414]	= "Xi*_cc^+";
      namemap[-4414]	= "Xi*_cc~^-";
      namemap[4422]	= "Xi_cc^++";
      namemap[-4422]	= "Xi_cc~^--";
      namemap[4424]	= "Xi*_cc^++";
      namemap[-4424]	= "Xi*_cc~^--";
      namemap[4432]	= "Omega_cc^+";
      namemap[-4432]	= "Omega_cc~^-";
      namemap[4434]	= "Omega*_cc^+";
      namemap[-4434]	= "Omega*_cc~^-";
      namemap[4444]	= "Omega*_ccc^++";
      namemap[-4444]	= "Omega*_ccc~^--";
      namemap[5101]	= "bd_0";
      namemap[-5101]	= "bd_0~";
      namemap[5103]	= "bd_1";
      namemap[-5103]	= "bd_1~";
      namemap[5112]	= "Sigma_b^-";
      namemap[-5112]	= "Sigma_b~^+";
      namemap[5114]	= "Sigma*_b^-";
      namemap[-5114]	= "Sigma*_b~^+";
      namemap[5122]	= "Lambda_b^0";
      namemap[-5122]	= "Lambda_b~^0";
      namemap[5132]	= "Xi_b^-";
      namemap[-5132]	= "Xi_b~^+";
      namemap[5142]	= "Xi_bc^0";
      namemap[-5142]	= "Xi_bc~^0";
      namemap[5201]	= "bu_0";
      namemap[-5201]	= "bu_0~";
      namemap[5203]	= "bu_1";
      namemap[-5203]	= "bu_1~";
      namemap[5212]	= "Sigma_b^0";
      namemap[-5212]	= "Sigma_b~^0";
      namemap[5214]	= "Sigma*_b^0";
      namemap[-5214]	= "Sigma*_b~^0";
      namemap[5222]	= "Sigma_b^+";
      namemap[-5222]	= "Sigma_b~^-";
      namemap[5224]	= "Sigma*_b^+";
      namemap[-5224]	= "Sigma*_b~^-";
      namemap[5232]	= "Xi_b^0";
      namemap[-5232]	= "Xi_b~^0";
      namemap[5242]	= "Xi_bc^+";
      namemap[-5242]	= "Xi_bc~^-";
      namemap[5301]	= "bs_0";
      namemap[-5301]	= "bs_0~";
      namemap[5303]	= "bs_1";
      namemap[-5303]	= "bs_1~";
      namemap[5312]	= "Xi'_b^-";
      namemap[-5312]	= "Xi'_b~^+";
      namemap[5314]	= "Xi*_b^-";
      namemap[-5314]	= "Xi*_b~^+";
      namemap[5322]	= "Xi'_b^0";
      namemap[-5322]	= "Xi'_b~^0";
      namemap[5324]	= "Xi*_b^0";
      namemap[-5324]	= "Xi*_b~^0";
      namemap[5332]	= "Omega_b^-";
      namemap[-5332]	= "Omega_b~^+";
      namemap[5334]	= "Omega*_b^-";
      namemap[-5334]	= "Omega*_b~^+";
      namemap[5342]	= "Omega_bc^0";
      namemap[-5342]	= "Omega_bc~^0";
      namemap[5401]	= "bc_0";
      namemap[-5401]	= "bc_0~";
      namemap[5403]	= "bc_1";
      namemap[-5403]	= "bc_1~";
      namemap[5412]	= "Xi'_bc^0";
      namemap[-5412]	= "Xi'_bc~^0";
      namemap[5414]	= "Xi*_bc^0";
      namemap[-5414]	= "Xi*_bc~^0";
      namemap[5422]	= "Xi'_bc^+";
      namemap[-5422]	= "Xi'_bc~^-";
      namemap[5424]	= "Xi*_bc^+";
      namemap[-5424]	= "Xi*_bc~^-";
      namemap[5432]	= "Omega'_bc^0";
      namemap[-5432]	= "Omega'_bc~^0";
      namemap[5434]	= "Omega*_bc^0";
      namemap[-5434]	= "Omega*_bc~^0";
      namemap[5442]	= "Omega_bcc^+";
      namemap[-5442]	= "Omega_bcc~^-";
      namemap[5444]	= "Omega*_bcc^+";
      namemap[-5444]	= "Omega*_bcc~^-";
      namemap[5503]	= "bb_1";
      namemap[-5503]	= "bb_1~";
      namemap[5512]	= "Xi_bb^-";
      namemap[-5512]	= "Xi_bb~^+";
      namemap[5514]	= "Xi*_bb^-";
      namemap[-5514]	= "Xi*_bb~^+";
      namemap[5522]	= "Xi_bb^0";
      namemap[-5522]	= "Xi_bb~^0";
      namemap[5524]	= "Xi*_bb^0";
      namemap[-5524]	= "Xi*_bb~^0";
      namemap[5532]	= "Omega_bb^-";
      namemap[-5532]	= "Omega_bb~^+";
      namemap[5534]	= "Omega*_bb^-";
      namemap[-5534]	= "Omega*_bb~^+";
      namemap[5542]	= "Omega_bbc^0";
      namemap[-5542]	= "Omega_bbc~^0";
      namemap[5544]	= "Omega*_bbc^0";
      namemap[-5544]	= "Omega*_bbc~^0";
      namemap[5554]	= "Omega*_bbb^-";
      namemap[-5554]	= "Omega*_bbb~^+";
      namemap[6101]	= "td_0";
      namemap[-6101]	= "td_0~";
      namemap[6103]	= "td_1";
      namemap[-6103]	= "td_1~";
      namemap[6112]	= "Sigma_t^0";
      namemap[-6112]	= "Sigma_t~^0";
      namemap[6114]	= "Sigma*_t^0";
      namemap[-6114]	= "Sigma*_t~^0";
      namemap[6122]	= "Lambda_t^+";
      namemap[-6122]	= "Lambda_t~^-";
      namemap[6132]	= "Xi_t^0";
      namemap[-6132]	= "Xi_t~^0";
      namemap[6142]	= "Xi_tc^+";
      namemap[-6142]	= "Xi_tc~^-";
      namemap[6152]	= "Xi_tb^0";
      namemap[-6152]	= "Xi_tb~^0";
      namemap[6201]	= "tu_0";
      namemap[-6201]	= "tu_0~";
      namemap[6203]	= "tu_1";
      namemap[-6203]	= "tu_1~";
      namemap[6212]	= "Sigma_t^+";
      namemap[-6212]	= "Sigma_t~^-";
      namemap[6214]	= "Sigma*_t^+";
      namemap[-6214]	= "Sigma*_t~^-";
      namemap[6222]	= "Sigma_t^++";
      namemap[-6222]	= "Sigma_t~^--";
      namemap[6224]	= "Sigma*_t^++";
      namemap[-6224]	= "Sigma*_t~^--";
      namemap[6232]	= "Xi_t^+";
      namemap[-6232]	= "Xi_t~^-";
      namemap[6242]	= "Xi_tc^++";
      namemap[-6242]	= "Xi_tc~^--";
      namemap[6252]	= "Xi_tb^+";
      namemap[-6252]	= "Xi_tb~^-";
      namemap[6301]	= "ts_0";
      namemap[-6301]	= "ts_0~";
      namemap[6303]	= "ts_1";
      namemap[-6303]	= "ts_1~";
      namemap[6312]	= "Xi'_t^0";
      namemap[-6312]	= "Xi'_t~^0";
      namemap[6314]	= "Xi*_t^0";
      namemap[-6314]	= "Xi*_t~^0";
      namemap[6322]	= "Xi'_t^+";
      namemap[-6322]	= "Xi'_t~^-";
      namemap[6324]	= "Xi*_t^+";
      namemap[-6324]	= "Xi*_t~^-";
      namemap[6332]	= "Omega_t^0";
      namemap[-6332]	= "Omega_t~^0";
      namemap[6334]	= "Omega*_t^0";
      namemap[-6334]	= "Omega*_t~^0";
      namemap[6342]	= "Omega_tc^+";
      namemap[-6342]	= "Omega_tc~^-";
      namemap[6352]	= "Omega_tb^0";
      namemap[-6352]	= "Omega_tb~^0";
      namemap[6401]	= "tc_0";
      namemap[-6401]	= "tc_0~";
      namemap[6403]	= "tc_1";
      namemap[-6403]	= "tc_1~";
      namemap[6412]	= "Xi'_tc^+";
      namemap[-6412]	= "Xi'_tc~^-";
      namemap[6414]	= "Xi*_tc^+";
      namemap[-6414]	= "Xi*_tc~^-";
      namemap[6422]	= "Xi'_tc^++";
      namemap[-6422]	= "Xi'_tc~^--";
      namemap[6424]	= "Xi*_tc^++";
      namemap[-6424]	= "Xi*_tc~^--";
      namemap[6432]	= "Omega'_tc^+";
      namemap[-6432]	= "Omega'_tc~^-";
      namemap[6434]	= "Omega*_tc^+";
      namemap[-6434]	= "Omega*_tc~^-";
      namemap[6442]	= "Omega_tcc^++";
      namemap[-6442]	= "Omega_tcc~^--";
      namemap[6444]	= "Omega*_tcc^++";
      namemap[-6444]	= "Omega*_tcc~^--";
      namemap[6452]	= "Omega_tbc^+";
      namemap[-6452]	= "Omega_tbc~^-";
      namemap[6501]	= "tb_0";
      namemap[-6501]	= "tb_0~";
      namemap[6503]	= "tb_1";
      namemap[-6503]	= "tb_1~";
      namemap[6512]	= "Xi'_tb^0";
      namemap[-6512]	= "Xi'_tb~^0";
      namemap[6514]	= "Xi*_tb^0";
      namemap[-6514]	= "Xi*_tb~^0";
      namemap[6522]	= "Xi'_tb^+";
      namemap[-6522]	= "Xi'_tb~^-";
      namemap[6524]	= "Xi*_tb^+";
      namemap[-6524]	= "Xi*_tb~^-";
      namemap[6532]	= "Omega'_tb^0";
      namemap[-6532]	= "Omega'_tb~^0";
      namemap[6534]	= "Omega*_tb^0";
      namemap[-6534]	= "Omega*_tb~^0";
      namemap[6542]	= "Omega'_tbc^+";
      namemap[-6542]	= "Omega'_tbc~^-";
      namemap[6544]	= "Omega*_tbc^+";
      namemap[-6544]	= "Omega*_tbc~^-";
      namemap[6552]	= "Omega_tbb^0";
      namemap[-6552]	= "Omega_tbb~^0";
      namemap[6554]	= "Omega*_tbb^0";
      namemap[-6554]	= "Omega*_tbb~^0";
      namemap[6603]	= "tt_1";
      namemap[-6603]	= "tt_1~";
      namemap[6612]	= "Xi_tt^+";
      namemap[-6612]	= "Xi_tt~^-";
      namemap[6614]	= "Xi*_tt^+";
      namemap[-6614]	= "Xi*_tt~^-";
      namemap[6622]	= "Xi_tt^++";
      namemap[-6622]	= "Xi_tt~^--";
      namemap[6624]	= "Xi*_tt^++";
      namemap[-6624]	= "Xi*_tt~^--";
      namemap[6632]	= "Omega_tt^+";
      namemap[-6632]	= "Omega_tt~^-";
      namemap[6634]	= "Omega*_tt^+";
      namemap[-6634]	= "Omega*_tt~^-";
      namemap[6642]	= "Omega_ttc^++";
      namemap[-6642]	= "Omega_ttc~^--";
      namemap[6644]	= "Omega*_ttc^++";
      namemap[-6644]	= "Omega*_ttc~^--";
      namemap[6652]	= "Omega_ttb^+";
      namemap[-6652]	= "Omega_ttb~^-";
      namemap[6654]	= "Omega*_ttb^+";
      namemap[-6654]	= "Omega*_ttb~^-";
      namemap[6664]	= "Omega*_ttt^++";
      namemap[-6664]	= "Omega*_ttt~^--";
      namemap[7101]	= "b'd_0";
      namemap[-7101]	= "b'd_0~";
      namemap[7103]	= "b'd_1";
      namemap[-7103]	= "b'd_1~";
      namemap[7112]	= "Sigma_b'^-";
      namemap[-7112]	= "Sigma_b'~^+";
      namemap[7114]	= "Sigma*_b'^-";
      namemap[-7114]	= "Sigma*_b'~^+";
      namemap[7122]	= "Lambda_b'^0";
      namemap[-7122]	= "Lambda_b'~^0";
      namemap[7132]	= "Xi_b'^-";
      namemap[-7132]	= "Xi_b'~^+";
      namemap[7142]	= "Xi_b'c^0";
      namemap[-7142]	= "Xi_b'c~^0";
      namemap[7152]	= "Xi_b'b^-";
      namemap[-7152]	= "Xi_b'b~^+";
      namemap[7162]	= "Xi_b't^0";
      namemap[-7162]	= "Xi_b't~^0";
      namemap[7201]	= "b'u_0";
      namemap[-7201]	= "b'u_0~";
      namemap[7203]	= "b'u_1";
      namemap[-7203]	= "b'u_1~";
      namemap[7212]	= "Sigma_b'^0";
      namemap[-7212]	= "Sigma_b'~^0";
      namemap[7214]	= "Sigma*_b'^0";
      namemap[-7214]	= "Sigma*_b'~^0";
      namemap[7222]	= "Sigma_b'^+";
      namemap[-7222]	= "Sigma_b'~^-";
      namemap[7224]	= "Sigma*_b'^+";
      namemap[-7224]	= "Sigma*_b'~^-";
      namemap[7232]	= "Xi_b'^0";
      namemap[-7232]	= "Xi_b'~^0";
      namemap[7242]	= "Xi_b'c^+";
      namemap[-7242]	= "Xi_b'c~^-";
      namemap[7252]	= "Xi_b'b^0";
      namemap[-7252]	= "Xi_b'b~^0";
      namemap[7262]	= "Xi_b't^+";
      namemap[-7262]	= "Xi_b't~^-";
      namemap[7301]	= "b's_0";
      namemap[-7301]	= "b's_0~";
      namemap[7303]	= "b's_1";
      namemap[-7303]	= "b's_1~";
      namemap[7312]	= "Xi'_b'^-";
      namemap[-7312]	= "Xi'_b'~^+";
      namemap[7314]	= "Xi*_b'^-";
      namemap[-7314]	= "Xi*_b'~^+";
      namemap[7322]	= "Xi'_b'^0";
      namemap[-7322]	= "Xi'_b'~^0";
      namemap[7324]	= "Xi*_b'^0";
      namemap[-7324]	= "Xi*_b'~^0";
      namemap[7332]	= "Omega'_b'^-";
      namemap[-7332]	= "Omega'_b'~^+";
      namemap[7334]	= "Omega*_b'^-";
      namemap[-7334]	= "Omega*_b'~^+";
      namemap[7342]	= "Omega_b'c^0";
      namemap[-7342]	= "Omega_b'c~^0";
      namemap[7352]	= "Omega_b'b^-";
      namemap[-7352]	= "Omega_b'b~^+";
      namemap[7362]	= "Omega_b't^0";
      namemap[-7362]	= "Omega_b't~^0";
      namemap[7401]	= "b'c_0";
      namemap[-7401]	= "b'c_0~";
      namemap[7403]	= "b'c_1";
      namemap[-7403]	= "b'c_1~";
      namemap[7412]	= "Xi'_b'c^0";
      namemap[-7412]	= "Xi'_b'c~^0";
      namemap[7414]	= "Xi*_b'c^0";
      namemap[-7414]	= "Xi*_b'c~^0";
      namemap[7422]	= "Xi'_b'c^+";
      namemap[-7422]	= "Xi'_b'c~^-";
      namemap[7424]	= "Xi*_b'c^+";
      namemap[-7424]	= "Xi*_b'c~^-";
      namemap[7432]	= "Omega'_b'c^0";
      namemap[-7432]	= "Omega'_b'c~^0";
      namemap[7434]	= "Omega*_b'c^0";
      namemap[-7434]	= "Omega*_b'c~^0";
      namemap[7442]	= "Omega'_b'cc^+";
      namemap[-7442]	= "Omega'_b'cc~^-";
      namemap[7444]	= "Omega*_b'cc^+";
      namemap[-7444]	= "Omega*_b'cc~^-";
      namemap[7452]	= "Omega_b'bc^0";
      namemap[-7452]	= "Omega_b'bc~^0";
      namemap[7462]	= "Omega_b'tc^+";
      namemap[-7462]	= "Omega_b'tc~^-";
      namemap[7501]	= "b'b_0";
      namemap[-7501]	= "b'b_0~";
      namemap[7503]	= "b'b_1";
      namemap[-7503]	= "b'b_1~";
      namemap[7512]	= "Xi'_b'b^-";
      namemap[-7512]	= "Xi'_b'b~^+";
      namemap[7514]	= "Xi*_b'b^-";
      namemap[-7514]	= "Xi*_b'b~^+";
      namemap[7522]	= "Xi'_b'b^0";
      namemap[-7522]	= "Xi'_b'b~^0";
      namemap[7524]	= "Xi*_b'b^0";
      namemap[-7524]	= "Xi*_b'b~^0";
      namemap[7532]	= "Omega'_b'b^-";
      namemap[-7532]	= "Omega'_b'b~^+";
      namemap[7534]	= "Omega*_b'b^-";
      namemap[-7534]	= "Omega*_b'b~^+";
      namemap[7542]	= "Omega'_b'bc^0";
      namemap[-7542]	= "Omega'_b'bc~^0";
      namemap[7544]	= "Omega*_b'bc^0";
      namemap[-7544]	= "Omega*_b'bc~^0";
      namemap[7552]	= "Omega'_b'bb^-";
      namemap[-7552]	= "Omega'_b'bb~^+";
      namemap[7554]	= "Omega*_b'bb^-";
      namemap[-7554]	= "Omega*_b'bb~^+";
      namemap[7562]	= "Omega_b'tb^0";
      namemap[-7562]	= "Omega_b'tb~^0";
      namemap[7601]	= "b't_0";
      namemap[-7601]	= "b't_0~";
      namemap[7603]	= "b't_1";
      namemap[-7603]	= "b't_1~";
      namemap[7612]	= "Xi'_b't^0";
      namemap[-7612]	= "Xi'_b't~^0";
      namemap[7614]	= "Xi*_b't^0";
      namemap[-7614]	= "Xi*_b't~^0";
      namemap[7622]	= "Xi'_b't^+";
      namemap[-7622]	= "Xi'_b't~^-";
      namemap[7624]	= "Xi*_b't^+";
      namemap[-7624]	= "Xi*_b't~^-";
      namemap[7632]	= "Omega'_b't^0";
      namemap[-7632]	= "Omega'_b't~^0";
      namemap[7634]	= "Omega*_b't^0";
      namemap[-7634]	= "Omega*_b't~^0";
      namemap[7642]	= "Omega'_b'tc^+";
      namemap[-7642]	= "Omega'_b'tc~^-";
      namemap[7644]	= "Omega*_b'tc^+";
      namemap[-7644]	= "Omega*_b'tc~^-";
      namemap[7652]	= "Omega'_b'tb^0";
      namemap[-7652]	= "Omega'_b'tb~^0";
      namemap[7654]	= "Omega*_b'tb^0";
      namemap[-7654]	= "Omega*_b'tb~^0";
      namemap[7662]	= "Omega'_b'tt^+";
      namemap[-7662]	= "Omega'_b'tt~^-";
      namemap[7664]	= "Omega*_b'tt^+";
      namemap[-7664]	= "Omega*_b'tt~^-";
      namemap[7703]	= "b'b'_1";
      namemap[-7703]	= "b'b'_1~";
      namemap[7712]	= "Xi'_b'b'^-";
      namemap[-7712]	= "Xi'_b'b'~^+";
      namemap[7714]	= "Xi*_b'b'^-";
      namemap[-7714]	= "Xi*_b'b'~^+";
      namemap[7722]	= "Xi'_b'b'^0";
      namemap[-7722]	= "Xi'_b'b'~^0";
      namemap[7724]	= "Xi*_b'b'^0";
      namemap[-7724]	= "Xi*_b'b'~^0";
      namemap[7732]	= "Omega'_b'b'^-";
      namemap[-7732]	= "Omega'_b'b'~^+";
      namemap[7734]	= "Omega*_b'b'^-";
      namemap[-7734]	= "Omega*_b'b'~^+";
      namemap[7742]	= "Omega'_b'b'c^0";
      namemap[-7742]	= "Omega'_b'b'c~^0";
      namemap[7744]	= "Omega*_b'b'c^0";
      namemap[-7744]	= "Omega*_b'b'c~^0";
      namemap[7752]	= "Omega'_b'b'b^-";
      namemap[-7752]	= "Omega'_b'b'b~^+";
      namemap[7754]	= "Omega*_b'b'b^-";
      namemap[-7754]	= "Omega*_b'b'b~^+";
      namemap[7762]	= "Omega'_b'b't^0";
      namemap[-7762]	= "Omega'_b'b't~^0";
      namemap[7764]	= "Omega*_b'b't^0";
      namemap[-7764]	= "Omega*_b'b't~^0";
      namemap[7774]	= "Omega*_b'b'b'^-";
      namemap[-7774]	= "Omega*_b'b'b'~^+";
      namemap[8101]	= "t'd_0";
      namemap[-8101]	= "t'd_0~";
      namemap[8103]	= "t'd_1";
      namemap[-8103]	= "t'd_1~";
      namemap[8112]	= "Sigma_t'^0";
      namemap[-8112]	= "Sigma_t'~^0";
      namemap[8114]	= "Sigma*_t'^0";
      namemap[-8114]	= "Sigma*_t'~^0";
      namemap[8122]	= "Lambda_t'^+";
      namemap[-8122]	= "Lambda_t'~^-";
      namemap[8132]	= "Xi_t'^0";
      namemap[-8132]	= "Xi_t'~^0";
      namemap[8142]	= "Xi_t'c^+";
      namemap[-8142]	= "Xi_t'c~^-";
      namemap[8152]	= "Xi_t'b^0";
      namemap[-8152]	= "Xi_t'b~^0";
      namemap[8162]	= "Xi_t't^+";
      namemap[-8162]	= "Xi_t't~^-";
      namemap[8172]	= "Xi_t'b'^0";
      namemap[-8172]	= "Xi_t'b'~^0";
      namemap[8201]	= "t'u_0";
      namemap[-8201]	= "t'u_0~";
      namemap[8203]	= "t'u_1";
      namemap[-8203]	= "t'u_1~";
      namemap[8212]	= "Sigma_t'^+";
      namemap[-8212]	= "Sigma_t'~^-";
      namemap[8214]	= "Sigma*_t'^+";
      namemap[-8214]	= "Sigma*_t'~^-";
      namemap[8222]	= "Sigma_t'^++";
      namemap[-8222]	= "Sigma_t'~^--";
      namemap[8224]	= "Sigma*_t'^++";
      namemap[-8224]	= "Sigma*_t'~^--";
      namemap[8232]	= "Xi_t'^+";
      namemap[-8232]	= "Xi_t'~^-";
      namemap[8242]	= "Xi_t'c^++";
      namemap[-8242]	= "Xi_t'c~^--";
      namemap[8252]	= "Xi_t'b^+";
      namemap[-8252]	= "Xi_t'b~^-";
      namemap[8262]	= "Xi_t't^++";
      namemap[-8262]	= "Xi_t't~^--";
      namemap[8272]	= "Xi_t'b'^+";
      namemap[-8272]	= "Xi_t'b'~^-";
      namemap[8301]	= "t's_0";
      namemap[-8301]	= "t's_0~";
      namemap[8303]	= "t's_1";
      namemap[-8303]	= "t's_1~";
      namemap[8312]	= "Xi'_t'^0";
      namemap[-8312]	= "Xi'_t'~^0";
      namemap[8314]	= "Xi*_t'^0";
      namemap[-8314]	= "Xi*_t'~^0";
      namemap[8322]	= "Xi'_t'^+";
      namemap[-8322]	= "Xi'_t'~^-";
      namemap[8324]	= "Xi*_t'^+";
      namemap[-8324]	= "Xi*_t'~^-";
      namemap[8332]	= "Omega'_t'^0";
      namemap[-8332]	= "Omega'_t'~^0";
      namemap[8334]	= "Omega*_t'^0";
      namemap[-8334]	= "Omega*_t'~^0";
      namemap[8342]	= "Omega_t'c^+";
      namemap[-8342]	= "Omega_t'c~^-";
      namemap[8352]	= "Omega_t'b^0";
      namemap[-8352]	= "Omega_t'b~^0";
      namemap[8362]	= "Omega_t't^+";
      namemap[-8362]	= "Omega_t't~^-";
      namemap[8372]	= "Omega_t'b'^0";
      namemap[-8372]	= "Omega_t'b'~^0";
      namemap[8401]	= "t'c_0";
      namemap[-8401]	= "t'c_0~";
      namemap[8403]	= "t'c_1";
      namemap[-8403]	= "t'c_1~";
      namemap[8412]	= "Xi'_t'c^+";
      namemap[-8412]	= "Xi'_t'c~^-";
      namemap[8414]	= "Xi*_t'c^+";
      namemap[-8414]	= "Xi*_t'c~^-";
      namemap[8422]	= "Xi'_t'c^++";
      namemap[-8422]	= "Xi'_t'c~^--";
      namemap[8424]	= "Xi*_t'c^++";
      namemap[-8424]	= "Xi*_t'c~^--";
      namemap[8432]	= "Omega'_t'c^+";
      namemap[-8432]	= "Omega'_t'c~^-";
      namemap[8434]	= "Omega*_t'c^+";
      namemap[-8434]	= "Omega*_t'c~^-";
      namemap[8442]	= "Omega'_t'cc^++";
      namemap[-8442]	= "Omega'_t'cc~^--";
      namemap[8444]	= "Omega*_t'cc^++";
      namemap[-8444]	= "Omega*_t'cc~^--";
      namemap[8452]	= "Omega_t'bc^+";
      namemap[-8452]	= "Omega_t'bc~^-";
      namemap[8462]	= "Omega_t'tc^++";
      namemap[-8462]	= "Omega_t'tc~^--";
      namemap[8472]	= "Omega_t'b'c ^+";
      namemap[-8472]	= "Omega_t'b'c ~^-";
      namemap[8501]	= "t'b_0";
      namemap[-8501]	= "t'b_0~";
      namemap[8503]	= "t'b_1";
      namemap[-8503]	= "t'b_1~";
      namemap[8512]	= "Xi'_t'b^0";
      namemap[-8512]	= "Xi'_t'b~^0";
      namemap[8514]	= "Xi*_t'b^0";
      namemap[-8514]	= "Xi*_t'b~^0";
      namemap[8522]	= "Xi'_t'b^+";
      namemap[-8522]	= "Xi'_t'b~^-";
      namemap[8524]	= "Xi*_t'b^+";
      namemap[-8524]	= "Xi*_t'b~^-";
      namemap[8532]	= "Omega'_t'b^0";
      namemap[-8532]	= "Omega'_t'b~^0";
      namemap[8534]	= "Omega*_t'b^0";
      namemap[-8534]	= "Omega*_t'b~^0";
      namemap[8542]	= "Omega'_t'bc^+";
      namemap[-8542]	= "Omega'_t'bc~^-";
      namemap[8544]	= "Omega*_t'bc^+";
      namemap[-8544]	= "Omega*_t'bc~^-";
      namemap[8552]	= "Omega'_t'bb^0";
      namemap[-8552]	= "Omega'_t'bb~^0";
      namemap[8554]	= "Omega*_t'bb^0";
      namemap[-8554]	= "Omega*_t'bb~^0";
      namemap[8562]	= "Omega_t'tb^+";
      namemap[-8562]	= "Omega_t'tb~^-";
      namemap[8572]	= "Omega_t'b'b ^0";
      namemap[-8572]	= "Omega_t'b'b ~^0";
      namemap[8601]	= "t't_0";
      namemap[-8601]	= "t't_0~";
      namemap[8603]	= "t't_1";
      namemap[-8603]	= "t't_1~";
      namemap[8612]	= "Xi'_t't^+";
      namemap[-8612]	= "Xi'_t't~^-";
      namemap[8614]	= "Xi*_t't^+";
      namemap[-8614]	= "Xi*_t't~^-";
      namemap[8622]	= "Xi'_t't^++";
      namemap[-8622]	= "Xi'_t't~^--";
      namemap[8624]	= "Xi*_t't^++";
      namemap[-8624]	= "Xi*_t't~^--";
      namemap[8632]	= "Omega'_t't^+";
      namemap[-8632]	= "Omega'_t't~^-";
      namemap[8634]	= "Omega*_t't^+";
      namemap[-8634]	= "Omega*_t't~^-";
      namemap[8642]	= "Omega'_t'tc^++";
      namemap[-8642]	= "Omega'_t'tc~^--";
      namemap[8644]	= "Omega*_t'tc^++";
      namemap[-8644]	= "Omega*_t'tc~^--";
      namemap[8652]	= "Omega'_t'tb^+";
      namemap[-8652]	= "Omega'_t'tb~^-";
      namemap[8654]	= "Omega*_t'tb^+";
      namemap[-8654]	= "Omega*_t'tb~^-";
      namemap[8662]	= "Omega'_t'tt^++";
      namemap[-8662]	= "Omega'_t'tt~^--";
      namemap[8664]	= "Omega*_t'tt^++";
      namemap[-8664]	= "Omega*_t'tt~^--";
      namemap[8672]	= "Omega_t'b't ^+";
      namemap[-8672]	= "Omega_t'b't ~^-";
      namemap[8701]	= "t'b'_0";
      namemap[-8701]	= "t'b'_0~";
      namemap[8703]	= "t'b'_1";
      namemap[-8703]	= "t'b'_1~";
      namemap[8712]	= "Xi'_t'b'^0";
      namemap[-8712]	= "Xi'_t'b'~^0";
      namemap[8714]	= "Xi*_t'b'^0";
      namemap[-8714]	= "Xi*_t'b'~^0";
      namemap[8722]	= "Xi'_t'b'^+";
      namemap[-8722]	= "Xi'_t'b'~^-";
      namemap[8724]	= "Xi*_t'b'^+";
      namemap[-8724]	= "Xi*_t'b'~^-";
      namemap[8732]	= "Omega'_t'b'^0";
      namemap[-8732]	= "Omega'_t'b'~^0";
      namemap[8734]	= "Omega*_t'b'^0";
      namemap[-8734]	= "Omega*_t'b'~^0";
      namemap[8742]	= "Omega'_t'b'c^+";
      namemap[-8742]	= "Omega'_t'b'c~^-";
      namemap[8744]	= "Omega*_t'b'c^+";
      namemap[-8744]	= "Omega*_t'b'c~^-";
      namemap[8752]	= "Omega'_t'b'b^0";
      namemap[-8752]	= "Omega'_t'b'b~^0";
      namemap[8754]	= "Omega*_t'b'b^0";
      namemap[-8754]	= "Omega*_t'b'b~^0";
      namemap[8762]	= "Omega'_t'b't^+";
      namemap[-8762]	= "Omega'_t'b't~^-";
      namemap[8764]	= "Omega*_t'b't^+";
      namemap[-8764]	= "Omega*_t'b't~^-";
      namemap[8772]	= "Omega'_t'b'b'^0";
      namemap[-8772]	= "Omega'_t'b'b'~^0";
      namemap[8774]	= "Omega*_t'b'b'^0";
      namemap[-8774]	= "Omega*_t'b'b'~^0";
      namemap[8803]	= "t't'_1";
      namemap[-8803]	= "t't'_1~";
      namemap[8812]	= "Xi'_t't'^+";
      namemap[-8812]	= "Xi'_t't'~^-";
      namemap[8814]	= "Xi*_t't'^+";
      namemap[-8814]	= "Xi*_t't'~^-";
      namemap[8822]	= "Xi'_t't'^++";
      namemap[-8822]	= "Xi'_t't'~^--";
      namemap[8824]	= "Xi*_t't'^++";
      namemap[-8824]	= "Xi*_t't'~^--";
      namemap[8832]	= "Omega'_t't'^+";
      namemap[-8832]	= "Omega'_t't'~^-";
      namemap[8834]	= "Omega*_t't'^+";
      namemap[-8834]	= "Omega*_t't'~^-";
      namemap[8842]	= "Omega'_t't'c^++";
      namemap[-8842]	= "Omega'_t't'c~^--";
      namemap[8844]	= "Omega*_t't'c^++";
      namemap[-8844]	= "Omega*_t't'c~^--";
      namemap[8852]	= "Omega'_t't'b^+";
      namemap[-8852]	= "Omega'_t't'b~^-";
      namemap[8854]	= "Omega*_t't'b^+";
      namemap[-8854]	= "Omega*_t't'b~^-";
      namemap[8862]	= "Omega'_t't't^++";
      namemap[-8862]	= "Omega'_t't't~^--";
      namemap[8864]	= "Omega*_t't't^++";
      namemap[-8864]	= "Omega*_t't't~^--";
      namemap[8872]	= "Omega'_t't'b'^+";
      namemap[-8872]	= "Omega'_t't'b'~^-";
      namemap[8874]	= "Omega*_t't'b'^+";
      namemap[-8874]	= "Omega*_t't'b'~^-";
      namemap[8884]	= "Omega*_t't't'^++";
      namemap[-8884]	= "Omega*_t't't'~^--";
      namemap[9990]	= "odderon";
      namemap[10022]	= "virtual-photon";
      namemap[10111]	= "a_0(1450)^0";
      namemap[10113]	= "b_1(1235)^0";
      namemap[10115]	= "pi_2(1670)^0";
      namemap[10211]	= "a_0(1450)^+";
      namemap[-10211]	= "a_0(1450)^-";
      namemap[10213]	= "b_1(1235)^+";
      namemap[-10213]	= "b_1(1235)^-";
      namemap[10215]	= "pi_2(1670)^+";
      namemap[-10215]	= "pi_2(1670)^-";
      namemap[10221]	= "f_0(1370)";
      namemap[10223]	= "h_1(1170)";
      namemap[10225]	= "eta_2(1645)";
      namemap[10311]	= "K*_0(1430)^0";
      namemap[-10311]	= "K*_0(1430)~^0";
      namemap[10313]	= "K_1(1270)^0";
      namemap[-10313]	= "K_1(1270)~^0";
      namemap[10315]	= "K_2(1770)^0";
      namemap[-10315]	= "K_2(1770)~^0";
      namemap[10321]	= "K*_0(1430)^+";
      namemap[-10321]	= "K*_0(1430)^-";
      namemap[10323]	= "K_1(1270)^+";
      namemap[-10323]	= "K_1(1270)^-";
      namemap[10325]	= "K_2(1770)^+";
      namemap[-10325]	= "K_2(1770)^-";
      namemap[10331]	= "f_0(1710)";
      namemap[10333]	= "h_1(1380)";
      namemap[10335]	= "eta_2(1870)";
      namemap[10411]	= "D*_0(2400)^+";
      namemap[-10411]	= "D*_0(2400)^-";
      namemap[10413]	= "D_1(2420)^+";
      namemap[-10413]	= "D_1(2420)^-";
      namemap[10421]	= "D*_0(2400)^0";
      namemap[-10421]	= "D*_0(2400)~^0";
      namemap[10423]	= "D_1(2420)^0";
      namemap[-10423]	= "D_1(2420)~^0";
      namemap[10431]	= "D*_s0(2317)^+";
      namemap[-10431]	= "D*_s0(2317)^-";
      namemap[10433]	= "D_s1(2536)^+";
      namemap[-10433]	= "D_s1(2536)^-";
      namemap[10441]	= "chi_c0(1P)";
      namemap[10443]	= "hc(1P)";
      namemap[10511]	= "B*_0^0";
      namemap[-10511]	= "B*_0~^0";
      namemap[10513]	= "B_1(L)^0";
      namemap[-10513]	= "B_1(L)~^0";
      namemap[10521]	= "B*_0^+";
      namemap[-10521]	= "B*_0^-";
      namemap[10523]	= "B_1(L)^+";
      namemap[-10523]	= "B_1(L)^-";
      namemap[10531]	= "B*_s0^0";
      namemap[-10531]	= "B*_s0~^0";
      namemap[10533]	= "B_s1(L)^0";
      namemap[-10533]	= "B_s1(L)~^0";
      namemap[10541]	= "B*_c0^+";
      namemap[-10541]	= "B*_c0^-";
      namemap[10543]	= "B_c1(L)^+";
      namemap[-10543]	= "B_c1(L)^-";
      namemap[10551]	= "chi_b0(1P)";
      namemap[10553]	= "h_b(1P)";
      namemap[10555]	= "eta_b2(1D)";
      namemap[11114]	= "Delta(1700)^-";
      namemap[11116]	= "Delta(1930)^-";
      namemap[11216]	= "Delta(1930)^0";
      namemap[12112]	= "N(1440)^0";
      namemap[12114]	= "Delta(1700)^0";
      namemap[12116]	= "N(1680)^0";
      namemap[12126]	= "Delta(1930)^+";
      namemap[12212]	= "N(1440)^+";
      namemap[12214]	= "Delta(1700)^+";
      namemap[12216]	= "N(1680)^+";
      namemap[12224]	= "Delta(1700)^++";
      namemap[12226]	= "Delta(1930)^++";
      namemap[13112]	= "Sigma(1660)^-";
      namemap[-13112]	= "Sigma~(1660)^-";
      namemap[13114]	= "Sigma(1670)^-";
      namemap[-13114]	= "Sigma~(1670)^-";
      namemap[13116]	= "Sigma(1915)^-";
      namemap[-13116]	= "Sigma~(1915)^-";
      namemap[13122]	= "Lambda(1405)^0";
      namemap[-13122]	= "Lambda~(1405)^0";
      namemap[13124]	= "Lambda(1690)^0";
      namemap[-13124]	= "Lambda~(1690)^0";
      namemap[13126]	= "Lambda(1830)^0";
      namemap[-13126]	= "Lambda~(1830)^0";
      namemap[13212]	= "Sigma(1660)^0";
      namemap[-13212]	= "Sigma~(1660)^0";
      namemap[13214]	= "Sigma(1670)^0";
      namemap[-13214]	= "Sigma~(1670)^0";
      namemap[13216]	= "Sigma(1915)^0";
      namemap[-13216]	= "Sigma~(1915)^0";
      namemap[13222]	= "Sigma(1660)^+";
      namemap[-13222]	= "Sigma~(1660)^+";
      namemap[13224]	= "Sigma(1670)^+";
      namemap[-13224]	= "Sigma~(1670)^+";
      namemap[13226]	= "Sigma(1915)^+";
      namemap[-13226]	= "Sigma~(1915)^+";
      namemap[13314]	= "Xi(1820)^-";
      namemap[-13314]	= "Xi(1820)~^+";
      namemap[13324]	= "Xi(1820)^0";
      namemap[-13324]	= "Xi(1820)~^0";
      namemap[14122]	= "Lambda_c(2593)^+";
      namemap[-14122]	= "Lambda_c~(2593)^-";
      namemap[14124]	= "Lambda_c(2625)^+";
      namemap[-14124]	= "Lambda_c~(2625)^-";
      namemap[20022]	= "Cerenkov-radiation";
      namemap[20113]	= "a_1(1260)^0";
      namemap[20213]	= "a_1(1260)^+";
      namemap[-20213]	= "a_1(1260)^-";
      namemap[20223]	= "f_1(1285)";
      namemap[20313]	= "K_1(1400)^0";
      namemap[-20313]	= "K_1(1400)~^0";
      namemap[20315]	= "K_2(1820)^0";
      namemap[-20315]	= "K_2(1820)~^0";
      namemap[20323]	= "K_1(1400)^+";
      namemap[-20323]	= "K_1(1400)^-";
      namemap[20325]	= "K_2(1820)^+";
      namemap[-20325]	= "K_2(1820)^-";
      namemap[20333]	= "f_1(1420)";
      namemap[20413]	= "D_1(H)^+";
      namemap[-20413]	= "D_1(H)^-";
      namemap[20423]	= "D_1(2430)^0";
      namemap[-20423]	= "D_1(2430)~^0";
      namemap[20433]	= "D_s1(2460)^+";
      namemap[-20433]	= "D_s1(2460)^-";
      namemap[20443]	= "chi_c1(1P)";
      namemap[20513]	= "B_1(H)^0";
      namemap[-20513]	= "B_1(H)~^0";
      namemap[20523]	= "B_1(H)^+";
      namemap[-20523]	= "B_1(H)^-";
      namemap[20533]	= "B_s1(H)^0";
      namemap[-20533]	= "B_s1(H)~^0";
      namemap[20543]	= "B_c1(H)^+";
      namemap[-20543]	= "B_c1(H)^-";
      namemap[20553]	= "chi_b1(1P)";
      namemap[20555]	= "Upsilon_2(1D)";
      namemap[21112]	= "Delta(1910)^-";
      namemap[21114]	= "Delta(1920)^-";
      namemap[21212]	= "Delta(1910)^0";
      namemap[21214]	= "N(1700)^0";
      namemap[22112]	= "N(1535)^0";
      namemap[22114]	= "Delta(1920)^0";
      namemap[22122]	= "Delta(1910)^+";
      namemap[22124]	= "N(1700)^+";
      namemap[22212]	= "N(1535)^+";
      namemap[22214]	= "Delta(1920)^+";
      namemap[22222]	= "Delta(1910)^++";
      namemap[22224]	= "Delta(1920)^++";
      namemap[23112]	= "Sigma(1750)^-";
      namemap[-23112]	= "Sigma~(1750)^-";
      namemap[23114]	= "Sigma(1940)^-";
      namemap[-23114]	= "Sigma~(1940)^-";
      namemap[23122]	= "Lambda(1600)^0";
      namemap[-23122]	= "Lambda~(1600)^0";
      namemap[23124]	= "Lambda(1890)^0";
      namemap[-23124]	= "Lambda~(1890)^0";
      namemap[23126]	= "Lambda(2110)^0";
      namemap[-23126]	= "Lambda~(2110)^0";
      namemap[23212]	= "Sigma(1750)^0";
      namemap[-23212]	= "Sigma~(1750)^0";
      namemap[23214]	= "Sigma(1940)^0";
      namemap[-23214]	= "Sigma~(1940)^0";
      namemap[23222]	= "Sigma(1750)^+";
      namemap[-23222]	= "Sigma~(1750)^+";
      namemap[23224]	= "Sigma(1940)^+";
      namemap[-23224]	= "Sigma~(1940)^+";
      namemap[30113]	= "rho(1700)^0";
      namemap[30213]	= "rho(1700)^+";
      namemap[-30213]	= "rho(1700)^-";
      namemap[30223]	= "omega(1650)";
      namemap[30313]	= "K*(1680)^0";
      namemap[-30313]	= "K*(1680)~^0";
      namemap[30323]	= "K*(1680)^+";
      namemap[-30323]	= "K*(1680)^-";
      namemap[30443]	= "psi(3770)";
      namemap[30553]	= "Upsilon_1(1D)";
      namemap[31114]	= "Delta(1600)^-";
      namemap[31214]	= "N(1720)^0";
      namemap[32112]	= "N(1650)^0";
      namemap[32114]	= "Delta(1600)^0";
      namemap[32124]	= "N(1720)^+";
      namemap[32212]	= "N(1650)^+";
      namemap[32214]	= "Delta(1600)^+";
      namemap[32224]	= "Delta(1600)^++";
      namemap[33122]	= "Lambda(1670)^0";
      namemap[-33122]	= "Lambda~(1670)^0";
      namemap[42112]	= "N(1710)^0";
      namemap[42212]	= "N(1710)^+";
      namemap[43122]	= "Lambda(1800)^0";
      namemap[-43122]	= "Lambda~(1800)^0";
      namemap[53122]	= "Lambda(1810)^0";
      namemap[-53122]	= "Lambda~(1810)^0";
      namemap[100111]	= "pi(1300)^0";
      namemap[100113]	= "rho(1450)^0";
      namemap[100211]	= "pi(1300)^+";
      namemap[-100211]	= "pi(1300)^-";
      namemap[100213]	= "rho(1450)^+";
      namemap[-100213]	= "rho(1450)^-";
      namemap[100221]	= "eta(1295)";
      namemap[100223]	= "omega(1420)";
      namemap[100311]	= "K(1460)^0";
      namemap[-100311]	= "K(1460)~^0";
      namemap[100313]	= "K*(1410)^0";
      namemap[-100313]	= "K*(1410)~^0";
      namemap[100321]	= "K(1460)^+";
      namemap[-100321]	= "K(1460)^-";
      namemap[100323]	= "K*(1410)^+";
      namemap[-100323]	= "K*(1410)^-";
      namemap[100325]	= "K_2(1980)^+";
      namemap[-100325]	= "K_2(1980)^-";
      namemap[100331]	= "eta(1475)";
      namemap[100333]	= "phi(1680)";
      namemap[100411]	= "D(2S)^+";
      namemap[-100411]	= "D(2S)^-";
      namemap[100413]	= "D*(2S)^+";
      namemap[-100413]	= "D*(2S)^+";
      namemap[100421]	= "D(2S)^0";
      namemap[-100421]	= "D(2S)~^0";
      namemap[100423]	= "D*(2S)^0";
      namemap[-100423]	= "D*(2S)~^0";
      namemap[100441]	= "eta_c(2S)";
      namemap[100443]	= "psi(2S)";
      namemap[100445]	= "chi_c2(2P)";
      namemap[100551]	= "eta_b(2S)";
      namemap[100553]	= "Upsilon(2S)";
      namemap[100555]	= "chi_b2(2P)";
      namemap[100557]	= "Upsilon_3(2D)";
      namemap[110551]	= "chi_b0(2P)";
      namemap[110553]	= "h_b(2P)";
      namemap[110555]	= "eta_b2(2D)";
      namemap[120553]	= "chi_b1(2P)";
      namemap[120555]	= "Upsilon_2(2D)";
      namemap[130553]	= "Upsilon_1(2D)";
      namemap[200551]	= "eta_b(3S)";
      namemap[200553]	= "Upsilon(3S)";
      namemap[200555]	= "chi_b2(3P)";
      namemap[210551]	= "chi_b0(3P)";
      namemap[210553]	= "h_b(3P)";
      namemap[220553]	= "chi_b1(3P)";
      namemap[300553]	= "Upsilon(4S)";
  }

  std::string name(int pdgid)
  {
    if ( namemap.find(pdgid) != namemap.end() )
      return namemap[pdgid];
    else
      return string("not defined");
  }
};
'''
gROOT.ProcessLine(PDGCode)
pdg = PDG()

