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
                                            color,
                                            return_plot=True))
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
