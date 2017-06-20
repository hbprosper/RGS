//----------------------------------------------------------------------------
//  File:    RGS.cc
//  Purpose: Implement the Random Grid Search algorithm. This code
//           can be called from Python and Root.
//           
//  Created: 18-Aug-2000 Harrison B. Prosper, Chandigarh, India. A rewrite of
//                       code written in early 1995.
//
//  Updated: 05-Apr-2002 HBP tidy up
//           17-May-2006 HBP use weightindex instead of a vector of weights
//           11-Aug-2012 HBP & Sezen - generalize
//           18-Jan-2015 HBP - add optional event selection
//           04-Apr-2015 HBP & Sezen - improve screen printout
//                       improve comments
//           06-Apr-2015 HBP - add functions to save either to a text file or
//                       to an ntuple.
//           28-May-2016 HBP - minor update
//           02-Jun-2016 HBP - add option overall weighting of events
//           19-Feb-2017 HBP - add histogram with counts (from 30,000 feet!)
//           17-Jun-2017 HBP - rename ladder to staircase, box to twosided
//----------------------------------------------------------------------------
#include <stdio.h>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <set>
#include <map>

#include "TRandom3.h"
#include "TLeaf.h"
#include "TTreeFormula.h"
#include "TH1F.h"
#include "RGS.h"

using namespace std;

// define DBRGS at the command line in order to tun on debugging
static int DEBUG = getenv("DBRGS") > (char*)0 ? atoi(getenv("DBRGS")) : 0;

// return RGS version number
string rgsversion() 
{
  return string("RGS v2.5");
}
//----------------------------------------------------------------------------
// A bunch of simple string-based utilities
//----------------------------------------------------------------------------
void error(string message)
{
  cerr << "** ERROR ** " << message << endl;
  exit(0);
}

// True if string "strg" contains string "str"
bool inString(string strg, string str)
{
  int j = strg.find(str,0);
  return (j>-1) && (j<(int)strg.size());
}

// Extract name of a file with neither its extension
// nor its path
string nameonly(string filename)
{
  int i = filename.rfind("/");
  int j = filename.rfind(".");
  if ( j < 0 ) j = filename.size();
  return filename.substr(i+1,j-i-1);
}

// Strip spaces and tabs from either end of a string
string strip(string line)
{
  int l = line.size();
  if ( l == 0 ) return string("");
  int n = 0;
  while (((line[n] == 0)    ||
	  (line[n] == ' ' ) ||
	  (line[n] == '\n') ||
	  (line[n] == '\t')) && n < l) n++;
  
  int m = l-1;
  while (((line[m] == 0)    ||
	  (line[m] == ' ')  ||
	  (line[m] == '\n') ||
	  (line[m] == '\t')) && m > 0) m--;
  return line.substr(n,m-n+1);
}

// Split string "str" into left and right parts using the delimeter "delim"
void bisplit(string str, string& left, string& right, string delim,
	     int direction=-1)
{
  left  = str;
  right = "";
  int i = 0;
  if ( direction > 0 )
    i = str.rfind(delim);
  else      
    i = str.find(delim);
  if ( i > 0 )
    {
      left  = str.substr(0, i);
      right = str.substr(i+delim.size(), str.size()-i-delim.size());
    }
}

//----------------------------------------------------------------------------
// Description: Read from a text file or a simple ROOT ntuple.
// If treename is not given
// Created: 03-May-2005 Harrison B. Prosper
//----------------------------------------------------------------------------
bool slurpTable(string filename,
		vector<string>& header, 
		vector<vector<double> >& data,
		int start,
		int count,
		string treename,
		string selection,
		string  weightname,
		double  fileweight,		
		double& tot,
		double& err,
		int&    weightindex)		
{
  cout << "\t" << filename << endl;
  header.clear();
  weightindex = -1;
  
  // If treename given assume this is an ntuple
  bool ntuple = treename != "";

  // Get maximum numer of rows or events available
  int nrow=0;     // number of rows read
  int maxrows=0;
  if ( ntuple )
    {
      // We assume this is a ROOT ntuple since a tree name has been
      // given
      TFile rfile(filename.c_str());
      if ( !rfile.IsOpen() )
	error("slurpTable - unable to open "+filename);
      
      TTree* tree = (TTree*)rfile.Get(treename.c_str());
      if ( !tree )
	error("slurpTable - unable to get tree "+treename);
      
      // Read "count" rows if count > 0, otherwise read all lines
      maxrows = tree->GetEntries();
    }
  else
    {
      // We assume this is a simple text file, with a header followed by
      // rows of data
      ifstream stream(filename.c_str());
      if ( ! stream.good() )
	{ 
	  error("slurpTable - unable to open "+filename);
	  return false;
	}
      // count rows
      string line;
      while ( getline(stream, line, '\n') ) maxrows++;
    }

  // Determine number of rows or events to read
  if ( count < 0 )
    {
      count = maxrows;
    }
  else if ( count == 0 )
    {
      error("What, ZERO entries requested? Thou lump of foul deformity!");
    } 
  else
    {
      count = count > 0 ? min(count, maxrows) : maxrows;
    }
  
  // I'm alive printout - print row count every "step" entries
  int step = 10000;
  tot = 0.0;
  err = 0.0;
  
  if ( ntuple )
    {
      // We assume this is a simple ROOT ntuple
      TFile rfile(filename.c_str());
      if ( !rfile.IsOpen() )
	error("slurpTable - unable to open "+filename);

      TTree* tree = (TTree*)rfile.Get(treename.c_str());
      if ( !tree )
	error("slurpTable - unable to get tree "+treename);
      
      tree->ResetBranchAddresses();
            
      // Optionally, allow for event selection when using trees
      TTreeFormula* keep = 0;
      bool selectEvent = selection != "";
      if ( selectEvent ) 
	keep = new TTreeFormula("selection", 
				selection.c_str(), tree);
      
      // Get branches
      TObjArray* branches = tree->GetListOfBranches();
      int nbranches = branches->GetEntries();
      vector<int>    ibuffer(nbranches);
      vector<float>  fbuffer(nbranches);
      vector<double> dbuffer(nbranches);
      vector<char>   vtype(nbranches);

      for(int i=0; i < nbranches; i++)
	{
	  // Assume simple ntuple with leaf name = branch name
	  TBranch* branch = (TBranch*)(branches->At(i));
	  header.push_back(branch->GetName());
	  
	  if ( header.back() == weightname ) weightindex = i;
	  
	  TLeaf* leaf = branch->GetLeaf(branch->GetName());

	  vtype[i] = leaf->GetTypeName()[0];
	  if      ( vtype[i] == 'I' )
	    tree->SetBranchAddress(header.back().c_str(), &ibuffer[i]);
	  else if ( vtype[i] == 'F' )
	    tree->SetBranchAddress(header.back().c_str(), &fbuffer[i]);
	  else if ( vtype[i] == 'D' )
	    tree->SetBranchAddress(header.back().c_str(), &dbuffer[i]);
	}

      // Loop "count" entries, starting at start.
      for(int row=start; row < maxrows; row++)
	{
	  tree->GetEntry(row);
	  // Apply selection if given
	  if ( selectEvent )
	    {
	      // CHECK!
	      if ( ! (keep->EvalInstance(row) > 0) ) continue;
	    }

	  // Convert floats and ints into doubles
	  for(unsigned int i=0; i < dbuffer.size(); i++)
	    {	     
	      if ( vtype[i] == 'I' )
		dbuffer[i] = ibuffer[i];
	      else if ( vtype[i] == 'F' )
		dbuffer[i] = fbuffer[i];
	    }
	  // Cache row
	  data.push_back(dbuffer);
	  
	  double w = fileweight;
	  if ( weightindex > -1 ) w *= dbuffer[weightindex];

	  tot += w;
	  err += w * w;

	  // Increment number of rows read
	  nrow++;
	  
	  if ( nrow % step == 0 )
	    cout << "\t\trows selected: " << nrow
		 << "\tread: " << row+1
		 << endl;

          if ( count < 0 ) continue;  // read all remaining rows
          if ( nrow >= count ) break; // read "count" rows

	}
      if ( keep ) delete keep;
    }
  else
    {
      // We assume this is a simple text file, with a header followed by
      // rows of data
      ifstream stream(filename.c_str());
      if ( ! stream.good() )
	{ 
	  error("slurpTable - unable to open "+filename);
	  return false;
	}

      // Read header
      string line;
      getline(stream, line, '\n');
      istringstream inp(line);
      int c = 0;
      while ( inp >> line )
	{
	  header.push_back(line);
	  if ( header.back() == weightname ) weightindex = c;
	  c++;
	}

      // Skip the first "start" lines
      int n=0;
      for(int i=0; i < start; i++)
	{
	  n++;
	  if ( !getline(stream, line, '\n') ) break;
	}
   
      vector<double> d(header.size());
      while ( getline(stream, line, '\n') )
        {	  
          istringstream inp2(line);      
          for(int i=0; i < (int)header.size(); i++) inp2 >> d[i];
          data.push_back(d);

	  double w = fileweight;
	  if ( weightindex > -1 ) w *= d[weightindex];

	  tot += w;
	  err += w * w;

	  if ( nrow % step == 0 ) cout << "\t\trows read: " << nrow << endl;
	  // Increment number of rows read
	  nrow++;	  
          if ( count < 0 ) continue; // read all remaining rows
          if ( nrow >= count ) break;
        }
      stream.close();
    }
  cout << "\t\t     selected: " << nrow << endl << endl;
  return true;
}

//-----------------------------
// CONSTRUCTORS
//-----------------------------

RGS::RGS()
  : _status(0)
{}
/**
   cutdatafilenames - The names of one or more files containing cuts
   start     - The starting row for reading from the file(s).
   numrows   - The number of rows to be cached. The number of
               rows read could be greater if a selection has been
               applied when using a ROOT ntuple.
   treename  - The name of the ROOT tree. If omitted, we assume the
               file to be read is a text file.
   selection - Optional selection string to be used to select rows. This
               works only with ROOT ntuples.
 */
RGS::RGS(vstring& cutdatafilenames, int start, int numrows, 
	 string treename,
	 string weightname,
	 string selection)
  : _status(0),
    _treename(treename),
    _weightname(weightname),
    _selection(selection),
    _weightindex(vector<int>()),
    _weight(vector<double>()),
    _totals(vector<double>()),
    _errors(vector<double>())
{
  // Definitions:
  //  Cut-point
  //    The AND of a set of cuts
  //
  //  Cut data file
  //    n-tuple or text file containing the variables that define the cuts.
  //
  //  Search file
  //    The file on which the cuts are to be applied
  //
  _init(cutdatafilenames, start, numrows, _treename, _selection);
}

// See description for other constructor
RGS::RGS(string cutdatafilename, int start, int numrows, 
	 string treename,
	 string weightname,
	 string selection)
  : _status(0),
    _treename(treename),
    _weightname(weightname),
    _selection(selection),
    _weightindex(vector<int>()),
    _weight(vector<double>()),    
    _totals(vector<double>()),
    _errors(vector<double>())    
{
  vstring cutdatafilenames(1, cutdatafilename);
  _init(cutdatafilenames, start, numrows, _treename, _selection);
}

//-----------------------------
// DESTRUCTOR
//-----------------------------
RGS::~RGS(){}


//-----------------------------
// METHODS
//-----------------------------

bool
RGS::good() { return _status == 0; }


// Read file and cache data to which cuts are to be applied
void 
RGS::add(string searchfilename, 
         int    start, 
         int    numrows,
	 string resultname,
	 double weight,
	 string selection)
{
  // Cache name of file (stripped of path and extension)
  _searchname.push_back(nameonly(searchfilename));

  // Cache name of prefix for count and fraction results
  // If resultname is blank, we shall assume "count" and "fraction"
  // to which we append the ordinal value of the search file, starting
  // at zero.
  _resultname.push_back(resultname);

  // Cache weight/file
  _weight.push_back(weight);
  
  // Create an empty buffer of vector<vector<double> >
  // that serves as a cache of the data to which cuts
  // are to be applied
  _searchdata.push_back(vvdouble());

  // Cache column number of weights
  _weightindex.push_back(-1); // index to weight field

  _totals.push_back(0);
  _errors.push_back(0);
  
  _status = 0;

  // Create an alias of _searchdata.back() (called "sdata"),  NOT a copy!
  vector<vector<double> >& sdata = _searchdata.back();


  // If no selection specified, we shall assume that the cached
  // selection should be used.
  if ( selection == "" ) selection = _selection;
  
  // Ok, suck up data!
  cout << "\nRGS: Reading search data from file:" << endl;
  vector<string> header;
  if ( ! slurpTable(searchfilename, header, sdata, start, numrows, 
		    _treename, 
		    selection,
		    _weightname,
		    _weight.back(),
		    _totals.back(),
		    _errors.back(),
		    _weightindex.back()) )
    error("RGS: unable to read file " + searchfilename);

  cout << "\tSearch data will be identified with " << resultname << " in the RGS results file." << endl;

}

// Read files and cache data to which cuts are to be applied
void 
RGS::add(vector<string>& searchfilenames, 
         int start, 
	 int numrows,
	 string resultname,
	 double weight,
	 string selection)
{
  // See description in add(...) method above.
  _searchname.push_back(nameonly(searchfilenames[0]));
  _resultname.push_back(resultname);  
  _searchdata.push_back(vvdouble());
  _weightindex.push_back(-1); // Index to weight field
  _totals.push_back(0);
  _errors.push_back(0);
  
  _status = 0;
  vector<vector<double> >& sdata = _searchdata.back();

  // If no selection specified, we shall assume that the cached
  // selection should be used.
  if ( selection == "" ) selection = _selection;
  
  // ok, suck up data to be searched
  cout << "\nRGS: Reading data from file(s):" << endl;
  for(unsigned int ifile=0; ifile < searchfilenames.size(); ifile++)
    {
      // Cache weight/file
      _weight.push_back(weight);
      vector<string> header;
      if ( ! slurpTable(searchfilenames[ifile], header,
			sdata, start, numrows, 
			_treename, 
			selection,
			_weightname,
			_weight.back(),
			_totals.back(),
			_errors.back(),
			_weightindex.back()) )
          error("RGS: unable to read file " + searchfilenames[ifile]);
    }
}

// Read and decode variables file, then call the RGS workhorse

void RGS::run(string varfilename, // variables file name
              int nprint)     // print every nprint, just for fun!
{

  cout << "\nRGS: Reading cut variables from file: " << endl;
  cout << "\t" << varfilename << endl;

  // A cut-point is and AND of cuts
  //
  // Cut directions/types
  //    >    0 
  //    <    1
  //    >|   2
  //    <|   3
  //    ==   4
  //    <>   5
  //
  // Syntax of variables file:
  //
  // 1. one-sided cut
  //      variable-name cut-direction
  //
  // 2. two-sided cut 
  //      variable-name <>
  //
  // 3. staircase cut
  //
  //      \staircase number of cut-points
  //          variable-name-1 cut-direction-1
  //          variable-name-2 cut-direction-2
  //                   : :
  //      \end

  // cutvar will contain the first field (the variable name)
  // cutdir will contain the remaining fields (either the cut-direction
  // or if cutvar is the keyword "\staircase", it will be the number of
  // cut-points that comprise the staircase, that is, the number of
  // cut-points to be ORed.
  //
  // Note: a staircase cut can comprise any combination of uni-directional cuts

  vstring cutvar;
  vstring cutdir;

  // read and decode the variables file, the file that defines the
  // variables comprising the cuts and the direction(s) of the cut.
  ifstream fin(varfilename.c_str());
  if ( ! fin.good() )
    error("RGS: unable to open file: " + varfilename);

  string line;
  while ( getline(fin, line) )
    {
      // skip blank lines and commented lines
      line = strip(line);
      if ( line == "" ) continue;
      if ( line.substr(0,1) == "/" ) continue;
      if ( line.substr(0,1) == "#" ) continue;
      if ( line.substr(0,1) == "%" ) continue;

      // split line into first field and the rest
      istringstream sin(line);
      string left, right;
      sin >> left;
      if ( left.substr(0,2) == "\e" ) 
        right="";
      else
        sin >> right;

      cutvar.push_back(left);
      cutdir.push_back(right);

      if ( DEBUG > 0 )
        cout << "varname(" << cutvar.back() << ") "
             << "cutdir(" << cutdir.back() << ") "
             << endl;
    }

  // ok, now do the real work
  run(cutvar, cutdir, nprint);
}

//-----------------------------
// The workhorse
//-----------------------------
void RGS::run(vstring&  cutvar,        // Variables defining cuts 
              vstring&  cutdir,        // Cut direction
              int nprint)
{
  if ( DEBUG > 0 )
    cout << endl << "RUN algorithm" << endl;

  // Make sure length of cutvar and cutdir are the same
  // Note: for staircase cuts, cutvar should contain the keyword 
  // (\staircase and \end) that delimit the staircase cuts.

  if ( cutvar.size() != cutdir.size() )
    error("RGS: length(cutvar) != length(cutdir)");

  // Clear indices from the cut variable name to its column number 
  // in the data to be searched.
  _index.clear();

  // The number of cut-points / cut (multiple cuts-points for a
  // a staircase cut, and 2 for a twosided cut)
  _cutpointcount.clear();

  // The indices of cut-points.
  // See how this is used below.
  _cutpointindex.clear();

  // codes defining nature of cut
  _cutcode.clear();

  _status = 0;

  // Initialize buffers for number of (possibly weighted)
  // events passing cuts
  // _totals.resize(_searchdata.size(),0.0);
  _errors.resize(_searchdata.size(),0.0);
  _counts.resize(_searchdata.size(),vdouble());
  for (int i = 0; i < (int)_counts.size(); i++)
    _counts[i].resize(_cutdata.size(), 0.0);

  // ----------------------------------------------------------
  // Decode cuts
  // ----------------------------------------------------------
  if ( DEBUG > 0 )
    cout << endl << "DECODE CUTS" << endl;

  // Default number of simultaneous cuts-points to consider.
  // However, by the end of the following loop, maxpoints will
  // have increased if we have at least one twosided or staircase cut.
  int maxpoints = 1; 

  bool staircaseActive = false; // Set to true while decoding a staircase cut

  for (unsigned int i = 0; i < cutdir.size(); i++)
    {
      int code  =-1;
      int pointcount = 1; // number of cut-points/staircase or twosided cut
      
      bool startOfStaircase =
	cutvar[i].substr(0, 2) == "\\s" ||
	cutvar[i].substr(0, 2) == "\\S" ||	
	cutvar[i].substr(0, 2) == "\\l" ||
	cutvar[i].substr(0, 2) == "\\L";

      bool endOfStaircase =
	cutvar[i].substr(0, 2) == "\\e" ||
	cutvar[i].substr(0, 2) == "\\E";

      // Check for the start of a staircase cut      
      if ( startOfStaircase )
        {
          // Found start of a staircase block.
          // NB: Nested staircases are not allowed
          if ( staircaseActive )
	    error("RGS: nested staircases are not allowed!");

          staircaseActive = true;
          code = STAIRCASE;
          istringstream sin(cutdir[i]);

          // Get number of cut-points/staircase
          sin >> pointcount;          
          if ( pointcount > maxpoints ) maxpoints = pointcount;
        }
      else if ( endOfStaircase )
        {
          // Found end of a staircase
          staircaseActive = false; // deactivate staircase state
          code = END;
        }
      else if ( inString(cutdir[i],"<>") )
        {
          // TWOSIDED cuts are not allowed in staircases
	  // (in current RGS version)
          if ( staircaseActive )
	    error("RGS: Two-sided cuts within staircasee not yet implemented!");

          code = TWOSIDED;
          pointcount = 2; // two points needed for a twosided cut
          if ( pointcount > maxpoints ) maxpoints = pointcount;
        }
      else if ( inString(cutdir[i],">") )
        {
          if   ( inString(cutdir[i],"|") )
            code = ABSGT;
          else
            code = GT;
        }
      else if ( inString(cutdir[i],"<") )
        {
          if   ( inString(cutdir[i],"|") )
            code = ABSLT;
          else
            code = LT;
        }
      else if ( inString(cutdir[i],"=") )
        {
          code = EQ;
        }

      // Save cut code and number of simultaneous cut-points
      // for current cut
      _cutcode.push_back(code);
      _cutpointcount.push_back(pointcount);

      if ( DEBUG > 0 )
	{
	  cout << "\tcutvar[" << i << "]: " << cutvar[i] << endl;
	  cout << "\t\tcode:       " << _cutcode.back() << endl;
	  cout << "\t\tpointcount: " << _cutpointcount.back() << endl;
	}
      
      // Get column index of current cut variable
      if       ( startOfStaircase )
        {
          _index.push_back(0); // Not used; so just set to any valid value
          //cout << "begin staircase cut" << endl;
        }
      else if ( endOfStaircase )
        {
          _index.push_back(0); // Not used; so just set to any valid value
          //cout << "end" << endl;
        }
      else
        {
          // Remember: index is the column number
          if ( _varmap.find(cutvar[i]) != _varmap.end() )
            _index.push_back(_varmap[cutvar[i]]); // index map into data
          else
	    error("RGS: cut variable "+cutvar[i]+" NOT found");
      
          // if ( code == TWOSIDED ) cout << "twosided cut" << endl;
          // cout << i << "\t" << _index.back() << "\t" 
          //      << cutvar[i] << "\t" << cutdir[i] 
          //      << " cut code " << code << endl;
        }
    } 
  if ( DEBUG > 0 )
    cout << "END DECODE CUTS" << endl << endl;
  
  // -------------------------------------------------------------------------
  // Here we augment each cut-point so that we can handle both twosided and
  // staircase cuts. For each cut-point, randomly select maxpoints-1 more 
  // cut-points.
  // -------------------------------------------------------------------------

  if ( DEBUG > 1 )
    cout << "SAMPLE CUT-POINTS" << endl << endl;

  TRandom3 ran;

  for (int cutpoint = 0; cutpoint < (int)_cutdata.size(); cutpoint++)
    {
      // Reserve space for indices of cut-points
      // maxpoints will be bigger than 1 if we have at least one twosided or
      // staircase cut.
      _cutpointindex.push_back(vector<int>(maxpoints));

      // The first position is the index of the original cut-point
      _cutpointindex[cutpoint][0] = cutpoint;

      // Randomly select maxpoints-1 additional cut-points.
      // We use a set to ensure no "collisions" between cut points,
      // that is, no duplicates.
      set<int> bucket;
      bucket.insert(cutpoint);

      if ( DEBUG > 1 )
        cout << cutpoint << "\t";

      int j = 1; // Initially there is only one cut-point/cut
      for(int i=0; i < 5*maxpoints; i++)
        {
          int jj = ran.Integer(_cutdata.size()-1);
          if ( bucket.find(jj) == bucket.end() )
            {
              // This randomly selected cut-point differs from
              // ones already selected, so add it to the selected set
              _cutpointindex[cutpoint][j] = jj;
              bucket.insert(jj);

              if ( DEBUG > 1 )
                cout << jj << "\t";

              j++;
              if ( j >= maxpoints ) break;
            }
        }
      if ( DEBUG > 1 )
        cout << endl;
    }
 
  if ( DEBUG > 1 )
    cout << "END SAMPLE CUT-POINTS" << endl << endl;

  // -------------------------------------------------------------------------
  // Loop over files to which cuts are to be applied
  // -------------------------------------------------------------------------
  cout << "\nRGS: Running the RGS algorithm..." << endl;
  for (int file = 0; file < (int)_searchdata.size(); file++)
    {
      // column index of weights (-1 if none given)
      int  weightindex = _weightindex[file];

      // NB: use a reference to avoid expensive copy of _searchdata[file]
      vvdouble& sdata  = _searchdata[file];

      cout << "\n\tRGS: Applying " << _cutdata.size() << " cuts to dataset "
	   << _searchname[file] 
           << " containing " << sdata.size() << " events" << endl;

      if ( nprint > 0 )
        {
          cout << "\tcuts to apply: " << endl;
          for (int cut = 0; cut < (int)cutvar.size(); cut++)
            {
              if ( _cutcode[cut] == STAIRCASE )
                {
                  cout << "\t\tstaircase" << endl;
                  continue;
                }
              else if ( _cutcode[cut] == END )
                {
                  cout << "\t\tend" << endl;
                  continue;
                }
              cout << "\t\t  " << _var[_index[cut]] 
                   << "\t" << cutdir[cut] << endl;
            }
          cout << endl;
        }

      // Check whether to use event weighting
      bool useEventWeight = weightindex > -1;

      // ------------------------------------
      // Loop over cut-points (==cutpoint)
      // ------------------------------------      
      int maxcuts = (int)_cutcode.size();

      for (int cutpoint = 0; cutpoint < (int)_cutdata.size(); ++cutpoint)
        {
#ifdef RGSDEBUG
          if ( DEBUG > 2 )
            cout << "\tCUT-POINT " << cutpoint 
                 << "\t========================== " << endl << endl;
#endif
	  // I'm alive printout every nprint cuts
          if ( (nprint > 0) &&
	       ((cutpoint % nprint == 0) ||
		(size_t)cutpoint==(_cutdata.size()-1)))
            cout << "\t\tapplying cut " << cutpoint << endl;
          
          // Loop over rows of file to be processed
          /////////////////////////////////////////
          
          _counts[file][cutpoint] = 0; // Count per file per cut-point
          
          for (int row = 0; row < (int)sdata.size(); row++)
            {
#ifdef RGSDEBUG
              if ( DEBUG > 2 )
                {
                  cout << "\tROW " << row << endl;
                }
#endif
              // Loop over cut values. 

              // If we are processing a staircase cut then on exit
              // from _staircasecut, cut should point to the end of
              // staircase 

              bool passed = true;

              int cut = 0;
              while ( cut < maxcuts )
                {
                  // -------------------------------------
                  // 1. get cut value for simple cut
                  // -------------------------------------
                  // column index into data vector
                  int jcut =_index[cut];

                  // datum value
                  float x  = sdata[row][jcut];

                  // cut value
                  float xcut =_cutdata[cutpoint][jcut];

#ifdef RGSDEBUG
                  if ( DEBUG > 2 )
                    {
                      switch (_cutcode[cut])
                        {
                        case GT:
                        case LT:
                        case ABSGT:
                        case ABSLT:
                        case EQ:
                          cout << "\t\t\t" 
                               << cutvar[cut] << "\t" << x << " "
                               << cutdir[cut] << " " << xcut << endl;
                          break;
                        default:
                          break;
                        }
                    }
#endif
                  // -------------------------------------
                  // 2. apply cut. If this is a twosided
                  // or a staircase cut, we need to loop over
                  // multiple cut-points. The looping is 
                  // handled by the methods _twosidedcut and 
                  // _staircasecut, respectively.
                  // -------------------------------------
                  switch (_cutcode[cut])
                    {
                    case GT:
                      passed = x > xcut;
                      break;
                      
                    case LT:
                      passed = x < xcut;
                      break;
                      
                    case ABSGT:
                      passed = abs(x) > abs(xcut);
                      break;
                      
                    case ABSLT:
                      passed = abs(x) < abs(xcut);
                      break;

                    case EQ:
                      passed = x == xcut;
                      break;

                    case TWOSIDED:
                      // Note: need to use index jcut, not cut
                      passed = _twosidedcut(x, cutpoint, jcut);
                      break;

                    case STAIRCASE:
                      // Upon exit, cut should point to end of staircase
                      passed = _staircasecut(sdata[row], cutpoint, cut);
                      break;
                      
                    default:
                      break;
                    }
                  
                  // If any cut fails, there is no point continuing.
                  if ( !passed ) break;

                  // IMPORTANT: remember to increment cut number
                  cut++;
                }
          
              // Keep a running sum of events that pass all cuts
              // of current cut-point
	      
	      double weight = _weight[file];
              if ( useEventWeight ) weight = weight * sdata[row][weightindex];
              
              // if ( cutpoint == 0 )
	      // 	{
	      // 	  _totals[file] += weight;
	      // 	  _errors[file] += weight*weight;
	      // 	}
          
              if ( passed ) _counts[file][cutpoint] += weight;

#ifdef RGSDEBUG
              if ( DEBUG > 2 )
                {
                  if ( passed )
                    cout << "\t\t\t\tPASSED" << endl << endl;
                  else
                    cout << "\t\t\t\tFAILED" << endl << endl;
                }
#endif
            }

        }
    }
}


// Twosided cut:    xcut_low < c < xcut_high
inline
bool RGS::_twosidedcut(float x, int cutpoint, int jcut)
{
  // get cut-points
  int cutpoint_1 =_cutpointindex[cutpoint][0];
  int cutpoint_2 =_cutpointindex[cutpoint][1];

  // get cut-values
  float xcut1   =_cutdata[cutpoint_1][jcut];
  float xcut2   =_cutdata[cutpoint_2][jcut];

  // get lower bound
  float xcutlow = min(xcut1, xcut2);

  bool passed = x > xcutlow;
  if ( passed ) 
    {
      // ok we've passed the lower bound, so check upper bound
      float xcuthigh = max(xcut1, xcut2);
      passed = x < xcuthigh;
    }
#ifdef RGSDEBUG
  if ( DEBUG > 2 )
    {
      float xcuthigh = max(xcut1, xcut2);
      cout << "\t\t   TWOSIDED" << endl;
      cout << "\t\t\t" 
           << _var[jcut] << "\t" << x << " > " << xcutlow << endl;
      cout << "\t\t\t" 
           << _var[jcut] << "\t" << x << " < " << xcuthigh << endl;
    }
#endif
  return passed;
}

// Staircase cut:    OR of cut-points
// Loop over cut-points for current cut
// IMPORTANT: on exit, the cut number should point to the 
// end of the staircase indicated by the cutcode END.

inline
bool RGS::_staircasecut(vdouble& datarow, int origcutpoint, int& cut)
{
#ifdef RGSDEBUG
  if ( DEBUG > 2 )
    {
      cout << "\t\t   STAIRCASE" << endl;
    }
#endif

  // Loop over cut-points of staircase
  // staircasepassed will be true if at least one cut-point of the
  // staircase returns true

  bool staircasepassed = false;

  // Recall: _cutcode contains all cut codes including the codes
  // for keywords \staircase and \end

  int maxcuts = (int)_cutcode.size();

  // Get number of cut points to use for current staircase
  int pointcount = _cutpointcount[cut];

#ifdef RGSDEBUG
  if ( DEBUG > 2 )
    {
      cout << "\t\t     cut-points = " << pointcount << endl;
    }
#endif

  // Remember first cut of staircase
  int firstcut = cut;
  firstcut++;

  for(int ii=0; ii < pointcount; ++ii)
    {
#ifdef RGSDEBUG
      if ( DEBUG > 2 )
        {
          cout << "\t\t       point " << ii << endl;
        }
#endif
      // Get cut-point index
      int cutpoint = _cutpointindex[origcutpoint][ii];

      bool endOfstaircase = false;
      bool passed = true;

      // loop over cuts, starting each time at first cut of staircase
      
      cut = firstcut;
      while ( cut < maxcuts )
        {
          // column index in data vector
          int jcut =_index[cut];

          // datum value
          float x  = datarow[jcut];

          // cut value
          float xcut =_cutdata[cutpoint][jcut];

#ifdef RGSDEBUG
          if ( DEBUG > 2 )
            {
              switch (_cutcode[cut])
                {
                case GT:
                  cout << "\t\t\t" 
                       << _var[jcut] << "\t" << x << " > " << xcut << endl;
                  break;
                case LT:
                  cout << "\t\t\t" 
                       << _var[jcut] << "\t" << x << " < " << xcut << endl;
                  break;
                case ABSGT:
                  cout << "\t\t\t" 
                       << _var[jcut] << "\t" << x << " >| " << xcut << endl;
                  break;
                case ABSLT:
                  cout << "\t\t\t" 
                       << _var[jcut] << "\t" << x << " <| " << xcut << endl;
                  break;
                case EQ:
                  cout << "\t\t\t" 
                       << _var[jcut] << "\t" << x << " == " << xcut << endl;
                  break;
                default:
                  break;
                }
            }
#endif
          switch (_cutcode[cut])
            {
            case GT:
              passed = x > xcut;
              break;
              
            case LT:
              passed = x < xcut;
              break;
              
            case ABSGT:
              passed = abs(x) > abs(xcut);
              break;
              
            case ABSLT:
              passed = abs(x) < abs(xcut);
              break;
              
            case EQ:
              passed = x == xcut;
              break;
              
            case END:
              endOfstaircase = true;
              break;

            default:
              break;
            }

          // If any cut fails, there is no point continuing
          if ( ! passed ) break;

          // break out of loop over cuts if we have reached end of staircase
          if ( endOfstaircase ) break;

          // IMPORTANT: Remember to increment cut number
          cut++;
        }

      // Make sure we are at the end of the staircase
      if ( ! endOfstaircase )
        {
          cut++;
          while ( cut < maxcuts )
            {
              if ( _cutcode[cut] == END ) break;
              cut++;
            }
        }
#ifdef RGSDEBUG
      if ( DEBUG > 2 )
        {
          if ( passed )
            cout << "\t\t\t\t\tpassed " << endl;
          else
            cout << "\t\t\t\t\tfailed " << endl;
        }
#endif
      // Take OR of cuts over cut-points
      staircasepassed = staircasepassed || passed;

      // If a cut succeeds, there is no need to continue
      if ( staircasepassed ) break;
    }

#ifdef RGSDEBUG
  if ( DEBUG > 2 )
    {
      cout << "\t\t   END" << endl;
    }
#endif
  return staircasepassed;
}

void
RGS::save(string resultfilename)
{
  bool saveToNtuple = inString(resultfilename, ".root");
  if ( saveToNtuple )
    _saveToNtupleFile(resultfilename);
  else
    _saveToTextFile(resultfilename);
}

void
RGS::_saveToTextFile(string resultfilename)
{
  cout << "\nRGS: Saving RGS results to text file: " << endl;
  cout << "\t" << resultfilename << endl;
  cout << "\tVariable \tsize\tcut" << endl;

  ofstream fout(resultfilename.c_str());

  // Declare a buffer of size, maxpoints x maxcuts, for writing out cut values
  int maxpoints = (int)_cutpointindex[0].size();
  int maxcuts   = (int)_cutcode.size();

  vector<vector<float> > cutvalue(maxcuts, vector<float>(maxpoints,0));

  // Create a header variable for each cut.
  // For twosided and staircase cuts, use fixed
  // length arrays of the appropriate size
  vector<string> varname; // name of cut variable
  vector<int> varsize;    // number of simultaneous cut-points
  vector<int> varcut;     // index into cutvalue vector
  int cut = 0;
  while ( cut < maxcuts )
    {
      int jcut = _index[cut];
      string var = _var[jcut];
      
      switch (_cutcode[cut])
        {
        case GT:
        case LT:
        case ABSGT:
        case ABSLT:
        case EQ:
          {
	    varname.push_back(var);
	    varsize.push_back(1);
	    varcut.push_back(cut);
            cout << "\t" << varname.back()
		 << "\t" << varsize.back()
		 << "\t" << varcut.back()
		 << "\t<== onesided"
		 << endl;	    
          }
          break;

        case TWOSIDED:
          {
	    varname.push_back(var);
	    varsize.push_back(2);
	    varcut.push_back(cut);
            cout << "\t" << varname.back()
		 << "\t" << varsize.back()
		 << "\t" << varcut.back()
		 << "\t<== twosided"	      
		 << endl;	    	    
          }
	  break;
	  
        case STAIRCASE:
          {
            // Get number of cut points to use for current staircase
            int pointcount = _cutpointcount[cut];
	    
            cut++; // Go to first cut in staircase
            while ( cut < maxcuts )
              {
                if ( _cutcode[cut] == END ) break;

		// Get name of cut variable
                jcut =_index[cut];
                var  =_var[jcut];
		
		varname.push_back(var);
		varsize.push_back(pointcount);
		varcut.push_back(cut);
		cout << "\t" << varname.back()
		     << "\t" << varsize.back()
		     << "\t" << varcut.back()
		     << "\t<== staircase"	      
		     << endl;		

                // IMPORTANT: Remember to increment cut number
                cut++;
              }
          }
          break;
                      
        default:
          break;
        }

      // IMPORTANT: increment loop counter
      cut++;
    }

  // Add one count and fraction variable for each file scanned
  for(unsigned int i=0; i < _counts.size(); i++)
    {
      char postfix[80];
      string name("");
      if ( _resultname[i] != "" )
	sprintf(postfix, "%s", _resultname[i].c_str());
      else
	sprintf(postfix, "%d", i);
      
      name = string("count")+string(postfix);
      varname.push_back(name);
      varsize.push_back(1);
      cout << "\t" << name << "\t" << varsize.back() << endl;		      

      name = string("fraction")+string(postfix);      
      varname.push_back(name);
      varsize.push_back(1);
      cout << "\t" << name << "\t" << varsize.back() << endl;		      
    }

  // Add a variable containing indices of cut points.
  varname.push_back("cutpointindex");
  varsize.push_back(_cutpointindex[0].size());
  cout << "\t" << varname.back() << "\t" << varsize.back() << endl;		      
  
  // Write variable names
  for(unsigned int i=0; i < varname.size(); i++)
    fout << varname[i] << " " << varsize[i] << "\t";
  fout << endl;
  
  // Now write out results
  for (unsigned int cutpoint=0; cutpoint < _cutdata.size(); cutpoint++)
    {
      cut = 0;
      while ( cut < maxcuts )
        {
          int jcut =_index[cut];
          
          switch (_cutcode[cut])
            {
            case GT:
            case LT:
            case ABSGT:
            case ABSLT:
            case EQ:
              cutvalue[cut][0] = _cutdata[cutpoint][jcut];
              break;

            case TWOSIDED:
              {
                // Get the two cut-points
                int cutpoint1  =_cutpointindex[cutpoint][0];
                int cutpoint2  =_cutpointindex[cutpoint][1];
                // Get the cut values
                cutvalue[cut][0] = _cutdata[cutpoint1][jcut];
                cutvalue[cut][1] = _cutdata[cutpoint2][jcut];
              }
              break;

            case STAIRCASE:
              {
                // Get number of cut points to use for current staircase
                int pointcount = _cutpointcount[cut];
		cut++;
                // Remember cut number of first cut in staircase
                int firstcut = cut;
		// Loop of points to be ORed and, for each, loop
		// over cuts to be ANDed. But make sure we store
		// cuts of the same kind contiguously. Birds of a feather
		// flock together!
                for(int i=0; i < pointcount; i++)
                  {
                    // Get cut-point using the index into the
		    // cutdata vector.
                    int cutpoint_i  =_cutpointindex[cutpoint][i];
		    // Remember to reset to first cut in staircase
		    // to be ANDed
                    cut = firstcut;
                    while ( cut < maxcuts )
                      {
                        jcut =_index[cut];
          
                        // Get the cut value
                        float xcut =_cutdata[cutpoint_i][jcut];
                        // Store in output buffer
                        cutvalue[cut][i] = xcut;
                        cut++;
                      }
                  }
              }
              break;
              
            default:
              break;              
            }
          cut++;
        }

      // Write cuts. Use varcut to "zero-suppress"
      // cutvalue vector, that is, to skip over entries
      // that are not in fact cut-point indices, but
      // are just placeholders for the staircase start
      // and end keywords.
      for(unsigned int j=0; j < varcut.size(); j++)
	{
	  int icut = varcut[j];
	  for(int i=0; i < varsize[j]; i++)
	    {
	      fout << cutvalue[icut][i] << "\t"; 
	    }
	}

      // write counts and fractions
      for(unsigned int i=0; i < _counts.size(); i++)
	{
	  float passcount = _counts[i][cutpoint];
	  float fraction  = passcount/total(i);
	  fout << passcount << "\t" << fraction << "\t"; 
	}

      // write cut-point indices
      for(unsigned int i=0; i < _cutpointindex[cutpoint].size(); i++)
	{
	  fout << _cutpointindex[cutpoint][i] << "\t";
	}
      
      fout << endl;
    }
  fout.close();
}


void
RGS::_saveToNtupleFile(string resultfilename)
{
  cout << endl << "RGS: Saving RGS results to ROOT file: " << endl;
  cout << "\t" << resultfilename << endl;

  TFile* file = new TFile(resultfilename.c_str(), "recreate"); 
  TTree* tree = new TTree("RGS", "RGS");

  // Declare a buffer of size, maxcuts x maxpoints, for writing out cut
  // values. NOTE: for staircase cuts, cutvalue will have entries that are
  // just placeholders for the staircase keywords. But this does not matter since
  // we pass the address of the appropriate subsets of cutvalue to Root.
  int maxpoints = (int)_cutpointindex[0].size();
  int maxcuts   = (int)_cutcode.size();
  vector<vector<float> > cutvalue(maxcuts, vector<float>(maxpoints,0));
  vector<int> cutpointindex(maxpoints);
  
  // Define buffers for counts and fractions
  // Note: _counts.size() = number of scanned files
  vector<float> counts(_counts.size());
  vector<float> fractions(_counts.size());

  // Create branches for each cut. For twosided and staircase cuts, use fixed
  // length arrays of the appropriate size
  int cut = 0;
  char name[80];
  char fmt[80];
  char postfix[80];  
  while ( cut < maxcuts )
    {
      int jcut =_index[cut];
      string var = _var[jcut];

      switch (_cutcode[cut])
        {
        case GT:
        case LT:
        case ABSGT:
        case ABSLT:
        case EQ:
          {
            sprintf(fmt, "%s/F", var.c_str());
            tree->Branch(var.c_str(), &cutvalue[cut][0], fmt);
            cout << "\t" << fmt << endl;
          }
          break;

        case TWOSIDED:
          {
            sprintf(fmt, "%s[2]/F", var.c_str());
            tree->Branch(var.c_str(), &cutvalue[cut][0], fmt);
            cout << "\t" << fmt << endl;
          }
	  break;
	  
        case STAIRCASE:
          {
            // Get number of cut points to use for current staircase
            int pointcount = _cutpointcount[cut];
            cut++; // go to first cut of staircase, i.e., the cut to be
	    // ANDed with the rest.
            while ( cut < maxcuts )
              {
                if ( _cutcode[cut] == END ) break;
                jcut =_index[cut];
                var  =_var[jcut];
                sprintf(fmt, "%s[%d]/F", var.c_str(), pointcount);
		// pass address of cutvalue to Root
                tree->Branch(var.c_str(), &cutvalue[cut][0], fmt);
                cout << "\t" << fmt << endl;

                // IMPORTANT: Remember to increment cut number
                cut++;
              }
          }
          break;
                      
        default:
          break;
        }

      // IMPORTANT: increment loop counter
      cut++;
    }
  // Create a branch for the cut indices 
  sprintf(fmt, "cutpointindex[%d]/I", (int)cutpointindex.size());
  tree->Branch("cutpointindex", &cutpointindex[0], fmt);
  cout << "\t" << fmt << endl;
  // Add one count variable for each file scanned

  for(unsigned int i=0; i < counts.size(); i++)
    {
      if ( _resultname[i] != "" )
	sprintf(postfix, "%s", _resultname[i].c_str());
      else
	sprintf(postfix, "%d", i);
      
      sprintf(name, "%s%s", "count", postfix);      
      sprintf(fmt,  "%s/F", name);
      tree->Branch(name, &counts[i], fmt);
      cout << "\t" << fmt << endl;

      sprintf(name, "%s%s", "fraction", postfix);      
      sprintf(fmt, "%s/F", name);
      tree->Branch(name, &fractions[i], fmt);
      cout << "\t" << fmt << endl;
    }

  // Now fill ntuple

  for (unsigned int cutpoint=0; cutpoint < _cutdata.size(); cutpoint++)
    {
      cut = 0;
      while ( cut < maxcuts )
        {
          int jcut =_index[cut];
          
          switch (_cutcode[cut])
            {
            case GT:
            case LT:
            case ABSGT:
            case ABSLT:
            case EQ:
              cutvalue[cut][0] = _cutdata[cutpoint][jcut];
              break;

            case TWOSIDED:
              {
                // Get the two cut-points
                int cutpoint_1  =_cutpointindex[cutpoint][0];
                int cutpoint_2  =_cutpointindex[cutpoint][1];
                // Get the cut values
                cutvalue[cut][0] =_cutdata[cutpoint_1][jcut];
                cutvalue[cut][1] =_cutdata[cutpoint_2][jcut];
              }
              break;

            case STAIRCASE:
              {
                // Get number of cut points to use for current staircase
                int pointcount = _cutpointcount[cut];
		cut++; // go to first cut of staircase, i.e., the cut to be
		// ANDed with the rest.	      
                // IMPORTANT: Remember first cut within staircase
                int firstcut = cut;

                for(int i=0; i < pointcount; i++)
                  {
                    // Get cut-point index, i.e., the index into
		    // the cutdata vector.
                    int cutpoint_i  =_cutpointindex[cutpoint][i];
		    
		    // Reset to first cut of staircase
                    cut = firstcut;
                    while ( cut < maxcuts )
                      {
                        jcut =_index[cut];
          
                        // Get the cut value
                        float xcut =_cutdata[cutpoint_i][jcut];
                        // Store in output buffer
                        cutvalue[cut][i] = xcut;
                        cut++;
                      }
                  }
              }
              break;
              
            default:
              break;              
            }
          cut++;
        }

      // store counts and fractions in
      // counts and fractions buffers, whose
      // addresses have are already known to Root.
      // (See Branch definitions above)
      for(unsigned int i=0; i < counts.size(); i++)
	{
	  counts[i]    =_counts[i][cutpoint];
	  fractions[i] = counts[i]/total(i);
	}

      for(unsigned int i=0; i < cutpointindex.size(); i++)
	{
	  cutpointindex[i] = _cutpointindex[cutpoint][i];
	}      
      file->cd();
      tree->Fill();
    }
  file->cd();

  // Fill a histogram with totals
  TH1F* htotal = new TH1F("total", "total",
			  _totals.size(), 0, _totals.size());
  for(size_t c=0; c < _totals.size(); c++)
    {
      htotal->SetBinContent(c+1, total(c));
      htotal->SetBinError(c+1, etotal(c));
    }
  //tree->AutoSave("SaveSelf");
  file->Write("", TObject::kOverwrite);
}

double
RGS::total(int index) 
{
  if ( index < 0 || index > (int)_totals.size()-1 )
    {
      _status = rBADINDEX;
      return 0;
    }
  return _totals[index];
}

double
RGS::etotal(int index) 
{
  if ( index < 0 || index > (int)_errors.size()-1 )
    {
      _status = rBADINDEX;
      return 0;
    }
  return sqrt(_errors[index]);
}


double
RGS::count(int index, int cutindex) 
{
  if ( index < 0 || index > (int)_totals.size()-1 )
    {
      _status = rBADINDEX;
      return 0;
    }
  if ( cutindex < 0 || cutindex > (int)_cutdata.size()-1 )
    {
      _status = rBADINDEX;
      return 0;
    }
  return _counts[index][cutindex];
}

vdouble
RGS::cuts(int index)
{ 
  if ( index < 0 || index > (int)_cutdata.size()-1 )
    {
      _status = rBADINDEX;
      return vdNULL;
    }
  vdouble cut(_cutcode.size());
  for(unsigned i=0; i < _cutcode.size(); i++)
    cut[i] = _cutdata[index][_index[i]];
  return cut;
}

vstring
RGS::cutvars()   
{ 
  vstring cutvar(_cutcode.size());
  for(unsigned i=0; i < _cutcode.size(); i++)
    cutvar[i] = _var[_index[i]];
  return cutvar;
}

int
RGS::ncuts() { return _cutdata.size(); }

vstring&
RGS::vars()   { return _var; }

int
RGS::ndata(int index)
{ 
  if ( index < 0 || index > (int)_searchdata.size()-1 )
    {
      _status = rBADINDEX;
      return 0;
    }
  return _searchdata[index].size();
}

vdouble&
RGS::data(int index, int event)
{ 
  if ( index < 0 || index > (int)_searchdata.size()-1 )
    {
      _status = rBADINDEX;
      return vdNULL;
    }
  return _searchdata[index][event]; 
}

//--------------------------------
// PRIVATE METHODS
//--------------------------------

void
RGS::_init(vstring& cutdatafilenames, int start, int numrows, 
	   string treename,
	   string selection)
{
  _cutdata.clear();
  _cutcode.clear();
  _searchname.clear();
  _searchdata.clear();
  _weightindex.clear();
  _varmap.clear();
  _totals.clear();
  _errors.clear();
  _counts.clear();

  _status = rSUCCESS;

  cout << endl << rgsversion() << endl << endl;
  cout << "\nRGS: Reading cut data from file(s):" << endl;
  vector<string> var;
  for(unsigned int ifile=0; ifile < cutdatafilenames.size(); ifile++)
    {
      var.clear();
      if ( ! slurpTable(cutdatafilenames[ifile], 
			var, _cutdata, start, numrows, 
			treename, selection) )
        {
          cout << "RGS **Error** unable to read cut-file "
	       << cutdatafilenames[ifile] << endl;
          _status = rFAILURE;
          return;
        }
    }
  for(int i = 0; i < (int)var.size(); i++)
    {
      _varmap[var[i]] = i;     // Map variable name to ordinal value
      _var.push_back(var[i]);  // Map ordinal value to variable name
    }
}
