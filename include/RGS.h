#ifndef RGS_H
#define RGS_H
//----------------------------------------------------------------------------
//  File:    RGS.h
//  Purpose: Implement the Random Grid Search algorithm. This code
//           can be called from Python and Root.
//
// Purpose: Declaration of RGS classes
// Created: 18-Aug-2000 Harrison B. Prosper, Chandigarh, India
// Updated  18-Jan-2015 HBP - add selection argument - Bari, Italy
//          06-Apr-2015 HBP - add functions to save either to a text file or
//                      to an ntuple.
//          19-Feb-2017 HBP - remove lumi argument from save and add
//                            error associated with total
//----------------------------------------------------------------------------
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iomanip>
#include "TFile.h"
#include "TTree.h"

typedef std::vector<int>     vint;
typedef std::vector<vint>    vvint;
typedef std::vector<float>   vfloat;
typedef std::vector<vfloat>  vvfloat;
typedef std::vector<double>  vdouble;
typedef std::vector<vdouble> vvdouble;
typedef std::vector<std::string>  vstring;
typedef std::map< std::string, int >  varmap;
static vdouble  vdNULL;
static double dNULL=0;
static int iNEG1=-1;

// ERROR CODES

const int rSUCCESS = 0;
const int rFAILURE =-1;
const int rBADINDEX=-2;
const int rEOF     =-6;
const int rPYTHONERROR =-7;

/** These codes define the available types of cuts.
 */ 
enum CUTCODE {GT, LT, ABSGT, ABSLT, EQ, TWOSIDED, STAIRCASE, END};

/** Read either a text file or a ROOT ntuple file that contains variable names and values.
    The text file should contain a header of variable names followed by rows of values.
    The ROOT ntuple file should contain branch (and leaf) names that are the same 
    as the variable names.
 */
bool slurpTable(std::string filename,
		std::vector<std::string>& header, 
		std::vector<std::vector<double> >& data,
		int start=0,
		int count=-1,
		std::string treename="",
		std::string selection="",
		std::string weightname="",
		double  fileweight=1,
		double& tot=dNULL,
		double& err=dNULL,
		int&    weightindex=iNEG1);

/** Return version number of RGS.
 */
std::string rgsversion();

/** Random Grid Search class.
 */
class RGS
{
public:
  RGS();

  /**
   */
  RGS(std::string cutdatafilename, int start=0, int numrows=0, 
      std::string treename="",
      std::string weightname="",
      std::string selection="");

  /**
   */
  RGS(std::vector<std::string>& cutdatafilename, int start=0, int numrows=0,
      std::string treename="",
      std::string weightname="",
      std::string selection="");

  virtual ~RGS();

  /// If true, all is well.
  bool  good();

  /// Add a data file to be searched.
  void  add(std::string searchfilename,
            int start=0, 
            int numrows=-1,
	    std::string resultname="",
	    double weight=1.0,
	    std::string selection="");
  
  /// Add one or more data files to be searched.
  void  add(std::vector<std::string>& searchfilename,
            int start=0, 
            int numrows=-1,
	    std::string resultname="",
	    double weight=1.0,
	    std::string selection="");

  /// Run the RGS algorithm for specified cut variables and cut directions.
  void  run(vstring&  cutvar,  // Variables defining cuts 
            vstring&  cutdir,  // Cut direction (cut-type)
            int nprint=500);

  /// Run the RGS algorithm for specified cut variables and cut directions.
  void  run(std::string  varfilename,  // file name of Variables file
            int nprint=500);

  /// Return the total count for the data file identified by dataindex.
  double    total(int dataindex);

  /// Return error in total count for the data file identified by dataindex.
  double    etotal(int dataindex);
  
  /// Return the count for the given data file and the given cut-point.
  double    count(int dataindex, int cutindex);

  /// Return all variables read from the cut file(s).
  vstring&  vars();

  /// Return number of cuts
  int       ncuts();

  /// Return cut values for cut-point identified by cutindex.
  vdouble   cuts(int cutindex);

  /// Return cut variable names
  vstring   cutvars();

  /// Return number of events for specified data file.
  int       ndata(int dataindex);

  /// Return values for data given data file and event.
  vdouble&  data(int dataindex, int event);

  /// Save counts to a text or root file
  void     save(std::string resultfilename);
  
private:
  int         _status;

  vvdouble    _cutdata;
  varmap      _varmap;
  vstring     _var;
  vint        _cutcode;
  vint        _cutpointcount;

  std::string _treename;
  std::string _weightname;
  std::string _selection;
  
  vstring                   _searchname;
  std::vector< vvdouble >   _searchdata;
  vint                      _weightindex;
  vint        _index;
  vdouble     _weight;  
  vdouble     _totals;
  vdouble     _errors;
  vvdouble    _counts;
  vstring _resultname;
  
  std::vector<std::vector<int> > _cutpointindex;

  void _init(vstring& cutdatafilename, int start=0, int numrows=0, 
	     std::string treename="",
	     std::string selection="");
  bool _twosidedcut(float x, int cutpoint, int jcut);
  bool _staircasecut(vdouble& datarow, int cutpoint, int& cut);
  void _saveToTextFile(std::string resultfilename);
  void _saveToNtupleFile(std::string resultfilename);
};

#endif
