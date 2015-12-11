/*****************************************************************************
  my_coverageMain.cpp

  Code modified from: https://github.com/ryanlayer/bedtools/blob/323b65f37ffb846b4216b0939700a0cb5f41182e/src/coverageBed/coverageMain.cpp

  Date of modification: 2015/11/04

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/

#include <RcppCommon.h>
#include "my_coverage.h"

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// function declarations
void coverage_help(void);

// [[Rcpp::export]]
void coverage_cpp(std::string bedAFile, std::string bedBFile, CharacterVector options = "") {

  // our configuration variables
  bool showHelp = false;

  // parm flags
  bool haveBedA       = false;
  bool haveBedB       = false;
  bool sameStrand     = false;
  bool diffStrand     = false;
  bool eachBase       = false;
  bool obeySplits     = false;
  bool countsOnly     = false;

  if(!bedAFile.empty()){
    haveBedA = true;
  }
  if(!bedBFile.empty()){
    haveBedB = true;
  }

  CharacterVector::iterator start = options.begin();
  CharacterVector::iterator end = options.end();

  if(std::find(start, end, "-s") != end){
    sameStrand = true;
  }
  if(std::find(start, end, "-S") != end) {
    diffStrand = true;
  }
  if(std::find(start, end, "-d") != end) {
    eachBase = true;
  }
  if(std::find(start, end, "-split") != end) {
    obeySplits = true;
  }
  if(std::find(start, end, "-counts") != end) {
    countsOnly = true;
  }

  // make sure we have both input files
  if (!haveBedA || !haveBedB) {
    Rcerr << endl << "*****" << endl << "*****ERROR: Need two inputs bed files. " << endl << "*****" << endl;
    showHelp = true;
  }

  if (sameStrand && diffStrand) {
    Rcerr << endl << "*****" << endl << "*****ERROR: Request either -s OR -S, not both." << endl << "*****" << endl;
    showHelp = true;
  }

  if (!showHelp) {
    BedCoverage *bg = new BedCoverage(bedAFile, bedBFile, sameStrand, diffStrand,
                                     obeySplits, eachBase, countsOnly);
    bg->CollectCoverageBed();
    delete bg;
  }
  else {
    coverage_help();
  }
}

void coverage_help(void) {

  Rcerr << "\nTool:    bedtools coverage (aka coverageBed)" << endl;
  Rcerr << "Summary: Returns the depth and breadth of coverage of features from A" << endl;
  Rcerr << "\t on the intervals in B." << endl << endl;

  Rcerr << "Usage:   " << "coverage_cpp( std::string bedAFile, std::string bedBFile, CharacterVector options = "") "<< endl << endl;

  Rcerr << "Options: " << endl;

  Rcerr << "\t-s\t"            << "Require same strandedness.  That is, only counts hits in A that" << endl;
  Rcerr                        << "\t\toverlap B on the _same_ strand." << endl;
  Rcerr                        << "\t\t- By default, overlaps are counted without respect to strand." << endl << endl;

  Rcerr << "\t-S\t"            << "Require different strandedness.  That is, only report hits in A" << endl;
  Rcerr                        << "\t\tthat overlap B on the _opposite_ strand." << endl;
  Rcerr                        << "\t\t- By default, overlaps are counted without respect to strand." << endl << endl;

  Rcerr << "\t-d\t"            << "Report the depth at each position in each B feature." << endl;
  Rcerr                        << "\t\tPositions reported are one based.  Each position" << endl;
  Rcerr                        << "\t\tand depth follow the complete B feature." << endl << endl;

  Rcerr << "\t-counts\t"       << "Only report the count of overlaps, don't compute fraction, etc." << endl << endl;

  Rcerr << "Default Output:  " << endl;
  Rcerr << "\t" << " After each entry in B, reports: " << endl;
  Rcerr << "\t   1) The number of features in A that overlapped the B interval." << endl;
  Rcerr << "\t   2) The number of bases in B that had non-zero coverage." << endl;
  Rcerr << "\t   3) The length of the entry in B." << endl;
  Rcerr << "\t   4) The fraction of bases in B that had non-zero coverage." << endl << endl;

  Rcpp::stop("error");
}
