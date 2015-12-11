/*****************************************************************************
  my_coverage.h

  Code modified from: https://github.com/ryanlayer/bedtools/blob/323b65f37ffb846b4216b0939700a0cb5f41182e/src/coverageBed/coverageBed.h
  Date of modification: 2015/11/04

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef COVERAGEBED_H
#define COVERAGEBED_H

#include "bedFile.h"

#include "BlockedIntervals.h"

#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <list>

using namespace std;

//************************************************
// Class methods and elements
//************************************************
class BedCoverage {

public:

    // constructor
    BedCoverage(string &bedAFile, string &bedBFile, bool sameStrand, bool diffStrand,
                bool obeySplits, bool eachBase, bool countsOnly);

    // destructor
    ~BedCoverage(void);

    void CollectCoverageBed();

private:

    // input files.
    string _bedAFile;
    string _bedBFile;

    // instance of a bed file class.
    BedFile *_bedA, *_bedB;

    // do we care about same or opposite strandedness when counting coverage?
    bool _sameStrand;
    bool _diffStrand;

    // should we split BED/BAM into discrete blocks?
    bool _obeySplits;

    // should discrete coverage be reported for each base in each feature?
    bool _eachBase;

    // should we just count overlaps and not try to describe the breadth?
    bool _countsOnly;

    // private function for reporting coverage information
    void ReportCoverage();

    // private function for reporting overlap counts
    void ReportCounts();
};
#endif /* COVERAGEBED_H */
