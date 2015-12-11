/*****************************************************************************
 my_coverage.cpp

 Code modified from: https://github.com/ryanlayer/bedtools/blob/323b65f37ffb846b4216b0939700a0cb5f41182e/src/coverageBed/coverageBed.h
 Date of modification: 2015/11/04

 Licensed under the GNU General Public License 2.0 license.
 ******************************************************************************/


#include "my_coverage.h"
using namespace std;

BedCoverage::BedCoverage(string &bedAFile, string &bedBFile, bool sameStrand,
                         bool diffStrand, bool obeySplits, bool eachBase, bool countsOnly)
{

  _bedAFile       = bedAFile;
  _bedBFile       = bedBFile;

  _bedA           = new BedFile(bedAFile);
  _bedB           = new BedFile(bedBFile);

  _sameStrand     = sameStrand;
  _diffStrand     = diffStrand;
  _obeySplits     = obeySplits;
  _eachBase       = eachBase;
  _countsOnly     = countsOnly;

}

// destroy
BedCoverage::~BedCoverage(void) {
  delete _bedA;
  delete _bedB;
}

void BedCoverage::CollectCoverageBed() {

  // load the "B" bed file into a map so
  // that we can easily compare "A" to it for overlaps
  _bedB->loadBedCovFileIntoMap();


  BED a;
  _bedA->Open();
  // process each entry in A
  while (_bedA->GetNextBed(a)) {
    if (_bedA->_status == BED_VALID) {
      // process the BED entry as a single block
      if (_obeySplits == false)
        _bedB->countHits(a, _sameStrand,
                         _diffStrand, _countsOnly);
      // split the BED into discrete blocksand process
      // each independently.
      else {
        bedVector bedBlocks;
        GetBedBlocks(a, bedBlocks);
        // use countSplitHits to avoid over-counting
        // each split chunk as distinct read coverage.
        _bedB->countSplitHits(bedBlocks, _sameStrand,
                              _diffStrand, _countsOnly);
      }
    }
  }
  _bedA->Close();

  // report the coverage (summary or histogram) for BED B.
  if (_countsOnly == true)
    ReportCounts();
  else
    ReportCoverage();
}

void BedCoverage::ReportCounts() {

  List counts;

  // process each chromosome
  masterBedCovMap::const_iterator chromItr = _bedB->bedCovMap.begin();
  masterBedCovMap::const_iterator chromEnd = _bedB->bedCovMap.end();
  for (; chromItr != chromEnd; ++chromItr)
  {
    // for each chrom, process each bin
    binsToBedCovs::const_iterator binItr = chromItr->second.begin();
    binsToBedCovs::const_iterator binEnd = chromItr->second.end();
    for (; binItr != binEnd; ++binItr)
    {
      // for each chrom & bin, compute and report
      // the observed coverage for each feature
      vector<BEDCOV>::const_iterator bedItr = binItr->second.begin();
      vector<BEDCOV>::const_iterator bedEnd = binItr->second.end();
      for (; bedItr != bedEnd; ++bedItr)
      {
        _bedB->reportBedTab(*bedItr);
        Rprintf("%d\n", bedItr->count);
      }
    }
  }
}

void BedCoverage::ReportCoverage() {

  map<unsigned long, unsigned long> allDepthHist;
  unsigned long totalLength = 0;

  // process each chromosome
  masterBedCovMap::const_iterator chromItr = _bedB->bedCovMap.begin();
  masterBedCovMap::const_iterator chromEnd = _bedB->bedCovMap.end();
  for (; chromItr != chromEnd; ++chromItr)
  {
    // for each chrom, process each bin
    binsToBedCovs::const_iterator binItr = chromItr->second.begin();
    binsToBedCovs::const_iterator binEnd = chromItr->second.end();
    for (; binItr != binEnd; ++binItr)
    {
      // for each chrom & bin, compute and report
      // the observed coverage for each feature
      vector<BEDCOV>::const_iterator bedItr = binItr->second.begin();
      vector<BEDCOV>::const_iterator bedEnd = binItr->second.end();
      for (; bedItr != bedEnd; ++bedItr)
      {
        int zeroDepthCount = 0; // number of bases with zero depth
        int depth          = 0; // tracks the depth at the current base

        // the start is either the first base in the feature OR
        // the leftmost position of an overlapping feature.
        // e.g. (s = start):
        // A    ----------
        // B    s    ------------
        int start = min(bedItr->minOverlapStart, bedItr->start);

        // track the number of bases in the feature covered by
        // 0, 1, 2, ... n features in A
        map<unsigned int, unsigned int> depthHist;
        map<unsigned int, DEPTH>::const_iterator depthItr;

        // compute the coverage observed at each base in
        // the feature marching from start to end.
        for (CHRPOS pos = start+1; pos <= bedItr->end; pos++)
        {
          // map pointer grabbing the starts and
          // ends observed at this position
          depthItr = bedItr->depthMap.find(pos);
          // increment coverage if starts observed at this position.
          if (depthItr != bedItr->depthMap.end())
            depth += depthItr->second.starts;
          // update coverage assuming the current position is
          // within the current B feature
          if ((pos > bedItr->start) && (pos <= bedItr->end)) {
            if (depth == 0) zeroDepthCount++;
            // update our histograms, assuming we are not
            // reporting "per-base" coverage.
            if (_eachBase == false) {
              depthHist[depth]++;
              allDepthHist[depth]++;
            }
            else if ((_eachBase == true) &&
                     (bedItr->zeroLength == false))
            {
              _bedB->reportBedTab(*bedItr);
              Rprintf("%d\t%d\n", pos-bedItr->start, depth);
            }
          }
          // decrement coverage if ends observed at this position.
          if (depthItr != bedItr->depthMap.end())
            depth = depth - depthItr->second.ends;
        }

        // handle the special case where the user wants "per-base" depth
        // but the current feature is length = 0.
        if ((_eachBase == true) && (bedItr->zeroLength == true)) {
          _bedB->reportBedTab(*bedItr);
          Rprintf("1\t%d\n",depth);
        }
        // Summarize the coverage for the current interval,
        // assuming the user has not requested "per-base" coverage.
        else if (_eachBase == false)
        {
          CHRPOS length     = bedItr->end - bedItr->start;
          if (bedItr->zeroLength == true) {
            length = 0;
          }
          totalLength       += length;
          int nonZeroBases   = (length - zeroDepthCount);

          float fractCovered = 0.0;
          if (bedItr->zeroLength == false) {
            fractCovered = (float) nonZeroBases / length;
          }

          // print a summary of the coverage
          _bedB->reportBedTab(*bedItr);
          Rprintf("%d\t%d\t%d\t%0.7f\n",
                 bedItr->count, nonZeroBases,
                 length, fractCovered);
        }
      }
    }
  }
}
