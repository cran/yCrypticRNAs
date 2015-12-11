/*****************************************************************************
  BlockedIntervals.cpp

  Code modified from https://github.com/arq5x/bedtools2/blob/master/src/utils/BlockedIntervals/BlockedIntervals.cpp
 
  Date of modification: 2015/11/04

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#include "BlockedIntervals.h"

void GetBedBlocks(const BED &bed, bedVector &bedBlocks) {

    // nothing to do if it is not a BED12 intervals
    if (bed.fields.size() != 12) {
        bedBlocks.push_back(bed);
        return;
    }

    int blockCount = atoi(bed.fields[9].c_str());
    if ( blockCount <= 0 ) {
        Rcpp::stop("Input error: found interval having <= 0 blocks.\n");
    }
    else if ( blockCount == 1 ) {
        //take a short-cut for single blocks
        bedBlocks.push_back(bed);
    }
    else {
        // get the comma-delimited strings for the BED12 block starts and block ends.
        string blockSizes(bed.fields[10]);
        string blockStarts(bed.fields[11]);

        vector<int> sizes;
        vector<int> starts;
        Tokenize(blockSizes, sizes, ',');
        Tokenize(blockStarts, starts, ',');

        if ( sizes.size() != (size_t) blockCount || starts.size() != (size_t) blockCount ) {
            Rcpp::stop("Input error: found interval with block-counts not matching starts/sizes on line. \n");
        }

        // add each BED block to the bedBlocks vector
        for (UINT i = 0; i < (UINT) blockCount; ++i) {
            CHRPOS blockStart = bed.start + starts[i];
            CHRPOS blockEnd   = bed.start + starts[i] + sizes[i];
            BED currBedBlock(bed.chrom, blockStart, blockEnd,
                             bed.name, bed.score, bed.strand, bed.fields, bed.other_idxs);
            bedBlocks.push_back(currBedBlock);
        }
    }
}

int GetTotalBlockLength(const bedVector &bedBlocks) {
    int total_size = 0;
    bedVector::const_iterator blockItr = bedBlocks.begin();
    bedVector::const_iterator blockEnd = bedBlocks.end();
    for (; blockItr != blockEnd; ++blockItr) {
        total_size += blockItr->end - blockItr->start;
    }
    return total_size;
}
