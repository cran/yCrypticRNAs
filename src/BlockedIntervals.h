/*****************************************************************************
  BlockedIntervals.h

  Code modified from https://github.com/arq5x/bedtools2/blob/master/src/utils/BlockedIntervals/BlockedIntervals.h
 
  Date of modification: 2015/11/04

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#include <vector>
#include "bedFile.h"

using namespace std;

/* break a BED12 record into discrete BED6 blocks. */
void GetBedBlocks(const BED &bed, bedVector &bedBlocks);

/* compute the total forprint of a set of BED blocks */
int GetTotalBlockLength(const bedVector &bedBlocks);
