/******************************************************************************
 *  Tagger: Run-length compressed metagenomic read classification with        *
 *          SMEM-finding and tagging                                          *
 *  Copyright (C) 2025 - Lore Depuydt <lore.depuydt@ugent.be> and             *
 *                       Luca Renders <luca.renders@ugent.be> and             *
 *                       Jan Fostier <jan.fostier@ugent.be>                   *
 *                                                                            *
 *  This program is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU Affero General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  This program is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Affero General Public License for more details.                       *
 *                                                                            *
 * You should have received a copy of the GNU Affero General Public License   *
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.     *
 ******************************************************************************/

#include "bmove.h"
#include "../alphabet.h" // for Alphabet
#include "../definitions.h"
#include "../indexhelpers.h"
#include "../logger.h" // for Logger
#include "../reads.h"

#include <algorithm>           // for max
#include <assert.h>            // for assert
#include <cstdint>             // for uint16_t
#include <ostream>             // for opera...
#include <sdsl/int_vector.hpp> // for int_v...
#include <stdexcept>           // for runti...
#include <string>              // for string
class MemoryMappedTextFile;
class Search;
class Substring;

using namespace std;

// ----------------------------------------------------------------------------
// PREPROCESSING ROUTINES
// ----------------------------------------------------------------------------

void BMove::fromFiles(const string& baseFile, bool verbose
#ifdef TAG_ARRAY_SUBSAMPLING
                      ,
                      length_t maxLF
#endif // TAG_ARRAY_SUBSAMPLING
) {

    readMetaAndCounts(baseFile, verbose);
    stringstream ss;
    // Read BMove specific files

    string LFextension = ".LFBP";

    if (verbose) {
        ss << "Reading " << baseFile << LFextension << "...";
        logger.logInfo(ss);
    }
    if (!move.load(baseFile)) {
        throw runtime_error("Error loading move file: " + baseFile +
                            LFextension);
    }

    if (verbose) {
        ss << "Reading " << baseFile << ".rev" << LFextension << "...";
        logger.logInfo(ss);
    }
    if (!moveR.load(baseFile + ".rev")) {
        throw runtime_error("Error loading reverse move file: " + baseFile +
                            ".rev" + LFextension);
    }

#ifdef TAG_ARRAY_SUBSAMPLING
    if (verbose) {
        ss << "Reading " << baseFile << ".tag." << maxLF << ".bv.bwt...";
        logger.logInfo(ss);
    }
    if (!samplesTaggingBV.read(baseFile + ".tag." + to_string(maxLF) +
                               ".bv.bwt")) {
        throw runtime_error("Cannot open file: " + baseFile + ".tag." +
                            to_string(maxLF) + ".bv.bwt");
    }

    maxLFSteps = maxLF;
#else

    length_t maxLF = 0;

#endif // TAG_ARRAY_SUBSAMPLING

    if (verbose) {
        ss << "Reading " << baseFile << ".tag." << maxLF << ".heads.bwt...";
        logger.logInfo(ss);
    }
    if (!readIntVector(baseFile + ".tag." + to_string(maxLF) + ".heads.bwt",
                       samplesTagging)) {
        throw runtime_error("Cannot open file: " + baseFile + ".tag." +
                            to_string(maxLF) + ".heads.bwt");
    }
    initialToeholdTagging =
        getTaggingSample(getTextSize() - 1, move.size() - 1);
}

// ----------------------------------------------------------------------------
// ROUTINES FOR ACCESSING DATA STRUCTURE
// ----------------------------------------------------------------------------

length_t BMove::getTaggingSample(length_t start_pos, length_t start_run) const {
#ifdef TAG_ARRAY_SUBSAMPLING
    size_t num_steps = 0;
    bool not_last_index_in_run = false;

    while (not_last_index_in_run || samplesTaggingBV[start_run] == 0) {
        move.findLF(start_pos, start_run);
        assert(start_run < move.size());
        not_last_index_in_run =
            move.getInputStartPos(start_run + 1) - 1 != start_pos;
        num_steps++;
        assert(num_steps <= maxLFSteps);
    }

    length_t sample_index = samplesTaggingBV.rank(start_run);

    // Assert that sample_index is within bounds for samplesTagging
    // assert(sample_index < samplesTagging.size());

    return samplesTagging[sample_index];
#else
    return samplesTagging[start_run];
#endif
}

length_t BMove::computeToeholdTagging(const MoveRange& range, const length_t c,
                                      const length_t previousToehold) const {
    length_t endRun = range.getEndRun();

    if (move.runHeadEquals(endRun, c)) {
        return previousToehold;
    }

    length_t previousPos;
    length_t previousRun;

    move.walkToPreviousRun(range, previousPos, previousRun, c);

    return getTaggingSample(previousPos, previousRun);
}

// ----------------------------------------------------------------------------
//  APPROXIMATE PATTERN MATCHING
// ----------------------------------------------------------------------------

void BMove::updateRangeSARevRuns(SARangePair& ranges) const {
    moveR.computeRunIndices(ranges.getRangeSARevMutable());
}

void BMove::findRangeWithExtraCharBackwardAuxiliary(
    length_t positionInAlphabet, SARange& parentBackwardRange,
    SARange& childBackwardRange) const {

    if (!parentBackwardRange.getRunIndicesValid()) {
        move.computeRunIndices(parentBackwardRange);
    }
    move.addChar(parentBackwardRange, childBackwardRange, positionInAlphabet);
}

bool BMove::findRangeWithExtraCharBackward(length_t posInAlpha,
                                           const SARangeBackwards& rangeOfP,
                                           SARangeBackwards& childRange) const {

    // first make the backward range by searching cP using B
    SARange range1, trivialRange = rangeOfP;
    findRangeWithExtraCharBackwardAuxiliary(posInAlpha, trivialRange, range1);

    if (range1.empty()) {
        childRange = SARangeBackwards(range1, 0, false, 0);
        return false;
    }

    // Check if range1 has the same width as the range of the parent
    if (trivialRange.width() == range1.width()) {
        childRange = SARangeBackwards(
            range1, rangeOfP.getToehold() - !rangeOfP.getToeholdRepresentsEnd(),
            rangeOfP.getToeholdRepresentsEnd(),
            rangeOfP.getOriginalDepth() + 1);
        return true;
    }

    length_t newToehold = 0;

    childRange = SARangeBackwards(range1, newToehold, false,
                                  rangeOfP.getOriginalDepth() + 1);
    return true;
}

bool BMove::findRangesWithExtraCharBackward(length_t positionInAlphabet,
                                            const SARangePair& rangesOfP,
                                            SARangePair& childRanges) const {

    // first make the backward range by searching cP using B
    SARange range1, trivialRange = rangesOfP.getRangeSA();
    findRangeWithExtraCharBackwardAuxiliary(positionInAlphabet, trivialRange,
                                            range1);

    // if the range is empty, return false
    if (range1.empty()) {
        childRanges = SARangePair(range1, range1, 0, false, 0);
        return false;
    }

    // Check if range1 has the same width as the range of the parent
    const SARange& otherRange = rangesOfP.getRangeSARev();
    if (trivialRange.width() == range1.width()) {
        // otherRange remains the same
        childRanges = SARangePair(range1, otherRange,
                                  rangesOfP.getToehold() -
                                      !rangesOfP.getToeholdRepresentsEnd(),
                                  rangesOfP.getToeholdRepresentsEnd(),
                                  rangesOfP.getOriginalDepth() + 1);
        return true;
    }

    // then make the less trivial range by counting the sizes of the ranges
    // of (dP) using B

    // first get the start of the range we are looking for of the parent
    length_t s = rangesOfP.getRangeSARev().getBegin();

    // get the start of the child within this range
    // find the number of occurrences of chars smaller than c in the parent
    // range
    length_t x = move.getCumulativeCounts(trivialRange, positionInAlphabet);

    // make the new range with width equal to that of the trivial range
    SARange range2 = SARange(s + x, s + x + range1.width(),
                             otherRange.getBeginRun(), otherRange.getEndRun());
    // set the run indices to be invalid
    range2.setRunIndicesValid(false);

    // aP occurs for some a != c, the toehold must be updated
    length_t newToehold = computeToeholdTagging(
        trivialRange, positionInAlphabet, rangesOfP.getToehold());

    // set the final SARangePair
    childRanges = SARangePair(range1, range2, newToehold, false,
                              rangesOfP.getOriginalDepth() + 1);
    return true;
}

bool BMove::findRangesWithExtraCharForward(length_t positionInAlphabet,
                                           const SARangePair& rangesOfP,
                                           SARangePair& childRanges) const {

    // first make the trivial range by searching (Pc)' using B'  if
    // searching forward we need to use B' so we need the reverse range
    SARange trivialRange = rangesOfP.getRangeSARev();

    // If the run indices of trivialRange are not valid (e.g., after a
    // directions switch), compute them
    if (!trivialRange.getRunIndicesValid()) {
        moveR.computeRunIndices(trivialRange);
    }

    // get the range of the child by adding ona character using move
    SARange range1;
    moveR.addChar(trivialRange, range1, positionInAlphabet);

    // if the range is empty, return false
    if (range1.empty()) {
        childRanges = SARangePair(range1, range1, 0, false, 0);
        return false;
    }

    // Check if range1 has the same width as the range of the parent
    const SARange& otherRange = rangesOfP.getRangeSA();
    if (trivialRange.width() == range1.width()) {
        // otherRange remains the same
        childRanges = SARangePair(otherRange, range1,
                                  rangesOfP.getToehold() +
                                      rangesOfP.getToeholdRepresentsEnd(),
                                  rangesOfP.getToeholdRepresentsEnd(),
                                  rangesOfP.getOriginalDepth() + 1);
        return true;
    }

    // then make the less trivial range by counting the size of the range of
    // (Pd)' using B' (forward)

    // first get the start of the range we are looking for of the parent
    length_t s = rangesOfP.getRangeSA().getBegin();

    // get the start of the child within this range
    // find the number of occurrences of chars smaller than c in the parent
    // range

    length_t x = moveR.getCumulativeCounts(trivialRange, positionInAlphabet);

    // make the new range with width equal to that of the trivial range
    SARange range2 = SARange(s + x, s + x + range1.width(),
                             otherRange.getBeginRun(), otherRange.getEndRun());
    // set the run indices to be invalid
    range2.setRunIndicesValid(false);

    length_t newToehold = 0;

    // set the final SARangePair
    childRanges = SARangePair(range2, range1, newToehold, true,
                              rangesOfP.getOriginalDepth() + 1);
    return true;
}

// ----------------------------------------------------------------------------
// TAGGING ROUTINES
// ----------------------------------------------------------------------------

bool BMove::taggingFunction(
    const std::string& readID, const std::string& read,
    std::vector<TextOcc>& occs, Counters& counters, const Strand strand,
    const uint16_t correctTag, const length_t minMEMLength, bool& allMEMsunique,
    std::unordered_map<uint16_t, length_t>& detectedTags,
    std::unordered_map<uint16_t, int64_t>& lastSeenPositions,
    length_t& totalWeight) const {

    bool mapped = false;

    assert(wordSize <= minMEMLength);

    length_t stepCount = 0;
    int64_t readSize = read.size();
    int64_t j = readSize;
    int64_t i = j;

    while (j >= static_cast<int64_t>(minMEMLength)) {
        int64_t k = j;
        if (i > j - static_cast<int64_t>(minMEMLength)) {
            int64_t jumpStart = j - minMEMLength;
            k = findMEMEndWithLookup(read, jumpStart, j, stepCount);
        }

        if (k < j) {
            j = k;
        } else {
            assert(j == k);

            mapped = true;

            SARange range;
            length_t toeHold;
            bool encounteredN = false;

            i = findMEMStartWithLookup(read, j, range, stepCount, encounteredN,
                                       readSize, toeHold);
            assert(static_cast<int64_t>(i) <=
                   j - static_cast<int64_t>(minMEMLength));

            processMatches(read, i, j, counters, occs, range, strand,
                           correctTag, readID, allMEMsunique, detectedTags,
                           lastSeenPositions, totalWeight, toeHold);

            if (encounteredN) {
                j = skipNs(read, i, readSize);
                continue;
            }

            if (i == 0) {
                break;
            }

            j = findMEMEnd(read, i - 1, j, stepCount);
            i--;
        }
    }

    counters.inc(Counters::TOTAL_STEP_COUNTS, stepCount);

    return mapped;
}

int64_t BMove::findMEMStartWithLookup(const std::string& read, int64_t end,
                                      SARange& range, length_t& stepCount,
                                      bool& encounteredN, int64_t readSize,
                                      length_t& toeHold) const {

    assert(end >= static_cast<int64_t>(wordSize));

    SARangePair ranges =
        lookUpInKmerTable({read, static_cast<length_t>(end - wordSize),
                           static_cast<length_t>(end)});

    range = ranges.getRangeSA();
    toeHold = ranges.getToehold();

    assert(!range.empty());

    int64_t k = end - 1 - wordSize;
    int64_t start = end - wordSize;
    for (; k >= 0; --k) {
        char c = read[k];
        if (c == 'N') {
            encounteredN = true;
            break;
        }
        length_t cIdx = sigma.c2i(c);
        SARange newRange;
        stepCount++;
        move.addChar(range, newRange, cIdx);
        if (newRange.empty()) {
            break;
        }
        toeHold = computeToeholdTagging(range, cIdx, toeHold);
        range = newRange;
        start = k;
    }
    return start;
}

int64_t BMove::findMEMEndWithLookup(const std::string& read, int64_t i,
                                    int64_t j, length_t& stepCount) const {

    assert(static_cast<int64_t>(read.size()) - i >=
           static_cast<int64_t>(wordSize));

    SARangePair ranges = lookUpInKmerTable(
        {read, static_cast<length_t>(i), static_cast<length_t>(i + wordSize)});

    SARange range = ranges.getRangeSARev();

    if (range.empty()) {
        return j - 1; // TODO: Do we lose a lot of time here?
    }

    assert(i >= 0);
    int64_t k = i + wordSize;
    for (; k < j; ++k) {
        char c = read[k];
        length_t cIdx = sigma.c2i(c);
        if (cIdx > sigma.size()) {
            break;
        }
        SARange newRange;
        stepCount++;
        moveR.addChar(range, newRange, cIdx);
        if (newRange.empty()) {
            break;
        }
        range = newRange;
    }
    return k;
}

int64_t BMove::findMEMEnd(const std::string& read, int64_t i, int64_t j,
                          length_t& stepCount) const {
    SARange range = getCompleteRangeForward();
    assert(i >= 0);
    int64_t k = i;
    for (; k < j; ++k) {
        char c = read[k];
        length_t cIdx = sigma.c2i(c);
        SARange newRange;
        stepCount++;
        moveR.addChar(range, newRange, cIdx);
        if (newRange.empty()) {
            break;
        }
        range = newRange;
    }
    return k;
}

uint16_t BMove::computeTagDetails(const Range& range,
                                  length_t& numberOfTagsBefore,
                                  length_t& numberOfTagsAfter,
                                  length_t& toeHold) const {
    return toeHold;
}

void BMove::updateDetectedTags(
    Counters& counters, int length, uint16_t detectedTag,
    const uint16_t correctTag,
    std::unordered_map<uint16_t, length_t>& detectedTags,
    length_t& totalWeight) const {
    detectedTags[detectedTag] += length;
    totalWeight += length;
}

uint16_t
BMove::computeTagMEM(Counters& counters, int64_t i, int64_t j,
                     const Range& range, const uint16_t correctTag,
                     bool& allMEMsunique,
                     std::unordered_map<uint16_t, length_t>& detectedTags,
                     std::unordered_map<uint16_t, int64_t>& lastSeenPositions,
                     length_t& totalWeight, length_t& toeHold) const {
    int64_t length = j - i;

    length_t numberOfTagsBefore = 0;
    length_t numberOfTagsAfter = 0;

    // Use helper to compute tag details
    uint16_t detectedTag = computeTagDetails(range, numberOfTagsBefore,
                                             numberOfTagsAfter, toeHold);

    if (detectedTag == UINT16_MAX) {
        // No unique tag
        allMEMsunique = false;
        return 0;
    }
    // Unique match
    // Check if the tag was already seen
    auto it = lastSeenPositions.find(detectedTag);
    if (it != lastSeenPositions.end()) {
        int64_t lastSeenPosition = it->second;
        length = min(j, lastSeenPosition) - i; // Update length if seen before
    }

    // Update detected tags and counters
    updateDetectedTags(counters, length, detectedTag, correctTag, detectedTags,
                       totalWeight);
    lastSeenPositions[detectedTag] = i;
    return detectedTag;
}

void BMove::processMatches(
    const std::string& read, int64_t i, int64_t j, Counters& counters,
    std::vector<TextOcc>& occs, const SARange& range, const Strand strand,
    const uint16_t correctTag, const std::string& readID, bool& allMEMsunique,
    std::unordered_map<uint16_t, length_t>& detectedTags,
    std::unordered_map<uint16_t, int64_t>& lastSeenPositions,
    length_t& totalWeight, length_t& toeHold) const {
    std::vector<TextOcc> temp;
    assert(exactMatchesOutput(read.substr(i, j - i), counters, temp));
    assert(i == 0 ||
           !exactMatchesOutput(read.substr(i - 1, j - i + 1), counters, temp));
    assert(static_cast<size_t>(j) == read.size() ||
           !exactMatchesOutput(read.substr(i, j - i + 1), counters, temp));

    // Compute the tag and update the necessary counters
    computeTagMEM(counters, i, j, range, correctTag, allMEMsunique,
                  detectedTags, lastSeenPositions, totalWeight, toeHold);

    counters.inc(Counters::TOTAL_UNIQUE_MATCHES);
}

int64_t BMove::skipNs(const std::string& read, int64_t pos,
                      int64_t readSize) const {
    pos--;
    while (pos >= 0 && read[pos] == 'N') {
        pos--;
    }
    pos++;
    return pos;
}