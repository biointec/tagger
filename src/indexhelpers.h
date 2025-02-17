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

#ifndef INDEX_HELPERS_H
#define INDEX_HELPERS_H

#include "definitions.h" // for length_t, Strand, FIRST_IN...
#include "reads.h"       // for ReadBundle

#include <algorithm>          // for max, copy, sort, unique
#include <array>              // for array
#include <assert.h>           // for assert
#include <cmath>              // for log10, round
#include <cstdint>            // for uint16_t, uint32_t, uint64_t, int16_t
#include <ext/alloc_traits.h> // for __alloc_traits<>::value_type
#include <fmt/core.h>         // for format
#include <fmt/format.h>       // for to_string
#include <memory>             // for allocator, allocator_traits<>::value_type
#include <sstream>            // for ostream
#include <stddef.h>           // for size_t
#include <string>             // for string, operator+, char_traits, basic_...
#include <unordered_map>      // for unordered_map
#include <utility>            // for move, pair
#include <vector>             // for vector

// ============================================================================
// HELPERS
// ============================================================================

// compute |a-b| in a safe manner
template <typename T> T abs_diff(T a, T b) {
    return a > b ? a - b : b - a;
}

// ============================================================================
// CLASS RANGE
// ============================================================================

/**
 * A range. A range is a pair of integers (begin, end) where begin is the
 * beginning of the range and end is the end of the range (non-inclusive).
 */
class Range {
  protected:
    length_t begin; // beginning of the range
    length_t end;   // end of the range (non-inclusive)

  public:
    /**
     * Constructor
     * @param b the beginning of the range
     * @param e  the end of the range (non-inclusive)
     */
    Range(length_t b, length_t e) : begin(b), end(e) {
    }

    /**
     * Default constructor, initializes an empty range
     */
    Range() : begin(0), end(0) {
    }

    /**
     * @returns the beginning of the range
     */
    length_t getBegin() const {
        return begin;
    }

    /**
     * @returns the end of the range (non-inclusive)
     */
    length_t getEnd() const {
        return end;
    }
    /**
     * Check if this range is empty.
     * @returns true if the range is empty, false otherwise
     */
    bool empty() const {
        return end <= begin;
    }

    /**
     * Gets the width of the range (end - begin)
     * @returns the width of this range
     */
    length_t width() const {
        return (empty()) ? 0 : end - begin;
    }

    /**
     * Operator overloading, two ranges are equal if their begin and end field
     * are equal.
     * @param o the range to compare to this
     * @returns true if this is equal to o
     */
    bool operator==(const Range& o) const {
        return o.getBegin() == begin && o.getEnd() == end;
    }
    /**
     * Operator overloading, to stream the range to an output stream.
     * @param os the output stream
     * @param r the range to print
     */
    friend std::ostream& operator<<(std::ostream& os, const Range& r);
};

std::ostream& operator<<(std::ostream& output, const Range& r);

// ============================================================================
// CLASS MOVERANGE
// ============================================================================

class MoveRange : public Range {
  private:
    length_t beginRun; // run index corresponding to the beginning of the range
    length_t
        endRun; // run index corresponding to the end of the range (inclusive)
    bool runIndicesValid;

  public:
    /**
     * Constructor
     * @param b, the beginning of the range
     * @param e, the end of the range (non-inclusive)
     * @param bRun, the run index corresponding to the beginning of the range
     * @param eRun, the run index corresponding to the end of the range
     * (inclusive)
     * @param runIndicesValid, true if the run indices are valid, false
     * otherwise
     */
    MoveRange(length_t b, length_t e, length_t bRun, length_t eRun,
              bool runIndicesValid = true)
        : Range(b, e), beginRun(bRun), endRun(eRun),
          runIndicesValid(runIndicesValid) {
    }

    /**
     * Default constructor, initializes an empty range
     */
    MoveRange() : Range(), beginRun(0), endRun(0), runIndicesValid(false) {
    }

    /**
     * @returns the run index corresponding to the beginning of the range
     */
    length_t getBeginRun() const {
        return beginRun;
    }

    /**
     * @returns the run index corresponding to the end of the range (inclusive)
     */
    length_t getEndRun() const {
        return endRun;
    }

    /**
     * Sets the run index corresponding to the beginning of the range
     * @param bRun the run index to set
     */
    void setBeginRun(length_t bRun) {
        beginRun = bRun;
    }

    /**
     * Sets the run index corresponding to the end of the range
     * @param eRun the run index to set
     */
    void setEndRun(length_t eRun) {
        endRun = eRun;
    }

    /**
     * @brief Set the Begin object
     *
     * @param b
     */
    void setBegin(length_t b) {
        begin = b;
    }

    /**
     * @brief Set the End object
     *
     * @param e
     */
    void setEnd(length_t e) {
        end = e;
    }

    /**
     * @brief Set the range to empty values
     *
     */
    void setEmpty() {
        begin = 0;
        end = 0;
        beginRun = 0;
        endRun = 0;
        runIndicesValid = false;
    }

    /**
     * @returns true if the run indices are valid, false otherwise
     */
    bool getRunIndicesValid() const {
        return runIndicesValid;
    }

    /**
     * Sets the run indices to valid or invalid
     * @param valid true if the run indices are valid, false otherwise
     */
    void setRunIndicesValid(bool valid) {
        this->runIndicesValid = valid;
    }

    /**
     * Operator overloading, two MoveRanges are equal if their begin, end,
     * beginRun and endRun fields are equal.
     * @param o the MoveRange to compare to this
     * @returns true if this is equal to o
     */
    bool operator==(const MoveRange& o) const {
        return o.getBegin() == begin && o.getEnd() == end &&
               o.getBeginRun() == beginRun && o.getEndRun() == endRun &&
               o.getRunIndicesValid() == runIndicesValid;
    }
    /**
     * Operator overloading, to stream the MoveRange to an output stream.
     * @param os the output stream
     * @param r the MoveRange to print
     */
    friend std::ostream& operator<<(std::ostream& os, const MoveRange& r);
};

/**
 * Operator overloading. Outputs the MoveRange to the output stream
 * @param output The output stream.
 * @param r The MoveRange to print.
 */
std::ostream& operator<<(std::ostream& output, const MoveRange& r);

// ============================================================================
// SARANGE DEFINITION
// ============================================================================

typedef MoveRange SARange; // the range in the suffix array in the move table

// ============================================================================
// CLASS TextOccurrence
// ============================================================================

/**
 * An occurrence in the text. An occurrence is a range in the text and a
 * distance to the mapped read. The distance can be either an edit distance or
 * a hamming distance. The occurrence can also contain an
 * assigned sequence, a line in SAM format and a flag representing whether this
 * occurrence is along the forward or reverse complemented strand.
 */
class TextOcc {
  private:
    Range range;       // the range in the text
    length_t distance; // the distance to this range (edit or hamming)

    length_t assignedSequenceID = -1; // the ID of the assigned sequence

    bool seqNameChecked = false; // indicates whether the sequence name was
    // found for this occurrence

    Strand strand = FORWARD_STRAND; // the strand on which this occurrence lies

    std::string outputLine = ""; // the output line for this text Occ

    length_t indexBegin; // the position in the indexed text where the
                         // occurrence starts

    /**
     * Get the flags for the SAM line of this occurrence (single-ended).
     * @param primaryAlignment Indicates whether this is the primary alignment
     * @returns the flags for the SAM line of this occurrence.
     */
    uint16_t getFlagsSE(bool primaryAlignment) const {
        uint16_t result = 0;

        // set rev complement flag
        result |= (isRevCompl() << 4);

        // set secondary alignment flag (256)
        result |= (!primaryAlignment << 8);

        return result;
    }

  public:
    /**
     * Constructor
     * @param range The range of this occurrence in the text.
     * @param distance The (edit or hamming) distance to the mapped read of
     * this occurrence.
     * @param strand The strand on which this occurrence lies.
     */
    TextOcc(Range range, length_t distance, Strand strand)
        : range(range), distance(distance), strand(strand),
          indexBegin(range.getBegin()) {
    }

    /**
     * Constructor for an invalid text occurrence
     */
    TextOcc() : range(0, 0) {
    }

    /**
     * @returns the range of this occurrence
     */
    Range& getRange() {
        return range;
    }

    /**
     * @returns the range of this occurrence
     */
    const Range& getRange() const {
        return range;
    }

    /**
     * @returns the distance of this occurrence
     */
    const length_t getDistance() const {
        return distance;
    }

    /**
     * Sets the distance of this occurrence
     * @param distance the distance to set
     */
    void setDistance(length_t distance) {
        this->distance = distance;
    }

    /**
     * @returns the output line of this occurrence
     */
    const std::string& getOutputLine() const {
        return outputLine;
    }

    void generateTaggingOutputLine(const std::string& seqID,
                                   const std::string& tag) {
        outputLine = fmt::format("{}\t{}\n", seqID, tag);
    }

    const length_t getAssignedSequenceID() const {
        return assignedSequenceID;
    }

    /**
     * Operator overloading for sorting the occurrences.
     * Occurrences are first sorted on their begin position, then on their
     * distance and their length
     * @param r the occurrence to compare to this (=right hand side of the
     * comparison)
     */
    bool operator<(const TextOcc& r) const {

        if (range.getBegin() != r.getRange().getBegin()) {
            return range.getBegin() < r.getRange().getBegin();
        } else {
            // begin is equal, better ed is smarter
            if (distance != r.getDistance()) {
                return distance < r.getDistance();
            } else if (range.width() != r.getRange().width()) {
                // shorter read is smaller...
                return range.width() < r.getRange().width();
            } else {
                return false;
            }
        }
    }

    /**
     * Operator overloading.  An occurrence is smaller than a value if its
     * distance is smaller than the value.
     */
    bool operator<(length_t value) const {
        return getDistance() < value;
    }

    /**
     * Operator overloading.  Two TextOcc are equal if their ranges and distance
     * are equal.
     * @param r the occurrence to compare to this (=right hand side of the
     * comparison)
     */
    bool operator==(const TextOcc& r) {
        return r.getRange() == range && r.getDistance() == distance;
    }

    /**
     * @returns the width of this occurrence.
     */
    length_t width() const {
        return range.width();
    }

    /**
     * @returns true if this is a valid occurrence
     */
    bool isValid() const {
        return !range.empty();
    }

    /**
     * @returns true if the occurrence is along the reverse complemented strand,
     * false otherwise.
     */
    bool isRevCompl() const {
        return strand == REVERSE_C_STRAND;
    }

    const length_t getEnd() const {
        return range.getEnd();
    }
    const length_t getBegin() const {
        return range.getBegin();
    }

    const Strand getStrand() const {
        return strand;
    }

    void setIndexBegin(length_t indexBegin) {
        this->indexBegin = indexBegin;
    }

    const length_t getIndexBegin() const {
        return indexBegin;
    }

    const length_t getIndexEnd() const {
        return indexBegin + range.width();
    }
};

// ============================================================================
// CLASS SA RANGE PAIR
// ============================================================================

class ToeholdInterface {
  private:
    length_t toehold; // the toehold, which is one occurrence of the current
                      // match in the text
    bool toeholdRepresentsEnd; // indicates whether the toehold represents the
                               // end of the match
    length_t originalDepth;    // the depth of the original match, needed for
                            // offsetting the toehold at the end since the FMPos
                            // depth is decremented with shifts, and there is no
                            // way of knowing how muc

  public:
    /**
     * Constructor. Creates a toehold interface.
     * @param toehold the toehold
     * @param toeholdRepresentsEnd indicates whether the toehold represents the
     * @param originalDepth the depth of the original match (??)
     */
    ToeholdInterface(length_t toehold, bool toeholdRepresentsEnd,
                     length_t originalDepth)
        : toehold(toehold), toeholdRepresentsEnd(toeholdRepresentsEnd),
          originalDepth(originalDepth) {
    }
    /**
     * @returns the toehold
     */
    length_t getToehold() const {
        return toehold;
    }

    /**
     * @brief See if the toehold represents the end of the match
     *
     * @return true if the toehold represents the end of the match
     * @return false otherwise
     */
    bool getToeholdRepresentsEnd() const {
        return toeholdRepresentsEnd;
    }

    /**
     * @returns the original depth of the match
     */
    length_t getOriginalDepth() const {
        return originalDepth;
    }

    bool operator==(const ToeholdInterface& o) const {
        return o.getToehold() == toehold &&
               o.getToeholdRepresentsEnd() == toeholdRepresentsEnd &&
               o.getOriginalDepth() == originalDepth;
    }
};

/**
 * A pair of ranges. The first range is range over the suffix array of the
 * text. The second range is the corresponding range over the suffix array
 * of the reversed text
 */
class SARangePair
    : public ToeholdInterface // include toehold info in b-move
{
  private:
    SARange rangeSA;    // the range over the suffix array
    SARange rangeSARev; // the range over the suffix array of the reversed text

  public:
    /**
     * Default constructor, creates two empty ranges
     */
    SARangePair()
        :
          ToeholdInterface(0, false, 0),
          rangeSA(SARange()), rangeSARev(SARange()) {
    }

    /**
     * Constructor. Creates a pair of ranges.
     * @param rangeSA the range over the suffix array
     * @param rangeSARev the range over the suffix array of the reversed text
     * @param toehold the toehold
     * @param toeholdRepresentsEnd indicates whether the toehold represents the
     */
    SARangePair(SARange rangeSA, SARange rangeSARev, length_t toehold,
                bool toeholdRepresentsEnd, length_t originalDepth)
        : ToeholdInterface(toehold, toeholdRepresentsEnd, originalDepth),
          rangeSA(rangeSA), rangeSARev(rangeSARev) {
    }

    /**
     * @returns the range over the suffix array
     */
    const SARange& getRangeSA() const {
        return rangeSA;
    }

    /**
     * @returns the range over the suffix array of the reversed text
     */
    const SARange& getRangeSARev() const {
        return rangeSARev;
    }

    /**
     * Provides mutable access to the SA range.
     * Note: Modifying the returned reference directly can lead to unintended
     * side effects. Use with caution.
     *
     * @returns a mutable reference to the SA Range.
     */
    SARange& getRangeSAMutable() {
        return rangeSA;
    }

    /**
     * Provides mutable access to the reverse SA range.
     * Note: Modifying the returned reference directly can lead to unintended
     * side effects. Use with caution.
     *
     * @returns a mutable reference to the reverse SA Range.
     */
    SARange& getRangeSARevMutable() {
        return rangeSARev;
    }

    /**
     * @returns true if the range over the suffix array is empty, false
     * otherwise
     */
    bool empty() const {
        return rangeSA.empty();
    }

    /**
     * @returns the width of the ranges, calculated by using the range over the
     * suffix array
     */
    length_t width() const {
        return rangeSA.width();
    }
    /**
     * Operator overloading
     * @param o the other pair (=right hand side of equation)
     * @returns true if this is equal to rhs
     */
    bool operator==(const SARangePair& o) const {
        // only the first range matters as the ranges imply each other
        return o.getRangeSA() == rangeSA
               && ToeholdInterface::operator==(o)
            ;
    }

    /**
     * Function for debug purposes. With bidirectional search we expect this to
     * be true. With unidirectional backwards search this should be false.
     * @returns true if the ranges are synchronized, false otherwise
     */
    bool isSynchronized() const {
        return rangeSA.width() == rangeSARev.width();
    }
};

class SARangeWithToehold : public SARange, public ToeholdInterface {
  public:
    SARangeWithToehold() : SARange(), ToeholdInterface(0, false, 0) {
    }

    SARangeWithToehold(SARange range, length_t toehold,
                       bool toeholdRepresentsEnd, length_t originalDepth)
        : SARange(range),
          ToeholdInterface(toehold, toeholdRepresentsEnd, originalDepth) {
    }

    const SARange& getRange() const {
        return *this;
    }

    SARange getRange() {
        return *this;
    }
};

// ============================================================================
// SARANGEBACKWARDS DEFINITION
// ============================================================================

typedef SARangeWithToehold SARangeBackwards; // the range in the suffix array in
                                             // the move table, with toehold

// ============================================================================
// CLASS FMPos
// ============================================================================
/**
 * A position in the bidirectional FM-index.
 */
class FMPos {
  protected:
    SARangePair ranges; // the ranges over the suffix arrays
    length_t depth; // the depth of the prefix of the suffixes of this position

  public:
    /**
     * Default constructor for empty position (= empty ranges and depth of
     * zero)
     */
    FMPos() : ranges(SARangePair()), depth(0) {
    }

    FMPos(SARangePair& ranges, length_t depth) : ranges(ranges), depth(depth) {
    }

    const SARangePair& getRanges() const {
        return ranges;
    }

    /**
     * Provides mutable access to the ranges.
     * Note: Modifying the returned reference directly can lead to unintended
     * side effects. Use with caution.
     *
     * @returns a mutable reference to the ranges.
     */
    SARangePair& getRangesMutable() {
        return ranges;
    }

    const length_t& getDepth() const {
        return depth;
    }

    void setRanges(SARangePair ranges) {
        this->ranges = ranges;
    }

    void setDepth(length_t depth) {
        this->depth = depth;
    }

    /**
     * Operator overloading, two FMPos are equal if their ranges and depth
     * are equal
     * @param rhs the FMPos to compare to this
     * @returns true if this is equal to rhs
     */
    bool operator==(const FMPos& rhs) const {
        return ranges == rhs.getRanges() && depth == rhs.getDepth();
    }
    /**
     * @returns true if the ranges are not empty, false otherwise
     */
    bool isValid() const {
        return !ranges.empty();
    }
};
// ============================================================================
// CLASS FMOcc
// ============================================================================

/**
 * An occurrence in the bidirectional FM-index
 */
class FMOcc {
  private:
    FMPos pos;         // The FM position of this occurrence
    length_t distance; // the distance (hamming or edit)
    length_t shift; // A right-sift to the corresponding positions in the text

  public:
    FMOcc() : pos(), distance(0), shift(0) {
    }
    /**
     * Make a bidirectional approximate match in the suffix array
     * @param ranges the ranges of this approximate match (range in SA and
     * in SA')
     * @param distance the (edit or hamming) distance of this approximate
     * match
     * @param depth the depth (=length) of this approximate match
     * @param shift The right shift to the corresponding positions in the
     * text, defaults to zero
     */
    FMOcc(SARangePair ranges, length_t distance, length_t depth,
          length_t shift = 0)
        : pos(ranges, depth), distance(distance), shift(shift) {
    }
    /**
     * Make a bidirectional approximate match in the suffix array
     * @param pos the position in the FMIndex of this approximate match
     * @param distance the (edit or hamming) distance of this approximate
     * match
     * @param shift The right shift to the corresponding positions in the
     * text, defaults to zero
     */
    FMOcc(FMPos pos, length_t distance, length_t shift = 0)
        : pos(pos), distance(distance), shift(shift) {
    }

    const SARangePair& getRanges() const {
        return pos.getRanges();
    }
    const length_t& getDistance() const {
        return distance;
    }

    const length_t& getDepth() const {
        return pos.getDepth();
    }

    length_t getWidth() const {
        return pos.getRanges().width();
    }

    const length_t& getShift() const {
        return shift;
    }

    void setRanges(SARangePair ranges) {
        pos.setRanges(ranges);
    }

    void setDistance(length_t distance) {
        this->distance = distance;
    }

    void setDepth(length_t depth) {
        pos.setDepth(depth);
    }

    /**
     * @returns true if the position is valid, false otherwise
     */
    bool isValid() const {
        return pos.isValid();
    }

    /**
     * Operator overloading to sort FMOcc
     * First the FMOcc are sorted on the begin of the range over the suffix
     * array of their position Then they are sorted on their distance
     * Lastly they are sorted on the width of their ranges
     * @param rhs the FMOcc to compare to this
     * @returns true if this is smaller than rhs
     */
    bool operator<(const FMOcc& rhs) const {
        if (pos.getRanges().getRangeSA().getBegin() !=
            rhs.getRanges().getRangeSA().getBegin()) {
            // begins of ranges are unequal, return smallest begin
            return getRanges().getRangeSA().getBegin() <
                   rhs.getRanges().getRangeSA().getBegin();
        } else if (distance != rhs.getDistance()) {
            // begin is equal, better ed is smarter
            return distance < rhs.getDistance();
        } else if (getRanges().width() != rhs.getRanges().width()) {
            // shorter read is smaller...
            return getRanges().width() < rhs.getRanges().width();
        } else {
            // prefer no shift
            return getShift() < rhs.getShift();
        }
    }
    /**
     * Operator overloading
     * Two FMOcc are equal if their ranges, distance, depth and shift are all
     * equal
     * @returns true if this is equal to rhs
     */
    bool operator==(const FMOcc& rhs) const {
        return getRanges() == rhs.getRanges() &&
               distance == rhs.getDistance() && getDepth() == rhs.getDepth() &&
               getShift() == rhs.getShift();
    }
    friend std::ostream& operator<<(std::ostream& os, const FMOcc& fmOcc);
};

// ============================================================================
// CLASS FMPosExt
// ============================================================================
/**
 * A single node in the bidirectional FM-index. Its depth is the depth from
 * the start match for a particular phase of a search
 */
class FMPosExt : public FMPos {
  private:
    char c;                // the character of this node
    bool reported = false; // has this particular node already reported?
  public:
    /**
     * Create a node of the search tree
     * @param character the character of this node
     * @param ranges the ranges over the suffix and reversed suffix array
     * that go to this node
     * @param row the row of this node in the alignment matrix = depth of
     * this node
     */
    FMPosExt(char character, SARangePair ranges, length_t row)
        : FMPos(ranges, row), c(character), reported(false) {
    }

    /**
     * Default constructor, this Node will have empty ranges
     */
    FMPosExt() : FMPos(), c(char(0)) {
    }

    /**
     * Sets the report flag to true
     */
    void report() {
        reported = true;
    }

    /**
     * Reports the match (with added depth) at this node,
     * @param occ the match will be stored here
     * @param startDepth the depth to add to the match
     * @param EDFound the found edit distance for this node
     * @param noDoubleReports false if this node is allowed to report more
     * than once, defaults to false
     * @param shift right shift of the match, defaults to zero
     */
    void report(FMOcc& occ, const length_t& startDepth, const length_t& EDFound,
                const bool& noDoubleReports = false, length_t shift = 0) {
        if (!reported) {
            occ = FMOcc(getRanges(), EDFound, depth + startDepth, shift);

            // if finalPiece, report only once
            if (noDoubleReports) {
                report();
            }
        }
    }

    /**
     * Gets the ranges of this node
     * @returns the ranges of this node
     */
    const SARangePair& getRanges() const {
        return ranges;
    }

    /**
     * Get the character of this node
     * @returns the character of this node
     */
    const char getCharacter() const {
        return c;
    }

    /**
     * Get the row of this node
     * @returns the row of this node
     */
    length_t getRow() const {
        return depth;
    }
};

// ============================================================================
// CLASS COUNTERS
// ============================================================================
/**
 * @brief Class to manage and manipulate various performance counters.
 */
class Counters {
  public:
    enum CounterType {
        NODE_COUNTER, // counts the number of nodes visited in the index
        TOTAL_REPORTED_POSITIONS, // counts the number of matches reported
                                  // (either via in-text verification or
                                  // in-index matching)
        SEARCH_STARTED,   // counts the number of times a search started
        DROPPED_UNIQUE_MATCHES, // counts the number of unique matches dropped
                                // because of one-line reporting

        // Aggregational counters
        NUMBER_OF_READS,           // counts the number of reads
        TOTAL_UNIQUE_MATCHES,      // counts the number of all unique matches
        MAPPED_READS,              // counts the number of mapped reads

        TOTAL_STEP_COUNTS,                  // counts the total number of steps in the tagging process
        TOTAL_READ_LENGTHS,                 // counts the total lengths of reads processed
#ifdef GROUND_TRUTH_CHECKS  // Tagging-related counters
        TOTAL_CORRECTLY_IDENTIFIED_READS,   // counts the number of correctly identified reads (first set)
        TOTAL_UNIDENTIFIED_READS,           // counts the number of unidentified reads (first set)
        TOTAL_INCORRECTLY_IDENTIFIED_READS, // counts the number of incorrectly identified reads (first set)
        TOTAL_CORRECTLY_IDENTIFIED_READS2,  // counts the number of correctly identified reads (second set)
        TOTAL_UNIDENTIFIED_READS2,          // counts the number of unidentified reads (second set)
        TOTAL_INCORRECTLY_IDENTIFIED_READS2,// counts the number of incorrectly identified reads (second set)
        TOTAL_CORRECTLY_IDENTIFIED_READS3,  // counts the number of correctly identified reads (third set)
        TOTAL_UNIDENTIFIED_READS3,          // counts the number of unidentified reads (third set)
        TOTAL_INCORRECTLY_IDENTIFIED_READS3,// counts the number of incorrectly identified reads (third set)
#endif

        COUNTER_TYPE_MAX // To denote the number of counters
    };

  private:
    std::array<uint64_t, COUNTER_TYPE_MAX> counters; // array with the counters

  public:

    /**
     * @brief Constructor to initialize counters. Resets all counters to 0.
     */
    Counters() {
        resetCounters();
    };

    // Reset all counters and clear unordered maps
    void resetCounters() {
        counters.fill(0);
    }

    /**
     * @brief Increments the specified counter by a given amount.
     *
     * @param type The type of counter to increment.
     * @param amount [optional] The amount by which to increment the counter
     * (default is 1).
     */
    void inc(CounterType type, uint64_t amount = 1) {
        counters[type] += amount;
    }

    /**
     * @brief Retrieves the current value of the specified counter.
     *
     * @param type The type of counter to retrieve.
     * @return The current value of the counter.
     */
    uint64_t get(CounterType type) const {
        return counters[type];
    }

    /**
     * @brief Adds counters from another Counters object to this one.
     *
     * @param o The Counters object from which to add counters.
     */
    void addCounters(const Counters& o) {
        for (int i = 0; i < COUNTER_TYPE_MAX; ++i) {
            counters[i] += o.counters[i];
        }
    }

    /**
     * @brief Reports statistics based on the counters and a given sequencing
     * mode.
     */
    void reportStatistics() const;
};

// ============================================================================
// CLASS Occurrences
// ============================================================================
// This class combines the occurrences in the text and the occurrences
// in the fm index into one data structure
/**
 * @class Occurrences
 * @brief Represents a collection of occurrences in the FM index and the text.
 *
 * The Occurrences class stores the in-index occurrences (FMOcc) and the in-text
 * occurrences (TextOcc). It provides methods to add occurrences, erase
 * duplicates, and retrieve the occurrences.
 */
class Occurrences {
  private:
    std::vector<TextOcc> inTextOcc; // the in-text occurrences
    std::vector<FMOcc> inFMOcc;     // the in-index occurrences

  public:
    /**
     * @brief Constructs an Occurrences object with optional reserve capacity.
     *
     * @param reserve The initial capacity to reserve for the occurrences
     * vectors.
     */
    Occurrences(const length_t reserve = 200) {
        inTextOcc.reserve(reserve);
        inFMOcc.reserve(reserve);
    }

    /**
     * @brief Adds an in-index occurrence to the collection.
     *
     * @param match The in-index occurrence to add.
     */
    void addFMOcc(const FMOcc& match) {
        inFMOcc.emplace_back(match);
    }

    /**
     * @brief Adds an in-text occurrence to the collection.
     *
     * @param occ The in-text occurrence to add.
     */
    void addTextOcc(const TextOcc& occ) {
        inTextOcc.emplace_back(occ);
    }

    /**
     * @brief Adds the given occurrences to the inTextOcc vector.
     *
     * @param occs The occurrences to add.
     */
    void addTextOccs(const std::vector<TextOcc>& occs) {
        inTextOcc.insert(inTextOcc.end(), occs.begin(), occs.end());
    }

    /**
     * @brief Sets the inTextOcc vector to the given occurrences.
     *
     * @param occs The occurrences to set.
     */
    void setTextOccs(std::vector<TextOcc>& occs) {
        inTextOcc = occs;
    }

    /**
     * @brief Sets the inTextOcc vector to the given occurrences.
     * The given occurrences are moved into the inTextOcc vector.
     * @param occs The occurrences to set.
     */
    void moveTextOccs(std::vector<TextOcc>& occs) {
        inTextOcc = std::move(occs);
    }

    /**
     * @brief Erases all duplicate in-index occurrences and sorts the
     * occurrences.
     */
    void eraseDoublesFM() {
#ifdef DEVELOPER_MODE
        stable_sort(inFMOcc.begin(), inFMOcc.end());
#else
        sort(inFMOcc.begin(), inFMOcc.end());
#endif

        inFMOcc.erase(unique(inFMOcc.begin(), inFMOcc.end()), inFMOcc.end());
    }

    /**
     * @brief Erases all duplicate in-text occurrences and sorts the
     * occurrences.
     */
    void eraseDoublesText() {
#ifdef DEVELOPER_MODE
        stable_sort(inTextOcc.begin(), inTextOcc.end());
#else
        sort(inTextOcc.begin(), inTextOcc.end());
#endif
        inTextOcc.erase(unique(inTextOcc.begin(), inTextOcc.end()),
                        inTextOcc.end());
    }

    /**
     * @brief Returns the in-index occurrences.
     *
     * @return The in-index occurrences.
     */
    const std::vector<FMOcc>& getFMOccurrences() const {
        return inFMOcc;
    };

    /**
     * @brief Returns the number of in-text occurrences.
     *
     * @return The number of in-text occurrences.
     */
    size_t textOccSize() const {
        return inTextOcc.size();
    }

    /**
     * @brief Returns the in-text occurrences.
     *
     * @return The in-text occurrences.
     */
    const std::vector<TextOcc>& getTextOccurrences() const {
        return inTextOcc;
    }

    /**
     * @brief Returns the maximum size between in-text and in-index occurrences.
     *
     * @return The maximum size between in-text and in-index occurrences.
     */
    length_t getMaxSize() const {
        return std::max(inTextOcc.size(), inFMOcc.size());
    }

    /**
     * Checks if the Occurrences object is empty.
     *
     * @return true if the Occurrences object is empty, false otherwise.
     */
    bool empty() const {
        return inTextOcc.empty() && inFMOcc.empty();
    }
};

#endif // end of include guard: INDEXHELPERS_H
