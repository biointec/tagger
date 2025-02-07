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

#ifndef INDEXINTERFACE_H
#define INDEXINTERFACE_H

#include "alphabet.h"
#include "definitions.h"       // for length_t, PairStatus, Strand
#include "indexhelpers.h"      // for SARangePair, Range, TextOcc
#include "reads.h"             // for ReadBundle
#include "substring.h"         // for Substring
#include "tkmer.h"             // for Kmer, KmerHash

#include <algorithm>                // for max
#include <parallel_hashmap/phmap.h> // phmap::parallel_flat_hash_map
#include <stddef.h>                 // for size_t
#include <stdint.h>                 // for uint16_t
#include <string>                   // for string, char_traits<>::pos_type
#include <unordered_map>            // for unordered_map
#include <utility>                  // for pair
#include <vector>                   // for vector

class MemoryMappedTextFile;
class Search;

class IndexInterface;
typedef bool (IndexInterface::*ExtraCharPtr)(length_t, const SARangePair&,
                                             SARangePair&) const;

/**
 * Class representing the bidirectional index and all operations on that
 * Index
 */
class IndexInterface {
  protected:
    // info about the text
    const std::string baseFile; // The base file of the reference text
    length_t textLength;        // the length of the text

    std::vector<length_t> counts; // the counts array of the reference genome

    Alphabet<ALPHABET> sigma; // the alphabet

    // direction variables
    thread_local static Direction dir; // the direction of the index
    thread_local static ExtraCharPtr
        extraChar; // pointer to extra char method (for direction)

    // sparse hash info
    const size_t wordSize = 10; // the size of the mers to be stored in a table
    phmap::parallel_flat_hash_map<Kmer, SARangePair, KmerHash>
        table; // hashtable that contains all wordSize-mers

    // ----------------------------------------------------------------------------
    // PREPROCESSING ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * Private helper function that reads in all the necessary files
     * @param baseFile the baseFile of the files that will be read in
     * @param verbose if true the steps will we written to cout
     */
    virtual void fromFiles(const std::string& baseFile, bool verbose
#ifdef TAG_ARRAY_SUBSAMPLING
                           ,
                           length_t maxLF = DEFAULT_MAX_LF
#endif // TAG_ARRAY_SUBSAMPLING
                           ) = 0;

    /**
     * Read the counts array, the build_tag info, the flavor and the
     * compilation-bits info of the index.
     * @param baseFile The base name of the files that contain the info about
     * the index.
     * @param verbose If true, the steps will be written to cout.
     */
    void readMetaAndCounts(const std::string& baseFile, bool verbose);

    /**
     * Read a binary file and stores content in array
     * @param filename File name
     * @param array Suffix array (contents will be overwritten)
     * @returns True if successful, false otherwise
     */
    static bool readArray(const std::string& filename,
                          std::vector<length_t>& array);

    /**
     * Populate the hash table
     * @param verbose if steps are written to cout
     */
    void populateTable(bool verbose);

    /**
     * @brief Update the runs corresponding to the SA range in the reverse
     *
     * @param ranges  SARangePair for which the SA range runs must be
     * updated
     */
    virtual void updateRangeSARevRuns(SARangePair& ranges) const = 0;

    /**
     * Finds the ranges of cP using the principle explained in the paper of
     * Lam.
     * @param positionInAlphabet the position in alphabet of the character
     * that is added in the front.
     * @param rangesOfP the ranges of pattern P.
     * @param childRanges the ranges corresponding to string cP, this
     * will be set during execution
     * @returns true if the range is not empty, false otherwise
     */
    virtual bool
    findRangesWithExtraCharBackward(length_t positionInAlphabet,
                                    const SARangePair& rangesOfP,
                                    SARangePair& childRanges) const = 0;

    /**
     * Finds the ranges of string Pc using the principle explained in the paper
     * of Lam
     * @param positionInAlphabet the position in alphabet of the character c
     * that is added in the back
     * @param rangesOfP the ranges of pattern P
     * @param childRanges the ranges corresponding to string Pc, this will be
     * set during execution
     * @returns true if the range is not empty, false otherwise
     */
    virtual bool
    findRangesWithExtraCharForward(length_t positionInAlphabet,
                                   const SARangePair& rangesOfP,
                                   SARangePair& childRanges) const = 0;

    /**
     * Finds the range of Pc using unidirectional backward matching. With
     * toeholds in the case of run-length compression
     * @param positionInAlphabet the position in alphabet of the character c
     * that is added in the back
     * @param rangeOfP the range over the suffix array of the text
     * corresponding to pattern P
     * @param childRange the range over the suffix array of text
     * corresponding to pattern Pc this will be set during execution
     * @returns true if the range is not empty, false otherwise
     */
    virtual bool
    findRangeWithExtraCharBackward(length_t positionInAlphabet,
                                   const SARangeBackwards& rangeOfP,
                                   SARangeBackwards& childRange

    ) const = 0;

    // ----------------------------------------------------------------------------
    // EXTEND ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * Creates the children nodes of the node with the given ranges and pushes
     * them to the stack.
     * @param ranges The ranges of the parent node to get the children of.
     * @param stack The stack to push the children on.
     * @param counters The performance counters.
     * @param row The row in the current matrix of the parentNode (defaults to
     * 0).
     */
    void extendFMPos(const SARangePair& ranges, std::vector<FMPosExt>& stack,
                     Counters& counters, length_t row = 0) const;

    // ----------------------------------------------------------------------------
    // SUBROUTINES FOR PATTERN MATCHING
    // ----------------------------------------------------------------------------

  public:
    // ----------------------------------------------------------------------------
    // INITIALIZATION ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * Constructor
     * @param baseFile The base name of the files that contain the info about
     * the index.
     * @param verbose If true, the steps will be written to cout. [default =
     * true]
     * @param wordSize The size of the mers to be stored in the hashtable. Used
     * for quick look-ups of exact seeds. [default = 10]
     */
    IndexInterface(const std::string& baseFile, bool verbose = true,
                   length_t wordSize = 10)
        : baseFile(baseFile), wordSize(wordSize) {
    }

    /**
     * @brief Destructor
     *
     */
    virtual ~IndexInterface() = default;

    /**
     * Get the complete range of this index corresponding to an exact match of
     * the empty string.
     * @returns A SARangePair with both ranges the complete range of the
     * index.
     */
    virtual SARangePair getCompleteRange() const = 0;

    // ----------------------------------------------------------------------------
    // ROUTINES FOR ACCESSING THE DATA STRUCTURE
    // ----------------------------------------------------------------------------

    /**
     * Get the length of the text
     */
    const length_t getTextSize() const {
        return textLength;
    }

    /**
     * @returns the wordsize of the mers stored in the table
     */
    length_t getWordSize() const {
        return wordSize;
    }

    /**
     * Looks up the SARangePair of the exact match of c p in the hashtable.
     * Warning: This function assumes p is of size wordSize.
     * @param p The substring to find the ranges of.
     * @returns the ranges corresponding to substring p, if no pair can be
     * found it returns empty ranges.
     */
    SARangePair lookUpInKmerTable(const Substring& p) const {

        auto it = table.find(Kmer(p.getText(), p.begin()));
        return (it == table.end() || p.containsN()) ? SARangePair()
                                                    : it->second;
    }
    // ----------------------------------------------------------------------------
    // ROUTINES FOR EXACT MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Calculates the exact matches to the string in the index. Creates
     * corresponding text-occurrences and adds them to the output vector.
     * Make sure to set the index in the correct mode before calling this
     * function to correctly label the found occurrences.
     * @param s The string to match in the reference genome.
     * @param counters The performance counters.
     * @param tOcc The vector to which the found text occurrences will be added.
     */
    bool exactMatchesOutput(const std::string& s, Counters& counters,

                            std::vector<TextOcc>& tOcc) const;

    // ----------------------------------------------------------------------------
    // ROUTINES FOR APPROXIMATE MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Sets the search direction of the fm-index
     * @param d the direction to search in, either FORWARD or BACKWARD
     * @param isUniDirectionalBackward if the search is unidirectional backward
     */
    void setDirection(Direction d) {
        dir = d;
        extraChar =((d == FORWARD)
                       ? &IndexInterface::findRangesWithExtraCharForward
                       : &IndexInterface::findRangesWithExtraCharBackward);
    }

    // ----------------------------------------------------------------------------
    // TAGGING ROUTINES
    // ----------------------------------------------------------------------------

    virtual bool
    taggingFunction(const std::string& readID, const std::string& read,
                    std::vector<TextOcc>& occs, Counters& counters,
                    const Strand strand, const uint16_t correctTag,
                    const length_t minMEMLength, bool& allMEMsunique,
                    std::unordered_map<uint16_t, length_t>& detectedTags,
                    std::unordered_map<uint16_t, int64_t>& lastSeenPositions,
                    length_t& totalWeight) const = 0;

    virtual uint16_t getCorrectTag(const std::string& refID) const = 0;
};

#endif
