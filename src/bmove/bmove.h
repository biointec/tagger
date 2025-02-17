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

#ifndef BMOVE_H
#define BMOVE_H

#include "../definitions.h"
#include "../indexhelpers.h"
#include "../indexinterface.h"
#include "../reads.h"
#include "moverepr.h"
#include "sparsebitvec.h"

#include <ios>
#include <sdsl/int_vector.hpp>
#include <string>
#include <vector>
class MemoryMappedTextFile;
class Substring;

class BMove : public IndexInterface {
  private:

    MoveLFReprBP move;
    MoveLFReprBP moveR;

    // Tagging data
    sdsl::int_vector<> samplesTagging;
#ifdef TAG_ARRAY_SUBSAMPLING
    SparseBitvec samplesTaggingBV;
    length_t maxLFSteps;
#endif
    length_t initialToeholdTagging;

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
                           ) override;

    /**
     * Read a binary file and stores content in sdsl int_vector
     * @param filename File name
     * @param intVector output int_vector (contents will be overwritten)
     * @returns True if successful, false otherwise
     */
    bool readIntVector(const std::string& filename,
                       sdsl::int_vector<>& intVector) {
        std::ifstream ifs(filename, std::ios::binary);
        if (!ifs) {
            return false;
        }
        intVector.load(ifs);
        ifs.close();
        return true;
    }

    // ----------------------------------------------------------------------------
    // ROUTINES FOR ACCESSING DATA STRUCTURE
    // ----------------------------------------------------------------------------

    /**
     * @brief Get the tagging sample using subsampling
     *
     * @param pos The position in the text to get the tagging sample for
     * @param run The run to get the tagging sample for
     * @return length_t - The tagging sample
     */
    length_t getTaggingSample(length_t pos, length_t run) const;

    /**
     * Get toehold for the given character to match.
     * @param range Range in the sa array in which to match char and get
     * toehold.
     * @param c Character to match.
     * @return Index in the SA where the toehold starts.
     */
    length_t computeToeholdTagging(const SARange& range, const length_t c,
                                   const length_t previousToehold) const;

    // ----------------------------------------------------------------------------
    //  APPROXIMATE PATTERN MATCHING
    // ----------------------------------------------------------------------------

    /**
     * @brief Update the runs corresponding to the SA range in the reverse text
     *
     * @param ranges SARangePair for which the SA range runs must be updated
     */
    virtual void updateRangeSARevRuns(SARangePair& ranges) const override;

    /**
     * Helper function for backwards adding of character (both uni- and
     * bidirectional). Finds the backwards (trivial) range for extending the
     * parent range with a character.
     * @param positionInAlphabet the position in the alphabet of the character.
     * @param parentBackwardRange the backward range of the parent
     * @param childBackwardRange the backward range of the child (output)
     *
     */
    void
    findRangeWithExtraCharBackwardAuxiliary(length_t positionInAlphabet,
                                            SARange& parentBackwardRange,
                                            SARange& childBackwardRange) const;

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
                                    SARangePair& childRanges) const override;

    /**
     * Finds the ranges of string Pc using the principle explained in the
     * paper of Lam
     * @param positionInAlphabet the position in alphabet of the character c
     * that is added in the back
     * @param rangesOfP the ranges of pattern P
     * @param childRanges the ranges corresponding to string Pc, this will
     * be set during execution
     * @returns true if the range is not empty, false otherwise
     */
    virtual bool
    findRangesWithExtraCharForward(length_t positionInAlphabet,
                                   const SARangePair& rangesOfP,
                                   SARangePair& childRanges) const override;

    /**
     * Finds the range of Pc using unidirectional backward matching
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
                                   SARangeBackwards& childRange) const override;

    // ----------------------------------------------------------------------------
    // TAGGING ROUTINES
    // ----------------------------------------------------------------------------

    int64_t findMEMStartWithLookup(const std::string& read, int64_t start,
                                   SARange& range, length_t& stepCount,
                                   bool& encounteredN, int64_t readSize,
                                   length_t& toeHold) const;

    int64_t findMEMEnd(const std::string& read, int64_t i, int64_t j,
                       length_t& stepCount) const;

    int64_t findMEMEndWithLookup(const std::string& read, int64_t i, int64_t j,
                                 length_t& stepCount) const;

    uint16_t computeTagDetails(const Range& range, length_t& numberOfTagsBefore,
                               length_t& numberOfTagsAfter,
                               length_t& toeHold) const;

    void
    updateDetectedTags(Counters& counters, int length, uint16_t detectedTag,
                       const uint16_t correctTag,
                       std::unordered_map<uint16_t, length_t>& detectedTags,
                       length_t& totalWeight) const;

    uint16_t
    computeTagMEM(Counters& counters, int64_t i, int64_t j, const Range& range,
                  const uint16_t correctTag, bool& allMEMsunique,
                  std::unordered_map<uint16_t, length_t>& detectedTags,
                  std::unordered_map<uint16_t, int64_t>& lastSeenPositions,
                  length_t& totalWeight, length_t& toeHold) const;

    void
    processMatches(const std::string& read, int64_t i, int64_t j,
                   Counters& counters, std::vector<TextOcc>& occs,
                   const SARange& range, const Strand strand,
                   const uint16_t correctTag, const std::string& readID,
                   bool& allMEMsunique,
                   std::unordered_map<uint16_t, length_t>& detectedTags,
                   std::unordered_map<uint16_t, int64_t>& lastSeenPositions,
                   length_t& totalWeight, length_t& toeHold) const;

    int64_t skipNs(const std::string& read, int64_t pos,
                   int64_t readSize) const;

  public:
    // ----------------------------------------------------------------------------
    // INITIALIZATION ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * Constructor
     * @param baseFile The base name of the files that contain the info
     * about the index.
     * @param verbose If true, the steps will be written to cout. [default =
     * true]
     * @param wordSize The size of the mers to be stored in the hashtable.
     * Used for quick look-ups of exact seeds. [default = 10]
     */
    BMove(const std::string& baseFile, bool verbose = true,
          length_t wordSize = 10
#ifdef TAG_ARRAY_SUBSAMPLING
          ,
          length_t maxLF = DEFAULT_MAX_LF
#endif // TAG_ARRAY_SUBSAMPLING
          )
        : IndexInterface(baseFile, verbose, wordSize) {

        // read in files
        fromFiles(baseFile, verbose
#ifdef TAG_ARRAY_SUBSAMPLING
                  ,
                  maxLF
#endif // TAG_ARRAY_SUBSAMPLING
        );

        // populate sparse hash table
        populateTable(verbose);
    }

    /**
     * @brief Destructor
     *
     */
    virtual ~BMove() override = default;

    // ----------------------------------------------------------------------------
    // ROUTINES FOR ACCESSING THE DATA STRUCTURE
    // ----------------------------------------------------------------------------

    /**
     * Get the complete range of this index corresponding to an exact match
     * of the empty string.
     * @returns A SARangePair with both ranges the complete range of the
     * index.
     */
    virtual SARangePair getCompleteRange() const override {
        return SARangePair(SARange(0, textLength, 0, move.size() - 1),
                           SARange(0, textLength, 0, moveR.size() - 1),
                           initialToeholdTagging, false, 0);
    }

    SARange getCompleteRangeBackward() const {
        return SARange(0, textLength, 0, move.size() - 1);
    }

    SARange getCompleteRangeForward() const {
        return SARange(0, textLength, 0, moveR.size() - 1);
    }

    // ----------------------------------------------------------------------------
    // ROUTINES FOR APPROXIMATE MATCHING
    // ----------------------------------------------------------------------------

    virtual bool
    taggingFunction(const std::string& readID, const std::string& read,
                    std::vector<TextOcc>& occs, Counters& counters,
                    const Strand strand, const uint16_t correctTag,
                    const length_t minMEMLength, bool& allMEMsunique,
                    std::unordered_map<uint16_t, length_t>& detectedTags,
                    std::unordered_map<uint16_t, int64_t>& lastSeenPositions,
                    length_t& totalWeight) const override;

#ifdef GROUND_TRUTH_CHECKS
    virtual uint16_t
    getCorrectTag(const std::string& refID) const override {
        // Find the last underscore
        size_t underscorePos = refID.find_last_of('_');
        assert(underscorePos != std::string::npos &&
               "Invalid refID: no underscore found");
        assert(underscorePos < refID.size() - 1 &&
               "Invalid refID: no numeric suffix");

        // Extract numeric part after last underscore
        const char* numericPart = refID.c_str() + underscorePos + 1;

        // Validate numeric part (ensure it's all digits)
        for (const char* p = numericPart; *p; ++p) {
            assert(std::isdigit(*p) &&
                   "Invalid refID: non-digit in numeric suffix");
        }

        // Convert to integer safely
        unsigned long value = std::strtoul(numericPart, nullptr, 10);
        assert(value <= UINT16_MAX && "Numeric value exceeds uint16_t range");

        return static_cast<uint16_t>(value);
    }
#endif // GROUND_TRUTH_CHECKS
};

#endif // BMOVE_H