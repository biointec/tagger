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

#ifndef SEARCHSTRATEGY_H
#define SEARCHSTRATEGY_H

#include "definitions.h"    // for length_t, MappingMode, Distance...
#include "indexhelpers.h"   // for TextOcc, PairedTextOccs, Occurrences
#include "indexinterface.h" // for the IndexInterface
#include "reads.h"          // for ReadPair, ReadBundle
#include "search.h"         // for Search, SearchScheme
#include "substring.h"      // for Substring

#include <algorithm> // for max, move, min
#include <cassert>   // for assert
#include <cstdint>   // for uint16_t
#include <ios>       // for ifstream, basic_ios
#include <iterator>  // for move_iterator, back_insert_iterator
#include <memory>    // for allocator, allocator_traits<>::value_type
#include <numeric>   // for accumulate
#include <random>    // for mt19937, uniform_int_distribution
#include <set>       // for set
#include <stdexcept> // for runtime_error, invalid_argument
#include <string>    // for string, operator+, char_traits, to_string
#include <utility>   // for pair, move
#include <vector>    // for vector, _Bit_iterator

#include <sys/stat.h> // for stat on POSIX systems
#include <unistd.h>   // For access on POSIX systems
#define PATH_SEPARATOR "/"

// ============================================================================
// CLASS SEARCHSTRATEGY
// ============================================================================

class SearchStrategy {

    // ATTRIBUTES
    //-------------------------------------------------------------------------------
  protected:
    IndexInterface&
        index;        // reference to the index of the text that is searched
    mutable std::mt19937 global_gen; // random number generator

    uint16_t getConsensusTag(
        const std::unordered_map<uint16_t, length_t>& detectedTags) const {
        uint16_t consensusTag = UINT16_MAX;
        length_t max_weighted_sum = 0;
        std::vector<uint16_t> max_candidates;

        for (const auto& pair : detectedTags) {
            if (pair.second > max_weighted_sum) {
                max_weighted_sum = pair.second;
                max_candidates.clear();
                max_candidates.push_back(pair.first);
            } else if (pair.second == max_weighted_sum) {
                max_candidates.push_back(pair.first);
            }
        }

        if (!max_candidates.empty()) {
            std::uniform_int_distribution<size_t> dist(
                0, max_candidates.size() - 1);
            consensusTag = max_candidates[dist(global_gen)];
        }

        return consensusTag;
    }

    uint16_t handleStrandSelection(
        Counters& counters,
        std::unordered_map<uint16_t, length_t>& detectedTagsFwd,
        length_t totalWeightFwd,
        std::unordered_map<uint16_t, length_t>& detectedTagsRevCompl,
        length_t totalWeightRevCompl, uint16_t correctTag,
        bool& allAlignmentsUnique) const {
        // Choose the strand with the higher total weight
        std::unordered_map<uint16_t, length_t>& chosenDetectedTags =
            (totalWeightFwd > totalWeightRevCompl) ? detectedTagsFwd
                                                   : detectedTagsRevCompl;

        // Get the consensus tag using the helper function
        uint16_t consensusTag = getConsensusTag(chosenDetectedTags);

# ifdef GROUND_TRUTH_CHECKS

        // Update counters for the first case
        uint16_t chosenTag = chosenDetectedTags.size() == 1
                                 ? chosenDetectedTags.begin()->first
                                 : UINT16_MAX;
        if (allAlignmentsUnique && chosenDetectedTags.size() == 1 &&
            chosenTag == correctTag) {
            counters.inc(Counters::TOTAL_CORRECTLY_IDENTIFIED_READS);
        } else if (!allAlignmentsUnique || chosenDetectedTags.size() != 1) {
            counters.inc(Counters::TOTAL_UNIDENTIFIED_READS);
        } else {
            counters.inc(Counters::TOTAL_INCORRECTLY_IDENTIFIED_READS);
        }

        allAlignmentsUnique = true; // Reset for second option
        // Update counters for the second case
        chosenTag = chosenDetectedTags.size() == 1
                        ? chosenDetectedTags.begin()->first
                        : UINT16_MAX;
        if (allAlignmentsUnique && chosenDetectedTags.size() == 1 &&
            chosenTag == correctTag) {
            counters.inc(Counters::TOTAL_CORRECTLY_IDENTIFIED_READS2);
        } else if (!allAlignmentsUnique || chosenDetectedTags.size() != 1) {
            counters.inc(Counters::TOTAL_UNIDENTIFIED_READS2);
        } else {
            counters.inc(Counters::TOTAL_INCORRECTLY_IDENTIFIED_READS2);
        }

        // Update counters for the third case - consensus case
        if (consensusTag == correctTag) {
            counters.inc(Counters::TOTAL_CORRECTLY_IDENTIFIED_READS3);
        } else if (consensusTag == UINT16_MAX) {
            counters.inc(Counters::TOTAL_UNIDENTIFIED_READS3);
        } else {
            counters.inc(Counters::TOTAL_INCORRECTLY_IDENTIFIED_READS3);
        }
# endif

        return consensusTag;
    }

  public:

    /**
     * Constructor
     * @param index the index to be used
     */
    SearchStrategy(IndexInterface& index) : index(index), global_gen(42) {
    }

    virtual ~SearchStrategy() {
    }

    std::vector<TextOcc> taggingFunction(const ReadBundle& readBundle,
                                         Counters& counters,
                                         const length_t minMEMLength) const {
        std::vector<TextOcc> occs;
#ifdef GROUND_TRUTH_CHECKS
        uint16_t correctTag = index.getCorrectTag(readBundle.getSeqID()); //todolore
        assert(correctTag < 9118); // todolore hardcoded
#else
        uint16_t correctTag = 0;
#endif
        bool allMEMsunique = true;
        bool mapped = false;

        // Forward strand detection
        std::unordered_map<uint16_t, length_t> detectedTagsFwd;
        std::unordered_map<uint16_t, int64_t> lastSeenPositionsFwd;
        length_t totalWeightFwd = 0;
        mapped |= index.taggingFunction(
            readBundle.getSeqID(), readBundle.getRead(), occs, counters,
            FORWARD_STRAND, correctTag, minMEMLength, allMEMsunique,
            detectedTagsFwd, lastSeenPositionsFwd, totalWeightFwd);

        // Reverse complement strand detection
        std::unordered_map<uint16_t, length_t> detectedTagsRevCompl;
        std::unordered_map<uint16_t, int64_t> lastSeenPositionsRevCompl;
        length_t totalWeightRevCompl = 0;
        mapped |= index.taggingFunction(
            readBundle.getSeqID(), readBundle.getRevComp(), occs, counters,
            REVERSE_C_STRAND, correctTag, minMEMLength, allMEMsunique,
            detectedTagsRevCompl, lastSeenPositionsRevCompl,
            totalWeightRevCompl);

        // Handle strand selection and counter updates
        uint16_t detectedTag =
            handleStrandSelection(counters, detectedTagsFwd, totalWeightFwd,
                                  detectedTagsRevCompl, totalWeightRevCompl,
                                  correctTag, allMEMsunique);


        if (mapped) {
            counters.inc(Counters::MAPPED_READS);
        }

        occs.emplace_back();
        assert(occs.size() == 1);
        std::string detectedTagStr =
            detectedTag == UINT16_MAX ? "*" : std::to_string(detectedTag);
        occs.back().generateTaggingOutputLine(readBundle.getSeqID(),
                                              detectedTagStr);

        return occs;
    }
};

#endif