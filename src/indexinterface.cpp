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

#include "indexinterface.h"
#include "definitions.h"
#include "indexhelpers.h"
#include "logger.h"

#include <assert.h>                 // for assert
#include <cstdint>                  // for uint16_t
#include <fmt/format.h>             // for fmt::to_string
#include <fstream>                  // for operator<<, ifstream, basic_ostream
#include <iterator>                 // for distance
#include <limits>                   // for numeric_limits
#include <memory>                   // for allocator_traits<>::value_type
#include <parallel_hashmap/phmap.h> // phmap::parallel_flat_hash_map
#include <stdexcept>                // for runtime_error
#include <type_traits>              // for __strip_reference_wrapper<>::__type

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace std;

// ============================================================================
// CLASS IndexInterface
// ============================================================================

thread_local Direction IndexInterface::dir = BACKWARD;
thread_local ExtraCharPtr IndexInterface::extraChar;

// ----------------------------------------------------------------------------
// PREPROCESSING ROUTINES
// ----------------------------------------------------------------------------

void IndexInterface::readMetaAndCounts(const string& baseFile, bool verbose) {

    // Open the meta file
    ifstream ifs(baseFile + ".meta");
    if (ifs) {
        // Read the build tag
        length_t tag;
        ifs >> tag;
        if (tag < TAGGER_BUILD_INDEX_TAG) {
            logger.logWarning(
                "The index was built with an older version of tagger_build. "
                "Proceed with caution or re-build the index.");
        } else if (tag > TAGGER_BUILD_INDEX_TAG) {
            logger.logDeveloper(
                "The index was built with a newer version of tagger_build. "
                "Proceed with caution or re-build the index.");
        }
        // read the compiled info
        size_t build_length_t_size;
        ifs >> build_length_t_size;
        if (build_length_t_size != sizeof(length_t)) {
            size_t thisSize = sizeof(length_t) * 8;
            build_length_t_size *= 8;
            throw runtime_error(
                "The index was built with a compiled version that uses " +
                to_string(build_length_t_size) +
                "-bit numbers, while the current programme was compiled "
                "using " +
                to_string(thisSize) +
                "-bit numbers. Recompile the programme with the correct "
                "THIRTY_TWO flag set or rebuild the index.");
        }

    } else {
        logger.logWarning(
            "The index was built with an older version of tagger_build or "
            "the meta data is missing! "
            "Proceed with caution or re-build the index.");
    }

    stringstream ss;
    if (verbose) {
        ss << "Reading " << baseFile << ".cct" << "...";
        logger.logInfo(ss);
    }

    // Read the counts table
    vector<length_t> charCounts(256, 0);
    if (!readArray(baseFile + ".cct", charCounts)) {
        throw runtime_error("Cannot open file: " + baseFile + ".cct");
    }

    length_t cumCount = 0; // Cumulative character counts
    for (size_t i = 0; i < charCounts.size(); i++) {
        if (charCounts[i] == 0)
            continue;
        cumCount += charCounts[i];
    }
    textLength = cumCount;
    sigma = Alphabet<ALPHABET>(charCounts);
}

bool IndexInterface::readArray(const string& filename,
                               vector<length_t>& array) {
    int fd = open(filename.c_str(), O_RDONLY);
    if (fd == -1) {
        // Handle error
        return false;
    }

    struct stat sb;
    if (fstat(fd, &sb) == -1) {
        close(fd);
        return false;
    }

    size_t fileSize = sb.st_size;
    if (fileSize % sizeof(length_t) != 0) {
        close(fd);
        return false;
    }

    size_t numElements = fileSize / sizeof(length_t);

    void* mapped = mmap(NULL, fileSize, PROT_READ, MAP_PRIVATE, fd, 0);
    if (mapped == MAP_FAILED) {
        close(fd);
        return false;
    }

    length_t* data = static_cast<length_t*>(mapped);
    try {
        array.assign(data, data + numElements);
    } catch (...) {
        munmap(mapped, fileSize);
        throw;
    }

    munmap(mapped, fileSize);
    close(fd);

    return true;
}

void IndexInterface::populateTable(bool verbose) {
    if (verbose) {
        stringstream ss;
        ss << "Populating FM-range table with " << wordSize << "-mers...";
        logger.logInfo(ss);
    }

    Kmer::setWordSize(wordSize);

    table.reserve(1 << (2 * wordSize)); // 2 << wordSize is 4^wordSize (DNA)
    setDirection(BACKWARD); // backward for toehold compatibility

    string word;
    vector<FMPosExt> stack;
    Counters counters;
    extendFMPos(getCompleteRange(), stack, counters);
    while (!stack.empty()) {
        auto curr = stack.back();
        stack.pop_back();

        word.resize(curr.getRow());
        word[curr.getRow() - 1] = curr.getCharacter();

        if (curr.getRow() == wordSize) { // max depth reached
            std::reverse(word.begin(),
                         word.end()); // Reverse for Kmer construction
            Kmer k(word);
            updateRangeSARevRuns(curr.getRangesMutable());
            table.insert(make_pair(k, std::move(curr.getRanges())));
            std::reverse(word.begin(),
                         word.end()); // Restore reverse order for reuse

        } else // add extra characters
            extendFMPos(curr.getRanges(), stack, counters, curr.getRow());
    }

    setDirection(FORWARD); // reset direction
}

// ----------------------------------------------------------------------------
// EXTEND ROUTINES
// ----------------------------------------------------------------------------

void IndexInterface::extendFMPos(const SARangePair& parentRanges,
                                 vector<FMPosExt>& stack, Counters& counters,
                                 length_t row) const {
    // iterate over the entire alphabet
    for (length_t i = 1; i < sigma.size() - 1; ++i) { // skip '$' and 'X'

        SARangePair pairForNewChar;

        // check if this character occurs in the specified range
        if ((this->*extraChar)(i, parentRanges, pairForNewChar)) {
            // push this range and character for the next iteration
            stack.emplace_back(sigma.i2c(i), pairForNewChar, row + 1);

            counters.inc(Counters::NODE_COUNTER);
        }
    }
}

// ----------------------------------------------------------------------------
// ROUTINES FOR EXACT PATTERN MATCHING
// ----------------------------------------------------------------------------

bool IndexInterface::exactMatchesOutput(const string& s, Counters& counters,
                                        vector<TextOcc>& tOcc) const {
    if (s.size() == 0) {
        return false;
    }

    SARangePair completeRanges = getCompleteRange();
    SARangeBackwards range = SARangeBackwards(
        completeRanges.getRangeSA(), completeRanges.getToehold(),
        completeRanges.getToeholdRepresentsEnd(),
        completeRanges.getOriginalDepth());

    length_t i = s.size();

    for (; i-- > 0;) {
        const auto& c = s[i];
        auto pos = sigma.c2i(c);
        if (pos == -1 || !findRangeWithExtraCharBackward(pos, range, range)) {
            // c is not in alphabet or no updated range is empty, exact match
            // impossible
            return false;
        }
        counters.inc(Counters::NODE_COUNTER);
    }

    tOcc.emplace_back(range, 0, FORWARD_STRAND);

    counters.inc(Counters::TOTAL_REPORTED_POSITIONS, tOcc.size());
    return tOcc.size() > 0;
}