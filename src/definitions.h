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

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <cstdint>

// ============================================================================
// DEFINITIONS OF CONSTANTS
// ============================================================================

#define VERSION_NUMBER_TAGGER 1     // The version number
#define SUB_VERSION_NUMBER_TAGGER 0 // The sub version number
#define TAGGER_BUILD_INDEX_TAG 1   // the tag for building the index
#define DEFAULT_SEED_LENGTH 0 // The default seed length for the FM-index
#define DEFAULT_MAX_LF 0        // The default maximum LF steps for tag subsampling
#define DEFAULT_L 25        // The default minimum SMEM length
// increment these tags if you change buildhelpers.cpp, buildfmindex.cpp or
// buildbmove.cpps

// ============================================================================
// AT-COMPILE-TIME DEFINITIONS
// ============================================================================

// The type of the length_t depending on compiler option
#if THIRTY_TWO
#define LENGTH_TYPE_NAME "32-bits"
typedef uint32_t length_t;
#else
#define LENGTH_TYPE_NAME "64-bits"
typedef uint64_t length_t;
#endif // THIRTY_TWO

// Check if the compiler is GCC or Clang and if the architecture is 64-bit for
// the use of uint128_t
#if (defined(__GNUC__) || defined(__clang__)) &&                               \
    (defined(__x86_64__) || __SIZEOF_POINTER__ == 8)
#define HAS_UINT128_T 1
#else
#define HAS_UINT128_T 0
#endif // (__GNUC__ || __clang__) && (__x86_64__ || __SIZEOF_POINTER__ == 8)

// ============================================================================
// ENUMS
// ============================================================================

/**
 * An enum for the direction of the search/a substring
 */
enum Direction { FORWARD, BACKWARD };

// An enum for the strand of the read
enum Strand { FORWARD_STRAND = 0, REVERSE_C_STRAND = 1 };

#endif // DEFINITIONS_H
