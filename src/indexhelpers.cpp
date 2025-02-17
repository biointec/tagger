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

#include "indexhelpers.h"
#include "logger.h"     // for logger
#include <cstdint>      // for uint16_t, uint64_t
#include <fmt/core.h>   // for format
#include <fmt/format.h> // for to_string
#include <iomanip>      // For std::setprecision

using namespace std;

// ============================================================================
// CLASS RANGE
// ============================================================================

ostream& operator<<(ostream& os, const Range& r) {
    os << "[" << r.begin << ", " << r.end << ")";
    return os;
}

// ============================================================================
// CLASS MOVERANGE
// ============================================================================

ostream& operator<<(ostream& os, const MoveRange& r) {
    os << "[" << r.begin << ", " << r.end << ") in runs [" << r.beginRun << ", "
       << r.endRun << "]";
    return os;
}

// ============================================================================
// CLASS FM OCC
// ============================================================================

ostream& operator<<(ostream& o, const FMOcc& m) {
    return o << "SARange: " << m.getRanges().getRangeSA()
             << "\tEdit distance: " << m.getDistance()
             << "\tdepth: " << m.getDepth();
}

// ============================================================================
// CLASS COUNTERS
// ============================================================================

void Counters::reportStatistics() const {
    std::stringstream ss;

    ss << "Average no. nodes: "
       << counters[NODE_COUNTER] / (counters[NUMBER_OF_READS] * 1.0);
    logger.logDeveloper(ss);

    ss << "Total no. Nodes: " << counters[NODE_COUNTER];
    logger.logDeveloper(ss);

    ss << "Average no. unique matches per read: "
       << counters[TOTAL_UNIQUE_MATCHES] / (counters[NUMBER_OF_READS] * 1.0);
    logger.logDeveloper(ss);

    ss << "Total no. matches: " << counters[TOTAL_UNIQUE_MATCHES];
    logger.logDeveloper(ss);

    ss << "Average no. matches per read "
       << counters[TOTAL_REPORTED_POSITIONS] /
              (counters[NUMBER_OF_READS] * 1.0);
    logger.logDeveloper(ss);

    ss << "Total no. reported matches: " << counters[TOTAL_REPORTED_POSITIONS];
    logger.logDeveloper(ss);

    ss << "Mapped reads: " << counters[MAPPED_READS];
    logger.logDeveloper(ss);

    ss << "Number of reads: " << counters[NUMBER_OF_READS];
    logger.logDeveloper(ss);

    ss << "Percentage reads mapped: "
       << (counters[MAPPED_READS] * 100.0) / counters[NUMBER_OF_READS] << "%";
    logger.logDeveloper(ss);

#ifdef GROUND_TRUTH_CHECKS
    logger.logInfo("Tagging statistics:");
    logger.logInfo(
        "Number of unmapped reads: " +
        std::to_string(counters[NUMBER_OF_READS] - counters[MAPPED_READS]));
    logger.logInfo(
        "Percentage of unmapped reads: " +
        std::to_string(
            ((counters[NUMBER_OF_READS] - counters[MAPPED_READS]) * 100.0) /
            counters[NUMBER_OF_READS]) +
        "%");

    logger.logInfo("If all MEMs must be unique:");
    logger.logInfo("Number of correctly identified reads: " +
                   std::to_string(counters[TOTAL_CORRECTLY_IDENTIFIED_READS]));
    logger.logInfo(
        "Percentage of correctly identified reads: " +
        std::to_string((counters[TOTAL_CORRECTLY_IDENTIFIED_READS] * 100.0) /
                       counters[NUMBER_OF_READS]) +
        "%");
    logger.logInfo(
        "Number of unidentified reads: " +
        std::to_string(counters[TOTAL_UNIDENTIFIED_READS] -
                       (counters[NUMBER_OF_READS] - counters[MAPPED_READS])));
    logger.logInfo(
        "Percentage of unidentified reads: " +
        std::to_string(((counters[TOTAL_UNIDENTIFIED_READS] -
                         (counters[NUMBER_OF_READS] - counters[MAPPED_READS])) *
                        100.0) /
                       counters[NUMBER_OF_READS]) +
        "%");
    logger.logInfo(
        "Number of incorrectly identified reads: " +
        std::to_string(counters[TOTAL_INCORRECTLY_IDENTIFIED_READS]));
    logger.logInfo(
        "Percentage of incorrectly identified reads: " +
        std::to_string((counters[TOTAL_INCORRECTLY_IDENTIFIED_READS] * 100.0) /
                       counters[NUMBER_OF_READS]) +
        "%");

    logger.logInfo("If not all MEMs must be unique:");
    logger.logInfo("Number of correctly identified reads: " +
                   std::to_string(counters[TOTAL_CORRECTLY_IDENTIFIED_READS2]));
    logger.logInfo(
        "Percentage of correctly identified reads: " +
        std::to_string((counters[TOTAL_CORRECTLY_IDENTIFIED_READS2] * 100.0) /
                       counters[NUMBER_OF_READS]) +
        "%");
    logger.logInfo(
        "Number of unidentified reads: " +
        std::to_string(counters[TOTAL_UNIDENTIFIED_READS2] -
                       (counters[NUMBER_OF_READS] - counters[MAPPED_READS])));
    logger.logInfo(
        "Percentage of unidentified reads: " +
        std::to_string(((counters[TOTAL_UNIDENTIFIED_READS2] -
                         (counters[NUMBER_OF_READS] - counters[MAPPED_READS])) *
                        100.0) /
                       counters[NUMBER_OF_READS]) +
        "%");
    logger.logInfo(
        "Number of incorrectly identified reads: " +
        std::to_string(counters[TOTAL_INCORRECTLY_IDENTIFIED_READS2]));
    logger.logInfo(
        "Percentage of incorrectly identified reads: " +
        std::to_string((counters[TOTAL_INCORRECTLY_IDENTIFIED_READS2] * 100.0) /
                       counters[NUMBER_OF_READS]) +
        "%");

    logger.logInfo("Using the consensus method:");
    logger.logInfo("Number of correctly identified reads: " +
                   std::to_string(counters[TOTAL_CORRECTLY_IDENTIFIED_READS3]));
    logger.logInfo(
        "Percentage of correctly identified reads: " +
        std::to_string((counters[TOTAL_CORRECTLY_IDENTIFIED_READS3] * 100.0) /
                       counters[NUMBER_OF_READS]) +
        "%");
    logger.logInfo(
        "Number of unidentified reads: " +
        std::to_string(counters[TOTAL_UNIDENTIFIED_READS3] -
                       (counters[NUMBER_OF_READS] - counters[MAPPED_READS])));
    logger.logInfo(
        "Percentage of unidentified reads: " +
        std::to_string(((counters[TOTAL_UNIDENTIFIED_READS3] -
                         (counters[NUMBER_OF_READS] - counters[MAPPED_READS])) *
                        100.0) /
                       counters[NUMBER_OF_READS]) +
        "%");
    logger.logInfo(
        "Number of incorrectly identified reads: " +
        std::to_string(counters[TOTAL_INCORRECTLY_IDENTIFIED_READS3]));
    logger.logInfo(
        "Percentage of incorrectly identified reads: " +
        std::to_string((counters[TOTAL_INCORRECTLY_IDENTIFIED_READS3] * 100.0) /
                       counters[NUMBER_OF_READS]) +
        "%");

    logger.logInfo("Average read length: " +
                   std::to_string(counters[TOTAL_READ_LENGTHS] /
                                  (counters[NUMBER_OF_READS])));
    logger.logInfo("Average number of steps per read: " +
                   std::to_string(counters[TOTAL_STEP_COUNTS] /
                                  (counters[NUMBER_OF_READS] * 2.0)));
#endif // GROUND_TRUTH_CHECKS
}