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
// CLASS TEXT OCC
// ============================================================================
void TextOcc::generateSAMSingleEnd(const string& seqID, const string& printSeq,
                                   const string& printQual, length_t nHits,
                                   length_t minScore, bool primaryAlignment,
                                   const vector<string>& seqNames,
                                   std::string tag) {

    // Estimate the length of the resulting string and reserve memory
    size_t estimated_length = seqID.length() + 100;

    outputLine.reserve(estimated_length);
    size_t _strand = (strand == FORWARD_STRAND) ? 0 : 1;

    outputLine = fmt::format("{}\t{}\t[{},{})\t{}\t{}\n",
                             seqID,            // read name
                             _strand,          // strand as integer
                             range.getBegin(), // range start (begin)
                             range.getEnd(),   // range end
                             distance,         // edit distance
                             tag               // tag
    );
}

void TextOcc::generateSAMSingleEndXA(const string& seqID,
                                     const string& printSeq,
                                     const string& printQual, length_t nHits,
                                     const vector<TextOcc>& otherMatches,
                                     const vector<string>& seqNames) {

    generateSAMSingleEnd(seqID, printSeq, printQual, nHits, distance, true,
                         seqNames);
    // remove \n from end of the line
    outputLine.pop_back();
    // add the X0 X1 and XA tag
    length_t x0 = nHits - 1; // number of co-optimal hits (nHits includes this)
    length_t x1 = otherMatches.size() - x0; // number of suboptimal hits

    // append the sam line
    outputLine += fmt::format("\tX0:i:{}\tX1:i:{}\tXA:Z:", x0, x1);
    for (const auto& m : otherMatches) {
        outputLine += m.asXA(seqNames);
    }
    outputLine += "\n";
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
// CLASS CLUSTER
// ============================================================================

FMOcc Cluster::getClusterCentra(uint16_t lowerBound, vector<FMPosExt>& desc,
                                vector<uint16_t>& initEds) {
    desc.reserve(eds.size());
    initEds.reserve(eds.size());
    FMOcc m;
    for (length_t i = 0; i <= lastCell; i++) {
        if (eds[i] > maxED || eds[i] < lowerBound) {
            continue;
        }
        bool betterThanParent = (i == 0) || eds[i] <= eds[i - 1];
        bool betterThanChild = (i == lastCell) || eds[i] <= eds[i + 1];

        if (betterThanParent && betterThanChild) {
            // this is a valid centre
            nodes[i].report(m, startDepth, eds[i], false, shift);

            // get all the descendants
            initEds.emplace_back(eds[i]);
            for (length_t j = i + 1; j <= lastCell; j++) {
                desc.emplace_back(nodes[j]);
                initEds.emplace_back(eds[j]);
            }

            // replace the clusters under the lower bound
            for (length_t k = 1; k < initEds.size(); k++) {
                if (initEds[k] < lowerBound && initEds[k] <= initEds[k - 1] &&
                    (k == initEds.size() - 1 || initEds[k] <= initEds[k + 1])) {
                    // k is a centre under the lower bound

                    length_t highestPoint = 0;
                    length_t lowestPoint = initEds.size() - 1;
                    // find highest point of this cluster
                    for (length_t l = k; l-- > 0;) {
                        if (initEds[l] != initEds[l + 1] + 1) {
                            highestPoint = l + 1;
                            break;
                        }
                    }
                    // find lowest point of this cluster
                    for (length_t l = k + 1; l < initEds.size(); l++) {
                        if (initEds[l] != initEds[l - 1] + 1) {
                            lowestPoint = l - 1;
                            break;
                        }
                    }

                    // highest and lowest cannot span entire
                    // initEds.size(), otherwise there would not be a
                    // valid cluster centre above the lower bound
                    if (highestPoint != 0 &&
                        lowestPoint != initEds.size() - 1) {
                        // Make /\ with ed values of this cluster
                        // do iE[hp] = ie[hp - 1] + 1 and iE[lp] = iE[lp
                        // + 1] +1 until entire cluster has been
                        // replaced
                        length_t lC = lowestPoint;
                        length_t hC = highestPoint;
                        bool highest = true;
                        // do not go over maxED + 1, to ensure
                        // continuity at the other end
                        while (lC > hC) {
                            if (highest) {
                                initEds[hC] =
                                    min(maxED + 1, initEds[hC - 1] + 1);
                                hC++;
                            } else {
                                initEds[lC] =
                                    min(maxED + 1, initEds[lC + 1] + 1);
                                lC--;
                            }
                            highest = !highest;
                        }
                        if (lC == hC) {
                            // change middle element of cluster
                            initEds[lC] =
                                min(initEds[lC + 1] + 1, initEds[lC - 1] + 1);
                        }

                    } else if (highestPoint == 0 &&
                               lowestPoint != initEds.size() - 1) {
                        // monotonous rise from lowestPoint to
                        // highestPoint
                        for (length_t l = lowestPoint; l-- > 0;) {
                            initEds[l] = initEds[l + 1] + 1;
                        }
                    } else if (highestPoint != 0 &&
                               lowestPoint == initEds.size() - 1) {
                        // monotonous rise from highestPoint to
                        // lowestPoint
                        for (length_t l = highestPoint; l < initEds.size();
                             l++) {
                            initEds[l] = initEds[l - 1] + 1;
                        }
                    }
                }
            }
            // stop searching
            break;
        }
    }

    return m;
}

// ============================================================================
// CLASS COUNTERS
// ============================================================================

void Counters::reportStatistics(const SequencingMode& sMode) const {
    std::stringstream ss;

    ss << "Average no. nodes: "
       << counters[NODE_COUNTER] / (counters[NUMBER_OF_READS] * 1.0);
    logger.logDeveloper(ss);

    ss << "Total no. Nodes: " << counters[NODE_COUNTER];
    logger.logDeveloper(ss);

    if (sMode == SINGLE_END) {
        ss << "Average no. unique matches per read: "
           << counters[TOTAL_UNIQUE_MATCHES] /
                  (counters[NUMBER_OF_READS] * 1.0);
        logger.logDeveloper(ss);

        ss << "Total no. matches: " << counters[TOTAL_UNIQUE_MATCHES];
        logger.logDeveloper(ss);

        ss << "Average no. matches per read "
           << counters[TOTAL_REPORTED_POSITIONS] /
                  (counters[NUMBER_OF_READS] * 1.0);
        logger.logDeveloper(ss);

        ss << "Total no. reported matches: "
           << counters[TOTAL_REPORTED_POSITIONS];
        logger.logDeveloper(ss);

        ss << "Mapped reads: " << counters[MAPPED_READS];
        logger.logDeveloper(ss);

        ss << "Number of reads: " << counters[NUMBER_OF_READS];
        logger.logDeveloper(ss);

        ss << "Percentage reads mapped: "
           << (counters[MAPPED_READS] * 100.0) / counters[NUMBER_OF_READS]
           << "%";
        logger.logDeveloper(ss);
    } else {
        // Paired end
        ss << "Average no. matches per pair: "
           << counters[TOTAL_UNIQUE_PAIRS] / (counters[NUMBER_OF_READS] / 2.0);
        logger.logInfo(ss);

        ss << "Total no. matches: " << counters[TOTAL_UNIQUE_PAIRS];
        logger.logInfo(ss);

        ss << "Mapped pairs: " << counters[MAPPED_PAIRS];
        logger.logInfo(ss);

        ss << "Percentage of pairs mapped: "
           << (counters[MAPPED_PAIRS] * 100.0) / (counters[NUMBER_OF_READS] / 2)
           << "%";
        logger.logInfo(ss);

        ss << "Discordantly mapped pairs: "
           << counters[DISCORDANTLY_MAPPED_PAIRS];
        logger.logInfo(ss);

        ss << "Percentage of discordantly mapped pairs: "
           << (counters[DISCORDANTLY_MAPPED_PAIRS] * 100.0) /
                  (counters[NUMBER_OF_READS] / 2)
           << "%";
        logger.logInfo(ss);

        ss << "No. unpaired reads that did match: "
           << counters[MAPPED_HALF_PAIRS];
        logger.logInfo(ss);

        ss << "Percentage pairs for which only 1 read matched: "
           << (counters[MAPPED_HALF_PAIRS] * 100.0) /
                  (counters[NUMBER_OF_READS] / 2)
           << "%";
        logger.logInfo(ss);

        ss << "Total read pairs both mapped but unpaired: "
           << counters[UNPAIRED_BUT_MAPPED_PAIRS];
        logger.logInfo(ss);
    }

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