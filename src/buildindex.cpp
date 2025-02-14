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

#include "alphabet.h"    // for Alphabet
#include "bitvec.h"      // for Bitvec, Bitref
#include "definitions.h" // for length_t, TAGGER_BUILD_INDEX_TAG, ...
#include "logger.h"      // for Logger, logger
#include "util.h"        // for fileExists..

#include <algorithm>  // for max, transform, equal, max_element, min_ele...
#include <array>      // for array
#include <cctype>     // for toupper
#include <cstdint>    // for int64_t, uint8_t, uint32_t
#include <functional> // for function
#include <iostream>   // for ifstream, ofstream, operator<<, ios, basic_...
#include <limits>     // for numeric_limits
#include <random>     // for minstd_rand
#include <stdexcept>  // for runtime_error
#include <stdlib.h>   // for size_t
#include <string>     // for string
#include <vector>     // for vector

#ifdef THIRTY_TWO
#include "libsais.h" // for libsais
#include "libsais64.h"
#else
#include "divsufsort64.h"
#endif

#ifdef HAVE_ZLIB
#include "seqfile.h"
#include <zlib.h> // for gzFile, gzopen, gzclose, gzread, gzwrite
#endif

#ifdef BIG_BWT_USABLE
#include "../external/Big-BWT/utils.h" // for SABYTES
#endif
#include "bmove/moverepr.h"     // for MoveLFRepr
#include "bmove/sparsebitvec.h" // for SparseB...
#include "sdsl/int_vector.hpp"  // for int_vector

using namespace std;

// ============================================================================
// PARSING
// ============================================================================

#include "parameters/buildparameters.h" // for BuildParameters

// ============================================================================
// FUNCTIONS FOR BOTH VANILLA AND RLC
// ============================================================================

bool foundU = false;

/**
 * @brief Replace non-ACGT characters with a random ACGT character.
 * @param original The original character.
 * @param gen The random number generator.
 * @param seed The seed for the random number generator. (dummy variable)
 * @param seedIndex The index of the seed. (dummy variable)
 * @return The replaced character.
 */
char replaceNonACGT(char original, std::minstd_rand& gen,
                    const std::string& seed, size_t& seedIndex) {
    const std::string validChars = "ACGT";
    if (original != 'A' && original != 'C' && original != 'G' &&
        original != 'T') {
        if (!foundU && original == 'U') {
            foundU = true;
            logger.logWarning("Found 'U' in reference text. Is this an RNA "
                              "reference? 'U' will be replaced with 'X'. "
                              "This warning will only be shown once.");
        }
        return 'X'; // Replace with 'N'
    }
    return original;
}

bool containsCaseInsensitive(const std::string& haystack,
                             const std::string& needle) {
    if (needle.empty())
        return true;
    if (haystack.size() < needle.size())
        return false;

    auto it =
        std::search(haystack.begin(), haystack.end(), needle.begin(),
                    needle.end(), [](char ch1, char ch2) {
                        return std::tolower(static_cast<unsigned char>(ch1)) ==
                               std::tolower(static_cast<unsigned char>(ch2));
                    });

    return it != haystack.end();
}

bool isGzipped(FileType fileType, const string& filename) {
    if (fileType == FileType::FASTA_GZ) {
#ifndef HAVE_ZLIB
        // If the file is gzipped but zlib support is not available, throw an
        // error
        throw runtime_error("Error: " + fastaFile +
                            " is a gzipped FASTA file, but tagger was not "
                            "compiled with zlib support.");
#else
        // return true if the file is gzipped and zlib support is
        // available
        return true;
#endif
    }
    return false;
}

SeqFile openFastaFile(const string& fastaFile) {
    // Determine the file type of the provided FASTA file
    FileType fileType;
    tie(fileType, ignore) = getFileType(fastaFile);

    // Validate that the file is either a FASTA or gzipped FASTA file
    if (fileType != FileType::FASTA && fileType != FileType::FASTA_GZ) {
        throw runtime_error("Error: " + fastaFile +
                            " is not a valid FASTA file.");
    }

    // Validate that the FASTA file exists
    if (!Util::fileExists(fastaFile)) {
        throw runtime_error("Error: " + fastaFile + " does not exist.");
    }

    // Create a SeqFile object and open the FASTA file
    SeqFile file(isGzipped(fileType, fastaFile));
    file.open(fastaFile);

    return file;
}

/**
 * @brief Concatenates and transforms the sequences from a FASTA file.
 *
 * This function reads a FASTA file and concatenates the sequences into a single
 * std::string. It also replaces non-ACGT characters with a random ACGT
 * character. The start positions of each sequence and the sequence names are
 * stored in separate std::vectors. The resulting concatenated std::string is
 * terminated with a '$' character.
 *
 * @param fastaFile The path to the FASTA file.
 * @param concatenation The resulting concatenated std::string. (output)
 * @param positions The std::vector to store the start positions of each
 * sequence. (output)
 * @param seqNames The std::vector to store the sequence names. (output)
 * @param replaceFunc The function to use for replacing non-ACGT characters.
 * @param expectedNumber The expected number of sequences in the FASTA file.
 * Default is 12000.
 */
void concatenateAndTransform(const std::string& fastaFile,
                             std::string& concatenation,
                             std::vector<length_t>& positions,
                             std::vector<std::string>& seqNames,
                             std::vector<length_t>& tags,
                             std::vector<std::string>& taggingCategoriesVec,
                             std::vector<length_t>& firstSeqIDPerFile,
                             std::function<char(char, size_t&)> replaceFunc,
                             length_t expectedNumber = 12000) {

    SeqFile file = openFastaFile(fastaFile);

    // Reserve memory for the concatenation string based on estimated file size
    concatenation.reserve(
        static_cast<size_t>(file.estimateASCIISize() + concatenation.size()));

    // Reserve memory for positions and seqNames vectors to avoid reallocations
    positions.reserve(expectedNumber + positions.size());
    seqNames.reserve(expectedNumber + seqNames.size());
    tags.reserve(expectedNumber + tags.size());

    // Record the index of the first sequence ID in the file
    firstSeqIDPerFile.push_back(seqNames.size());

    // Temporary variables for processing sequence data
    std::string sequence;      // To hold the current sequence
    bool sequenceName = false; // Flag to indicate if a seq name was present
    std::string line;          // To hold each line read from the file
    length_t startPosition = concatenation.size(); // start position in concat.
    size_t seedIndex = 0; // Index for seed to replace non-ACGT characters

    while (file.good()) {
        file.getLine(line);

        if (line.empty() || (line.size() == 1 && line[0] == '\n'))
            continue; // Skip empty lines

        // pop back the new line character
        if (line.back() == '\n')
            line.pop_back();

        if (line[0] == '>') { // Sequence name line

            // Process the previous sequence
            if (!sequence.empty()) {
                // Process the previous sequence
                for (char& c : sequence) {
                    c = replaceFunc(std::toupper(c), seedIndex);
                    concatenation += c;
                }

                concatenation += 'X'; // Add a separator character

                positions.emplace_back(
                    startPosition); // push back the start position of this
                                    // sequence
                startPosition +=
                    sequence
                        .length(); // create start position for next sequence
                sequence.clear();  // Clear the sequence for the next one
            }
            std::string description = line.substr(1);
            // get the sequence name (before the first space)
            std::string name = description.substr(0, description.find(' '));

            // Check if the sequence name is already present
            if (std::find(seqNames.begin(), seqNames.end(), name) !=
                seqNames.end()) {
                throw std::runtime_error("Error: Sequence name " + name +
                                         " in file " + fastaFile +
                                         " is not unique!");
            }

            seqNames.emplace_back(name); // Store sequence name
            sequenceName = true;

            if (taggingCategoriesVec.size() > 0) {
                // Iterate over the possible tagging categories
                bool found = false;
                for (size_t i = 0; i < taggingCategoriesVec.size(); i++) {
                    const std::string& tag = taggingCategoriesVec[i];
                    // Check if the description contains the tag without case
                    // sensitivity
                    if (containsCaseInsensitive(description, tag)) {
                        if (found) {
                            logger.logError("Reference " + description +
                                            " contains multiple tags from the "
                                            "tagging categories file");
                            throw std::runtime_error(
                                "Reference " + description +
                                " contains multiple tags from "
                                "the tagging categories file");
                        }
                        tags.push_back(i);
                        found = true;
                    }
                }
                if (!found) {
                    logger.logError("Reference " + description +
                                    " does not contain any tags from the "
                                    "tagging categories file");
                    throw std::runtime_error("Reference " + description +
                                             " does not contain any tags from "
                                             "the tagging categories file");
                }
            }
        } else {
            // Append to the current sequence
            sequence += line;
        }
    }

    // Process the last sequence in the file
    for (char& c : sequence) {
        c = replaceFunc(std::toupper(c), seedIndex);
        concatenation += c;
    }
    positions.emplace_back(startPosition);

    if (!sequenceName) {
        seqNames.emplace_back(fastaFile);
    }

    // Close input file
    file.close();
}

/**
 * @brief Read the contents of a text file into a std::string buffer.
 *
 */
void readText(const string& filename, string& buf) {
    ifstream ifs(filename);
    if (!ifs)
        throw runtime_error("Cannot open file: " + filename);

    ifs.seekg(0, ios::end);
    buf.resize(ifs.tellg());
    ifs.seekg(0, ios::beg);
    ifs.read((char*)buf.data(), buf.size());
}

/**
 * @brief Perform a sanity check on the suffix array.
 * @param T The text.
 * @param sa The suffix array.
 */
void sanityCheck(const string& T, vector<length_t>& sa) {
    logger.logInfo("Performing sanity checks...");
    // check T for correctness
    if (T.back() == '\n')
        throw runtime_error("T should end with a \'$\' character, "
                            "not with a newline");

    if (T.back() != '$')
        throw runtime_error("T should end with a \'$\' character");

    if (sa.size() != T.size())
        throw runtime_error("Text and suffix array contain a "
                            "different number of elements");

    // briefly check the suffix array
    length_t min = *min_element(sa.begin(), sa.end());
    length_t max = *max_element(sa.begin(), sa.end());

    if (min == 1 && max == T.size()) { // rebase to [0..T.size()-1]
        for (auto& el : sa)
            el--;
        min--;
        max--;
    }

    if (min != 0 || max != T.size() - 1)
        throw runtime_error("Suffix array must contain numbers between "
                            "[0 and " +
                            to_string(T.size() - 1) + "]");

    // check if all numbers in the suffix array are present
    Bitvec bv(sa.size());
    for (length_t i : sa)
        bv[i] = true;

    for (size_t i = 0; i < bv.size(); i++)
        if (!bv[i])
            throw runtime_error("Suffix " + to_string(i) +
                                " seems "
                                "to be missing from suffix array");

    // extra check:
    //      we could check T to see if the SA correctly sorts suffixes of T

    logger.logInfo("\tSanity checks OK");
}

/**
 * @brief Create and write the header information for output SAM files.
 * @param baseFN The base filename.
 * @param seqNames The sequence names.
 * @param positions The start positions of each sequence.
 */
void createAndWriteHeaderInfo(const string& baseFN,
                              const vector<string>& seqNames,
                              const vector<length_t>& positions) {
    logger.logInfo("Writing SAM header info for reference text to " + baseFN +
                   ".headerSN.bin...");
    // create the header lines with these reference sequences
    std::ofstream headerStream(baseFN + ".headerSN.bin", ios::binary);
    for (length_t i = 0; i < seqNames.size(); i++) {
        headerStream << "@SQ\tSN:" << seqNames[i]
                     << "\tLN:" << positions[i + 1] - positions[i] << "\n";
    }
    headerStream.close();
}

/**
 * @brief Write the start positions and sequence names to disk.
 * @param baseFN The base filename.
 * @param positions The start positions of each sequence.
 * @param seqNames The sequence names.
 */
void writePositionsAndSequenceNames(const string& baseFN,
                                    const vector<length_t>& positions,
                                    const vector<string>& seqNames,
                                    const vector<length_t>& firstSeqIDPerFile) {
    logger.logInfo("Write positions to " + baseFN + ".pos...");
    // Write the positions to disk
    std::ofstream ofs2(baseFN + ".pos", ios::binary);
    ofs2.write((char*)positions.data(), positions.size() * sizeof(length_t));
    ofs2.close();
}

/**
 * @brief Check if the text size is within the limits of the length_t type.
 * @param tSize The size of the text.
 * @param checkBothLimits Check both the lower and upper limits. Default is
 * true. If false, only the upper limit is checked.
 * @throws runtime_error if the text size is too large.
 */
void checkTextSize(const size_t tSize, bool checkBothLimits = true) {
    if (tSize > std::numeric_limits<length_t>::max()) {
        throw std::runtime_error(
            "The size of the reference (" + std::to_string(tSize) +
            ") is too large. Maximum size for the current "
            "compiler option is: " +
            std::to_string(std::numeric_limits<length_t>::max()) +
            "\nPlease recompile tagger with a larger word size.");
    }

    if (checkBothLimits && (sizeof(length_t) * 8 == 64) &&
        (tSize <= std::numeric_limits<uint32_t>::max())) {
        logger.logWarning(
            "Program was compiled with 64-bit words, but the text size (" +
            std::to_string(tSize) +
            ") fits in a 32-bit word. Consider recompiling the program to "
            "improve performance and save memory.");
    }
}

/**
 * @brief Count the frequency of each character in the text.
 * @param T The text.
 * @param charCounts The vector to store the character counts.
 */
void countChars(const string& T, vector<length_t>& charCounts) {
    // initialize the counts to 0 and 256 elements
    charCounts = vector<length_t>(256, 0);
    for (char c : T)
        charCounts[(unsigned char)c]++;
}

/**
 * @brief Write the character counts to disk.
 * @param baseFN The base filename.
 * @param charCounts The character counts.
 */
void writeCharCounts(const string& baseFN, const vector<length_t>& charCounts) {
    ofstream ofs(baseFN + ".cct", ios::binary);
    ofs.write((char*)charCounts.data(), charCounts.size() * sizeof(length_t));
    ofs.close();
    logger.logInfo("Wrote file " + baseFN + ".cct");
}

/**
 * @brief Create the alphabet based on the character counts.
 * @param T The text.
 * @param charCounts The character counts.
 * @param sigma The alphabet. (output)
 */
void createAlphabet(const string& T, const vector<length_t>& charCounts,
                    Alphabet<ALPHABET>& sigma) {
    // count the number of unique characters in T
    int nUniqueChar = 0;
    for (length_t count : charCounts)
        if (count > 0)
            nUniqueChar++;

    logger.logInfo("Text has length " + std::to_string(T.size()));
    logger.logInfo("Text has " + std::to_string(nUniqueChar) +
                   " unique characters");

    if (T.size() == 1) {
        logger.logError(
            "Reference text is empty except for sentinel character. Please "
            "provide a valid reference with the -f or -F flag.");
        exit(EXIT_FAILURE);
    }

    if (nUniqueChar > ALPHABET) {
        logger.logError("FATAL: the number of unique characters in the "
                        "text exceeds the alphabet size. Please recompile "
                        "tagger using a higher value for ALPHABET");
        exit(EXIT_FAILURE);
    }

    if (nUniqueChar < ALPHABET) {
        logger.logWarning("the number of unique characters in the "
                          "text is less than the ALPHABET size specified when "
                          "tagger was compiled. Performance may be affected.");
    }

    sigma = Alphabet<ALPHABET>(charCounts);
}

#ifdef THIRTY_TWO
void createSuffixArrayLibsais(const string& T, vector<uint32_t>& SA) {
    logger.logInfo("Generating the suffix array using libsais...");

    // Convert std::string to const uint8_t*
    const uint8_t* tPtr = reinterpret_cast<const uint8_t*>(T.c_str());

    if (T.size() > INT32_MAX) {
        // We need libsais64 to handle this
        // create a temporary vector to store the 64-bit suffix array
        SA.clear();
        vector<int64_t> SA64(T.size());
        logger.logDeveloper("Calling libsais64");
        int64_t result = libsais64(tPtr, SA64.data(),
                                   static_cast<int64_t>(T.size()), 0, nullptr);
        if (result != 0) {
            throw runtime_error("libsais64 failed with error code " +
                                to_string(result));
        }
        SA.reserve(T.size());
        // convert the 64-bit suffix array to 32-bit
        for (size_t i = 0; i < T.size(); i++) {
            SA.emplace_back(static_cast<uint32_t>(SA64[i]));
        }
    } else {
        // Resize SA to the size of the input string to hold the suffix array
        SA.resize(T.size());
        logger.logDeveloper("Calling libsais");
        int32_t result = libsais(tPtr, reinterpret_cast<int32_t*>(SA.data()),
                                 static_cast<int32_t>(T.size()), 0, nullptr);

        if (result != 0) {
            throw runtime_error("libsais failed with error code " +
                                to_string(result));
        }
    }

    logger.logDeveloper("Libsais successful");
}
#else
void createSuffixArrayDivsufsort64(const string& T, vector<uint64_t>& SA) {
    logger.logInfo("Generating the suffix array using divsufsort...");

    // Resize SA to the size of the input string to hold the suffix array
    SA.resize(T.size());

    // Convert std::string to const uint8_t*
    const sauchar_t* tPtr = reinterpret_cast<const sauchar_t*>(T.c_str());

    logger.logDeveloper("Calling divsufsort64");
    saint_t result = divsufsort64(tPtr, reinterpret_cast<int64_t*>(SA.data()),
                                  static_cast<int64_t>(T.size()));

    if (result != 0) {
        throw std::runtime_error("divsufsort failed with error code " +
                                 std::to_string(result));
    }

    logger.logDeveloper("Divsufsort successful");
}
#endif

/**
 * @brief Create the suffix array using the libsais library.
 * @param T The text.
 * @param SA The suffix array. (output)
 */
void createSuffixArray(const string& T, vector<length_t>& SA) {
#ifdef THIRTY_TWO
    createSuffixArrayLibsais(T, SA);
#else
    createSuffixArrayDivsufsort64(T, SA);

#endif

    logger.logInfo("Suffix array generated successfully!");
}

/**
 * @brief Create the BWT of the reversed text.
 * @param baseFN The base filename.
 * @param T The text.
 * @param revSA The suffix array of the reversed text.
 * @param sigma The alphabet.
 * @param rBWT The BWT of the reversed text. (output)
 */
void createRevBWT(const string& baseFN, const string& T,
                  const vector<length_t>& revSA,
                  const Alphabet<ALPHABET>& sigma, string& rBWT) {

    // build the reverse BWT
    logger.logInfo("Generating BWT of reversed text...");
    rBWT.resize(T.size());
    for (size_t i = 0; i < revSA.size(); i++)
        rBWT[i] = (revSA[i] > 0) ? T[T.size() - revSA[i]] : T.front();
}

/**
 * @brief Write a std::string to a binary file
 * @param name The filename.
 * @param str The std::string to write.
 */
void writeStringToBinaryFile(const std::string& name, const std::string& str) {
    std::ofstream outFile(name, std::ios::binary);
    if (!outFile) {
        std::cerr << "Error opening file for writing: " << name << std::endl;
        return;
    }

    // First, write the size of the string
    length_t size = str.size();
    outFile.write(reinterpret_cast<const char*>(&size), sizeof(size));
    // Then, write the string data itself
    outFile.write(str.c_str(), size);
    outFile.close();
}

/**
 * @brief Preprocess the FASTA files.
 * @param fastaFiles The paths to the FASTA files.
 * @param baseFN The base filename.
 * @param T The text. (output)
 * @param seedLength The length of the seed for seeded replacement.
 * @param noWriting Whether to write the text to disk. Default is false.
 */
length_t preprocessFastaFiles(const std::vector<std::string>& fastaFiles,
                              const std::string& baseFN, string& T,
                              length_t seedLength, bool noWriting = false,
                              const std::string& taggingCategoriesFile = "") {

    // Seeded random number generation
    std::minstd_rand gen(42);

    // Generate random seed for seeded replacement
    std::string seed;
    size_t seedIndex = 0;
    for (length_t i = 0; i < seedLength; i++) {
        seed += replaceNonACGT('N', gen, seed, seedIndex);
    }

    // Reset gen to original seed
    gen.seed(42);

    std::function<char(char, size_t&)> replaceFunc;

    replaceFunc = [&gen, &seed](char c, size_t& seedIndex) -> char {
        return replaceNonACGT(c, gen, seed, seedIndex);
    };

    logger.logInfo("Using X replacement for non-ACGT characters...");

    std::vector<length_t> positions; // the start positions of each
    std::vector<string> seqNames;    // the name of each sequence
    std::vector<length_t> tags;      // the tag of each sequence
    std::vector<length_t> firstSeqIDPerFile;

    vector<string> taggingCategoriesVec;
    if (taggingCategoriesFile != "") {
        logger.logInfo("Reading tagging categories file: " +
                       taggingCategoriesFile);
        std::ifstream tagFile(taggingCategoriesFile);
        if (!tagFile.is_open()) {
            logger.logError("Error opening file: " + taggingCategoriesFile);
            throw std::runtime_error("Error opening file: " +
                                     taggingCategoriesFile);
        }
        std::string line;
        while (std::getline(tagFile, line)) {
            if (line.empty())
                continue; // Skip empty lines
            taggingCategoriesVec.push_back(line);
        }
        tagFile.close();

        if (taggingCategoriesVec.empty()) {
            logger.logError("Tagging categories file is empty");
            throw std::runtime_error("Tagging categories file is empty");
        }

        // Max 2**16 categories
        if (taggingCategoriesVec.size() > 65536) {
            logger.logError("Too many tagging categories. Maximum is 65536");
            throw std::runtime_error(
                "Too many tagging categories. Maximum is 65536");
        }

        logger.logInfo("Read " + std::to_string(taggingCategoriesVec.size()) +
                       " tagging categories");
    }

    for (const auto& fastaFile : fastaFiles) {
        logger.logInfo("Preprocessing FASTA file " + fastaFile);
        concatenateAndTransform(fastaFile, T, positions, seqNames, tags,
                                taggingCategoriesVec, firstSeqIDPerFile,
                                replaceFunc);
        checkTextSize(T.size(), false);
    }
    // add the end of the text
    positions.emplace_back(T.size());
    // Ensure the concatenation ends with a dollar sign
    if (T.back() != '$') {
        T += '$';
    }
    logger.logInfo("Read and concatenated " + std::to_string(seqNames.size()) +
                   " sequence(s) from " + std::to_string(fastaFiles.size()) +
                   " file(s)");
    checkTextSize(T.size());

    if (!noWriting) {
        logger.logInfo("Writing concatenated uppercase sequence to disk...");
        writeStringToBinaryFile(baseFN + ".txt.bin", T);
    }

    if (taggingCategoriesVec.size() > 0) {
        assert(tags.size() == seqNames.size());

        // Remove the last extension
        logger.logInfo("Writing tags to " + baseFN + ".reference.tags...");
        std::ofstream ofs(baseFN + ".reference.tags", ios::binary);
        ofs.write((char*)tags.data(), tags.size() * sizeof(length_t));
        ofs.close();
        logger.logInfo("Wrote tags to " + baseFN + ".reference.tags");

        logger.logInfo("Writing tagging categories to " + baseFN +
                       ".catagories.tags...");
        std::ofstream ofs2(baseFN + ".categories.tags");
        for (const auto& tag : taggingCategoriesVec) {
            ofs2 << tag << std::endl;
        }
        ofs2.close();
        logger.logInfo("Wrote tagging categories to " + baseFN +
                       ".categories.tags");
    }

    createAndWriteHeaderInfo(baseFN, seqNames, positions);
    writePositionsAndSequenceNames(baseFN, positions, seqNames,
                                   firstSeqIDPerFile);
    return seqNames.size();
}

/**
 * @brief Write the meta information to a file.
 * @param baseFN The base filename.
 */
void writeMetaInfo(const string& baseFN) {
    ofstream metaFile(baseFN + ".meta");

    // Write the build tag to a file
    metaFile << TAGGER_BUILD_INDEX_TAG << endl;
    // Write the 64 or 32-bit compiled info
    metaFile << sizeof(length_t) << endl;
    metaFile.close();
}

/**
 * @brief Generate the BWT of the text.
 * @param T The text.
 * @param SA The suffix array.
 * @param BWT The BWT. (output)
 */
void generateBWT(const string& T, const vector<length_t>& SA, string& BWT) {
    // build the BWT
    logger.logInfo("Generating BWT...");
    BWT.resize(T.size());
    for (size_t i = 0; i < SA.size(); i++)
        BWT[i] = (SA[i] > 0) ? T[SA[i] - 1] : T.back();
}

/**
 * @brief Write the character counts and create the alphabet.
 * @param baseFN The base filename.
 * @param T The text.
 * @param sigma The alphabet. (output)
 * @param charCounts The character counts.
 */
void writeCharCountsAndCreateAlphabet(const string& baseFN, const string& T,
                                      Alphabet<ALPHABET>& sigma,
                                      vector<length_t>& charCounts) {
    // count the frequency of each characters in T
    countChars(T, charCounts);
    // write the character counts table
    writeCharCounts(baseFN, charCounts);
    // Create the alphabet
    createAlphabet(T, charCounts, sigma);
}

/**
 * @brief Create the suffix array with a sanity check.
 * @param SA The suffix array. (output)
 * @param T The text.
 */
void createSAWithSanityCheck(vector<length_t>& SA, const string& T) {
    createSuffixArray(T, SA);

    // perform a sanity check on the suffix array
    sanityCheck(T, SA);
}

/**
 * @brief Create the suffix array of the reversed text with a sanity check.
 * @param revSA The suffix array of the reversed text. (output)
 * @param T The text.
 */
void createRevSAWithSanityCheck(vector<length_t>& revSA, string& T) {

    if (T.size() < (2ull << 32)) {
        std::string revT = T;
        std::reverse(revT.begin(), revT.end());

        // create the suffix array of the reverse text
        createSuffixArray(revT, revSA);
        revT.clear();
    } else {
        std::reverse(T.begin(), T.end());
        createSuffixArray(T, revSA);
        std::reverse(T.begin(), T.end());
    }

    // perform a sanity check on the suffix array
    sanityCheck(T, revSA);
}

/**
 * @brief Create the suffix array and BWT
 * @param T The text.
 * @param SA The suffix array. (output)
 * @param BWT The BWT. (output)
 */
void createSuffixArrayAndBWT(const string& T, vector<length_t>& SA,
                             string& BWT) {
    createSAWithSanityCheck(SA, T);
    generateBWT(T, SA, BWT);
}

// ============================================================================
// FUNCTIONALITY FOR RLC
// ============================================================================
/**
 * Get run index of a given position.
 * @param position The position for which to find the run index.
 * @param size The amount of input-output intervals.
 * @param dPair Vector with pairs of input-output interval start positions.
 * @returns The run index of the given position.
 */
length_t getRunIndex(const length_t position, const length_t size,
                     const MoveLFReprBP& rows) {
    length_t rightBoundary = size - 1;
    length_t leftBoundary = 0;
    // Iteratively make the possible range smaller by binary search, until only
    // 1 interval remains.
    while (rightBoundary - leftBoundary >= 1) {
        // Use the middle of the possible range as a test value.
        length_t testIndex = ((rightBoundary + leftBoundary) / 2) + 1;

        // Eliminate half of the possible range by comparing the value to the
        // test value.
        if (rows.getInputStartPos(testIndex) <= position) {
            leftBoundary = testIndex;
        } else {
            rightBoundary = testIndex - 1;
        }
    }

    assert(position >= rows.getInputStartPos(leftBoundary));
    assert(leftBoundary == size - 1 ||
           position < rows.getInputStartPos(leftBoundary + 1));

    return leftBoundary;
}

/**
 * Fill dIndex array.
 * @param tIn Balanced tree: stores pairs of I_out with increasing input
 * interval index.
 * @param rows The vector to fill, rows[i] = input-output interval i +
 * index j of input interval containing q_i.
 * @param size The amount of input-output intervals.
 * @param bwtSize The BWT size.
 */
void fillRows(const map<length_t, pair<uint8_t, length_t>>& tIn,
              MoveLFReprBP& rows, const length_t size, const length_t bwtSize) {
    rows.initialize(tIn.size(), bwtSize);

    length_t i = 0;
    for (auto it = tIn.begin(); it != tIn.end(); ++it) {
        rows.setRowValues(i, it->second.first, it->first, it->second.second, 0);
        assert(rows.getRunHead(i) == it->second.first);
        assert(rows.getInputStartPos(i) == it->first);
        assert(rows.getOutputStartPos(i) == it->second.second);
        i++;
    }

    i = 0;
    for (auto it = tIn.begin(); it != tIn.end(); ++it) {
        length_t outputRunIndex = getRunIndex(it->second.second, size, rows);
        rows.setOutputStartRun(i, outputRunIndex);
        assert(rows.getOutputStartRun(i) == outputRunIndex);
        i++;
    }

    rows.setRowValues(size, 0, bwtSize, bwtSize, size);
}

/**
 * Creates the Move structure.
 * @param baseFN the base filename to write the move structure to.
 * @param BWT The burrows-wheeler transform in string format (uncompressed).
 * @param charCounts Accumulated number of characters in lex order.
 * @param sigma The alphabet.
 */
length_t createAndWriteMove(const string& baseFN, const string& BWT,
                            const vector<length_t>& charCounts,
                            const Alphabet<ALPHABET>& sigma,
                            SparseBitvec& sparseEndOfRun) {

    // Get the alphabet size
    size_t S = sigma.size();

    // Set bwt size
    length_t bwtSize = BWT.size();

    // balanced tree: stores pairs of I_out with increasing input interval index
    map<length_t, pair<uint8_t, length_t>> tIn;

    // Create accumulated charCounts
    vector<length_t> charCountsAcc(S, 0);
    length_t total = 0;
    for (size_t i = 0; i < S; i++) {
        char c = sigma.i2c(i);
        charCountsAcc[i] = total;
        total += charCounts[c];
    }

    vector<bool> isEndOfRun(bwtSize + 1, false);

    // fill tIn, tOut
    vector<length_t> charsVisited(S, 0);
    length_t prevC = S;
    length_t zeroCharPos = bwtSize;
    for (length_t i = 0; i < bwtSize; i++) {

        char _c = BWT[i];
        length_t c = sigma.c2i(_c);
        assert(c < S);
        if (c == 0) {
            zeroCharPos = i;
        }
        if (prevC != c) {
            length_t lf = charCountsAcc[c] + charsVisited[c];
            tIn[i] = make_pair(c, lf);
            isEndOfRun[i == 0 ? bwtSize - 1 : i - 1] = true;
        }

        charsVisited[c]++;
        prevC = c;
    }

    sparseEndOfRun = SparseBitvec(isEndOfRun);
    isEndOfRun.clear();

    // Set dPair and dIndex size
    length_t arraySize = tIn.size();
    length_t size = arraySize;
    logger.logInfo("There are " + to_string(size) + " runs.");

    // Fill dPair and dIndex
    MoveLFReprBP rows;
    rows.setZeroCharPos(zeroCharPos);
    fillRows(tIn, rows, size, BWT.size());

    string moveFileName = baseFN + ".LFBP";
    rows.write(moveFileName);
    logger.logInfo("\tWrote file " + moveFileName);

    return size;
}

void writeIntVectorBinary(const std::string& filename,
                          const std::vector<uint16_t>& array) {
    // convert to int_vector
    uint8_t width =
        (uint8_t)ceil(log2(*max_element(array.begin(), array.end())));
    sdsl::int_vector<> intVector(array.size(), 0, width);
    for (size_t i = 0; i < array.size(); i++) {
        intVector[i] = array[i];
    }
    std::ofstream ofs(filename, std::ios::binary);
    intVector.serialize(ofs);
    ofs.close();
}

/**
 * Builds samplesFirst and samplesLast from the suffix array (SA) and BWT.
 * Fills samplesFirst and samplesLast at character changes in BWT for length_t.
 *
 * @param samplesFirst The samplesFirst array to fill.
 * @param samplesLast The samplesLast array to fill.
 * @param SA The suffix array.
 * @param BWT The Burrows-Wheeler transform.
 */
void buildSamples(vector<length_t>& samplesFirst, vector<length_t>& samplesLast,
                  const vector<length_t>& SA, const string& BWT) {

    samplesFirst.emplace_back(SA[0]);
    for (size_t pos = 0; pos < BWT.size() - 1; pos++) {
        if (BWT[pos] != BWT[pos + 1]) {
            samplesLast.emplace_back(SA[pos]);
            samplesFirst.emplace_back(SA[pos + 1]);
        }
    }
    samplesLast.emplace_back(SA[BWT.size() - 1]);
}

/**
 * Generate predecessor structures
 * @param movePhi MovePhi array
 * @param movePhiInv MovePhiInv array
 * @param [out] predFirst predecessor bitvector of the samplesFirst array
 * @param [out] predLast predecessor bitvector of the samplesLast array
 * @param textLength length of the text
 */
void generatePredecessors(const vector<length_t>& movePhiInv,
                          SparseBitvec& predLast, const length_t textLength) {

    vector<bool> predLastBV(textLength + 1, false);

    for (const length_t& row : movePhiInv) {
        predLastBV[row > 0 ? row - 1 : textLength - 1] = true;
    }
    predLast = SparseBitvec(predLastBV);
}

/**
 * Generate predecessor structures
 * @param movePhi MovePhi array
 * @param [out] predFirst predecessor bitvector of the samplesFirst array
 * @param textLength length of the text
 */
void generatePredecessors(const MovePhiReprBP& move, SparseBitvec& pred,
                          const length_t textLength) {

    vector<bool> predFirstBV(textLength + 1, false);

    for (length_t i = 0; i < move.size(); i++) {
        predFirstBV[move.getInputStartPos(i) > 0 ? move.getInputStartPos(i) - 1
                                                 : textLength - 1] = true;
    }

    pred = SparseBitvec(predFirstBV);
}

/**
 * Generate predToRun arrays
 * @param samplesFirst samplesFirst array
 * @param samplesLast samplesLast array
 * @param [out] firstToRun mapping between rank of ones in predFirst bitvector
 * and run indices
 * @param [out] lastToRun mapping between rank of ones in predLast bitvector and
 * run indices
 */
void generatePredToRun(const vector<length_t>& samplesLast,
                       vector<length_t>& lastToRun, length_t textLength) {

    lastToRun.resize(samplesLast.size());
    iota(lastToRun.begin(), lastToRun.end(), 0);
    sort(lastToRun.begin(), lastToRun.end(),
         [&samplesLast, &textLength](length_t a, length_t b) {
             return (samplesLast[a] > 0 ? samplesLast[a] - 1 : textLength - 1) <
                    (samplesLast[b] > 0 ? samplesLast[b] - 1 : textLength - 1);
         });
}

#ifdef BIG_BWT_USABLE

void readSuffixArrayFile(const std::string& baseFN,
                         const std::string& extension,
                         vector<length_t>& samples, vector<length_t>& toRun,
                         length_t size, length_t nrOfRuns, bool reverse,
                         bool makeToRun) {
    logger.logInfo("Reading suffix array samples from " + baseFN + extension +
                   "...");

    string fileName = baseFN + extension;

    // Open the file
    FILE* file = fopen(fileName.c_str(), "rb");

    samples = vector<length_t>(nrOfRuns, 0);
    if (makeToRun) {
        toRun = vector<length_t>(nrOfRuns, 0);
    }
    std::vector<std::pair<length_t, length_t>> pos_run_pairs(nrOfRuns);

    uint64_t* buf = new uint64_t[1];

    for (length_t i = 0; i < nrOfRuns; ++i) {
        // Read the first SABYTES bytes into buf from file
        if (fread(buf, SABYTES, 1, file) != 1)
            throw runtime_error((fileName + " read failed."));

        // Read the next SABYTES bytes into buf from file
        if (fread(buf, SABYTES, 1, file) != 1)
            throw runtime_error((fileName + " read failed."));

        // Calculate sa_val from buf: take the first byte, ensure it's within
        // range, and adjust as needed based on the data size
        uint64_t sa_val = buf[0] % (1UL << (8 * SABYTES));
        if (reverse) {
            // Only for the reverse suffix array, sa_val must be corrected.
            // For the reverse case, the value is too small since Big-BWT
            // puts the sentinel character at the end of the reverse text as
            // well.
            sa_val = (sa_val < size - 1) ? (sa_val + 1) : 0;
        }

        // Store sa_val in samples array at index i
        assert(sa_val >= 0 && sa_val < size);
        samples[i] = (length_t)sa_val;

        // Store {sa_val, i} pair in pos_run_pairs array
        pos_run_pairs[i] = {sa_val > 0 ? (length_t)sa_val - 1 : size - 1, i};
    }

    delete[] buf;

    if (!reverse) {

        logger.logInfo(
            "Creating the mapping between the predecessor bits and the "
            "runs...");

        std::sort(pos_run_pairs.begin(), pos_run_pairs.end());

        std::vector<length_t> positions;
        for (length_t i = 0; i < nrOfRuns; ++i) {
            positions.push_back(pos_run_pairs[i].first);
            if (makeToRun) {
                toRun[i] = pos_run_pairs[i].second;
            }
        }
    }

    fclose(file);
}

#endif

void performSubSampling(const std::vector<int>& maxLFValues,
                        const std::vector<length_t>& samplesLast,
                        const std::vector<uint16_t>& tagRunHeadsBWT,
                        const std::string& baseFN,
                        const std::vector<length_t>& startPos) {

    uint64_t r = samplesLast.size();
    assert(r == tagRunHeadsBWT.size());

    // Create a vector of indices (0, 1, 2, ..., r-1)
    std::vector<uint64_t> sorted_indices(r);
    for (uint64_t i = 0; i < r; i++) {
        sorted_indices[i] = i;
    }

    // Sort the indices based on the sample values
    std::sort(sorted_indices.begin(), sorted_indices.end(),
              [&samplesLast](const uint64_t& a, const uint64_t& b) -> bool {
                  return samplesLast[a] < samplesLast[b];
              });

    // Iterate over each maxLF value
    for (const auto& maxLFint : maxLFValues) {
        length_t maxLF = static_cast<length_t>(maxLFint);
        std::string filePrefix = baseFN + ".tag." + std::to_string(maxLF);

        // Create bitvector to mark sampled BWT runs
        std::vector<bool> subsampled_runIdx_bv(r, 0);

        uint64_t last_subsampled_value = samplesLast[sorted_indices[0]];
        subsampled_runIdx_bv[sorted_indices[0]] = 1;
        uint64_t subsampled_count = 1;

        // Track the next boundary to cross
        size_t nextBoundaryIdx = 1; // Start from the second boundary
        length_t nextBoundary = startPos[nextBoundaryIdx];

        // Iterate through sorted indices and select based on LF condition or
        // boundary
        for (uint64_t i = 1; i < r; i++) {
            uint64_t current_idx = sorted_indices[i];
            // Check if we crossed a boundary
            if (samplesLast[current_idx] >= nextBoundary) {
                // Force subsampling at the boundary
                subsampled_runIdx_bv[current_idx] = 1;
                last_subsampled_value = samplesLast[current_idx];
                subsampled_count++;

                // Move to the next boundary
                if (nextBoundaryIdx + 1 < startPos.size()) {
                    nextBoundary = startPos[++nextBoundaryIdx];
                }
            }
            // Check if distance exceeds maxLF
            else if (samplesLast[current_idx] - last_subsampled_value > maxLF) {
                subsampled_runIdx_bv[current_idx] = 1;
                last_subsampled_value = samplesLast[current_idx];
                subsampled_count++;
            }
        }

        std::vector<uint16_t> subsampled_values(subsampled_count);
        uint64_t j = 0;
        for (uint64_t i = 0; i < r; i++) {
            if (subsampled_runIdx_bv[i] == 1) {
                subsampled_values[j++] = tagRunHeadsBWT[i];
            }
        }

        assert(j == subsampled_count);

        SparseBitvec subsampled_runIdx_bv_sparse(subsampled_runIdx_bv);

        // Log and write outputs
        logger.logInfo("Maximum number of LF steps between two tag samples: " +
                       std::to_string(maxLF));
        logger.logInfo("Subsampled " + std::to_string(subsampled_count) +
                       " tags out of " + std::to_string(r) + " tags");

        subsampled_runIdx_bv_sparse.write(filePrefix + ".bv.bwt");

        writeIntVectorBinary(filePrefix + ".heads.bwt", subsampled_values);
    }
}

void createTaggingData(const string& baseFN, const length_t& numberOfReferences,
                       const length_t& bwtSize, const SparseBitvec& predLast,
                       const vector<length_t>& samplesFirst,
                       const vector<length_t>& samplesLast,
                       const vector<length_t>& lastToRun,
                       const vector<length_t>& startPos,
                       SparseBitvec& endOfBWTRun,
                       const vector<int>& maxLFValues) {
    // Try to open the tagging files
    ifstream ifs(baseFN + ".reference.tags");
    if (ifs.is_open()) {
        logger.logInfo("Creating tagging data...");

        vector<length_t> referenceTags;
        referenceTags.resize(numberOfReferences);
        ifs.read((char*)&referenceTags[0],
                 numberOfReferences * sizeof(length_t));
        ifs.close();

        vector<bool> tagRunEnd(bwtSize + 1, false);
        vector<uint16_t> tagRunHeads;
        tagRunHeads.reserve(bwtSize);

        vector<uint16_t> tagRunHeadsBWT;
        tagRunHeadsBWT.reserve(endOfBWTRun.rank(bwtSize));

        length_t currentPos = bwtSize - 1;

        auto it = upper_bound(startPos.begin(), startPos.end(), currentPos);
        length_t referenceIdx = distance(startPos.begin(), --it);

        // assert(0 <= referenceIdx &&
        //        referenceIdx < numberOfReferences);
        uint16_t currentTag = 0; // Sentinel character will never be accessed,
                                 // so tag does not matter here

        if (endOfBWTRun[0]) {
            tagRunHeadsBWT.push_back(currentTag);
        }

        for (length_t i = 1; i < bwtSize; i++) {

            if (bwtSize > 100 && i % (bwtSize / 100) == 0) {
                logger.logDeveloper("Progress: " + std::to_string(i) + " / " +
                                    std::to_string(bwtSize) + " (" +
                                    std::to_string((i * 100) / bwtSize) + "%)");
            }
            // Find the rank of the predecessor of pos in a circular manner.
            length_t predRank = predLast.predecessorRankCircular(currentPos);
            // Select the predecessor position using its rank.
            length_t pred = predLast.select(predRank);

            // Calculate the distance from the predecessor to the current
            // position.
            length_t delta =
                pred < currentPos ? currentPos - pred : currentPos + 1;

            assert(predRank < lastToRun.size());

            // Ensure that phiInverse(SA[n-1]) is not called; the predecessor
            // rank must be valid.
            assert(lastToRun[predRank] < samplesFirst.size() - 1);

            // Get the next sample from the samplesFirst array.
            length_t prev_sample =
                samplesFirst[lastToRun[predRank] + 1] > 0
                    ? samplesFirst[lastToRun[predRank] + 1] - 1
                    : bwtSize - 1;

            // Calculate and return the new position, modulo the text length.
            length_t nextPos = (prev_sample + delta) % bwtSize;

            // Get next tag
            it = upper_bound(startPos.begin(), startPos.end(), nextPos);
            referenceIdx = distance(startPos.begin(), --it);
            // assert(0 <= referenceIdx &&
            //        referenceIdx < numberOfReferences);
            uint16_t nextTag = referenceTags[referenceIdx];

            if (currentTag != nextTag) {
                tagRunEnd[i - 1] = true;
                tagRunHeads.push_back(currentTag);
                currentTag = nextTag;
            }

            if (endOfBWTRun[i]) {
                tagRunHeadsBWT.push_back(currentTag);
            }

            currentPos = nextPos;
        }

        tagRunEnd[bwtSize - 1] = true;
        tagRunHeads.push_back(currentTag);

        // Report the number of tag runs
        logger.logInfo("Number of tag runs: " + to_string(tagRunHeads.size()));

        SparseBitvec tagRunEndBV(tagRunEnd);
        assert(tagRunEndBV.rank(bwtSize) == tagRunHeads.size());

        SparseBitvec sampledBWTEndRunTagsBV;
        vector<uint16_t> sampledBWTEndRunTags;
        performSubSampling(maxLFValues, samplesLast, tagRunHeadsBWT, baseFN,
                           startPos);

        assert(tagRunHeadsBWT.size() == endOfBWTRun.rank(bwtSize));

    } else {
        logger.logDeveloper(
            "No tagging data found. Skipping tagging data creation.");
    }
}

#ifdef BIG_BWT_USABLE

// Optimized helper function to replace multiple characters and log their counts
void replaceSentinel(std::string& str) {
    logger.logInfo("Replacing special characters in the BWT...");

    size_t countNull = 0; // Count for '\0'
    size_t countOne = 0;  // Count for '\1'
    size_t countTwo = 0;  // Count for '\2'

    // Traverse the string once to count and replace
    for (auto& ch : str) {
        if (ch == '\0') {
            countNull++; // Count '\0'
            ch = '$';    // Replace with '$'
        } else if (ch == '\1') {
            countOne++; // Count '\1'
            ch = '$';   // Replace with '$'
        } else if (ch == '\2') {
            countTwo++; // Count '\2'
            ch = '$';   // Replace with '$'
        }
    }

    // Log results after the single pass
    logger.logDeveloper("Found " + std::to_string(countNull) +
                        " \\0 characters.");
    logger.logDeveloper("Found " + std::to_string(countOne) +
                        " \\1 characters.");
    logger.logDeveloper("Found " + std::to_string(countTwo) +
                        " \\2 characters.");

    size_t totalCount = countNull + countOne + countTwo;

    // Error check
    if (totalCount != 1) {
        logger.logError("Found " + std::to_string(totalCount) +
                        " special characters in the BWT. Expected 1.");
        exit(1);
    }
}

void processSamplesAndPreds(const string& baseFN,
                            const vector<length_t>& samplesFirst,
                            const vector<length_t>& samplesLast,
                            const string& BWT, length_t numberOfReferences,
                            const vector<length_t>& startPos,
                            SparseBitvec& endOfRun,
                            const vector<int>& maxLFValues) {
    vector<length_t> lastToRun;

    // Generate the predToRun arrays
    logger.logInfo(
        "Mapping the predecessor bits to their corresponding runs...");
    generatePredToRun(samplesLast, lastToRun, BWT.size());

    SparseBitvec predFirst;
    SparseBitvec predLast;

    logger.logInfo("Generating the predecessor bitvector for the samples...");
    generatePredecessors(samplesLast, predLast, BWT.size());

    createTaggingData(baseFN, numberOfReferences, BWT.size(), predLast,
                      samplesFirst, samplesLast, lastToRun, startPos, endOfRun,
                      maxLFValues);

    lastToRun.clear();
}

void createIndex(string& T, const BuildParameters& params,
                 const Alphabet<ALPHABET>& sigma,
                 const vector<length_t>& charCounts,
                 const length_t& numberOfReferences,
                 const vector<length_t>& startPos) {
    // aliasing
    const auto& baseFN = params.baseFN;

    {
        // create SA and BWT
        vector<length_t> SA;
        string BWT;
        createSuffixArrayAndBWT(T, SA, BWT);
        // Create samplesLast and samplesFirst
        logger.logInfo("Sampling the suffix array values at run boundaries...");
        vector<length_t> samplesFirst;
        vector<length_t> samplesLast;
        buildSamples(samplesFirst, samplesLast, SA, BWT);

        // Clear the suffix array
        SA.clear();

        SparseBitvec endOfRun;

        // Create the Move structure
        logger.logInfo("Creating the move table...");
        createAndWriteMove(baseFN, BWT, charCounts, sigma, endOfRun);

        processSamplesAndPreds(baseFN, samplesFirst, samplesLast, BWT,
                               numberOfReferences, startPos, endOfRun,
                               params.maxLFValues);

        BWT.clear();
    }

    { // read or create the reverse suffix array
        vector<length_t> revSA;
        createRevSAWithSanityCheck(revSA, T);

        // build the reverse BWT
        string revBWT;
        createRevBWT(baseFN, T, revSA, sigma, revBWT);

        SparseBitvec endOfRun;

        // Create the Move structure
        logger.logInfo("Creating the reverse move table...");
        createAndWriteMove(baseFN + ".rev", revBWT, charCounts, sigma,
                           endOfRun);

        // Clear the reverse BWT
        revBWT.clear();
    }
}

void createIndexPFP(const string& baseFN, const length_t& numberOfReferences,
                    const vector<length_t>& startPos,
                    const vector<int>& maxLFValues) {

    // build the BWT
    string BWT;
    // Read the BWT from disk
    logger.logInfo("Reading " + baseFN + ".bwt...");
    readText(baseFN + ".bwt", BWT);

    replaceSentinel(BWT);

    length_t bwtSize = BWT.size();

    // count the frequency of each characters in T
    vector<length_t> charCounts;
    countChars(BWT, charCounts);

    // write the character counts table
    writeCharCounts(baseFN, charCounts);

    // Create the alphabet
    Alphabet<ALPHABET> sigma;
    createAlphabet(BWT, charCounts, sigma);

    { // Create the Move structure

        SparseBitvec endOfRun;

        logger.logInfo("Creating the move table...");
        length_t nrOfRuns =
            createAndWriteMove(baseFN, BWT, charCounts, sigma, endOfRun);

        BWT.clear();

        vector<length_t> samplesFirst;
        vector<length_t> firstToRun;
        readSuffixArrayFile(baseFN, ".ssa", samplesFirst, firstToRun, bwtSize,
                            nrOfRuns, false, false);

        vector<length_t> samplesLast;
        vector<length_t> lastToRun;
        readSuffixArrayFile(baseFN, ".esa", samplesLast, lastToRun, bwtSize,
                            nrOfRuns, false, true);

        { // Original phi operation
            SparseBitvec predFirst;
            SparseBitvec predLast;

            // Generate the predecessor bitvectors
            logger.logInfo(
                "Generating the predecessor bitvectors for the samples...");
            generatePredecessors(samplesLast, predLast, bwtSize);

            createTaggingData(baseFN, numberOfReferences, bwtSize, predLast,
                              samplesFirst, samplesLast, lastToRun, startPos,
                              endOfRun, maxLFValues);

            firstToRun.clear();
        }
    }

    logger.logInfo("Switching to reversed text...");

    {
        // build the BWT
        string revBWT;
        // Read the BWT from disk
        logger.logInfo("Reading " + baseFN + ".rev.bwt...");
        readText(baseFN + ".rev.bwt", revBWT);

        replaceSentinel(revBWT);

        SparseBitvec endOfRun;

        // Create the Move structure
        logger.logInfo("Creating the move table...");
        createAndWriteMove(baseFN + ".rev", revBWT, charCounts, sigma,
                           endOfRun);

        // Clear the reverse BWT
        revBWT.clear();
    }

    writeMetaInfo(baseFN);
}

int indexConstructingAfterPFP(const BuildParameters& params) {
    logger.logInfo(
        "Starting index construction after prefix-free parsing step...");

    std::array<string, 4> requiredExtensionsPFP = {".bwt", ".esa", ".ssa",
                                                   ".rev.bwt"};

    for (auto& ext : requiredExtensionsPFP) {
        string filename = params.baseFN + ext;

        // check if file with filename exists
        ifstream ifs(filename);
        if (!ifs) {
            logger.logError("Missing file: " + filename);
            logger.logError("Please run the prefix-free parsing step first.");
            return EXIT_FAILURE;
        }
    }

    vector<length_t> startPos;
    startPos.reserve(100000);

    ifstream ifs(params.baseFN + ".pos");
    // Read in without knowing the size beforehand
    while (true) {
        length_t pos;
        ifs.read((char*)&pos, sizeof(length_t));
        if (ifs.eof()) {
            break;
        }
        startPos.push_back(pos);
    }
    ifs.close();

    length_t numberOfReferences = startPos.size() - 1;

    createIndexPFP(params.baseFN, numberOfReferences, startPos,
                   params.maxLFValues);

    logger.logInfo("Index construction completed successfully!");
    logger.logInfo("Exiting... bye!");

    return EXIT_SUCCESS;
}

int preprocessingOnly(const BuildParameters& params) {
    logger.logInfo("Preprocessing only...");

    std::string T; // the (concatenated) text
    preprocessFastaFiles(params.fastaFiles, params.baseFN, T, params.seedLength,
                         true, params.taggingCategories);

    // Remove the sentinel character
    T.pop_back();

    // Write the preprocessed text to disk
    logger.logInfo("Writing concatenated uppercase sequence to disk...");
    std::ofstream ofs(params.baseFN);
    ofs << T;
    ofs.close();

    // Reverse the text
    logger.logInfo("Reversing text...");
    std::reverse(T.begin(), T.end());

    // Write the reversed text to disk
    logger.logInfo("Writing reversed text to disk...");
    ofs = std::ofstream(params.baseFN + ".rev");
    ofs << T;
    ofs.close();

    logger.logInfo("Fasta input preprocessed successfully!");
    logger.logInfo("Exiting... bye!");

    return EXIT_SUCCESS;
}
#endif
// ============================================================================
// INPUT FUNCTIONALITY FOR VANILLA AND RLC
// ============================================================================

/**
 * Process the input fasta files and create the index.
 * @param params the build parameters
 *
 */
void processFastaFiles(const BuildParameters& params) {
    std::string T; // the (concatenated) text

    bool noWriting = true;

    length_t numberOfReferences = preprocessFastaFiles(
        params.fastaFiles, params.baseFN, T, params.seedLength, noWriting,
        params.taggingCategories);

    // read the start positions
    vector<length_t> startPos;
    ifstream ifs(params.baseFN + ".pos");
    if (!ifs.is_open()) {
        logger.logError("Could not open file " + params.baseFN + ".pos");
        throw runtime_error("Could not open file " + params.baseFN + ".pos");
    }
    startPos.resize(numberOfReferences + 1);
    ifs.read((char*)&startPos[0], (numberOfReferences + 1) * sizeof(length_t));
    ifs.close();

    // create alphabet and charCounts
    Alphabet<ALPHABET> sigma;
    vector<length_t> charCounts;
    writeCharCountsAndCreateAlphabet(params.baseFN, T, sigma, charCounts);
    createIndex(T, params, sigma, charCounts, numberOfReferences, startPos);
    // write meta info
    writeMetaInfo(params.baseFN);
}

int main(int argc, char* argv[]) {

    logger.logInfo("Welcome to tagger's index construction!");

    BuildParameters params;
    if (!BuildParameters::parse(argc, argv, params)) {
        BuildParameters::showUsage();
        return EXIT_FAILURE;
    }

    logger.logInfo("Alphabet size is " + std::to_string(ALPHABET - 1) + " + 1");

#ifdef BIG_BWT_USABLE
    if (params.preprocessOnly) {
        return preprocessingOnly(params);
    }

    if (params.pfp) {
        return indexConstructingAfterPFP(params);
    }
#endif

    try {
        for (auto& fastaFile : params.fastaFiles) {
            if (!BuildParameters::validFastaExtension(fastaFile)) {
                throw runtime_error("Invalid fasta file extension for file " +
                                    fastaFile);
            }
        }
        processFastaFiles(params);
    } catch (const std::exception& e) {
        logger.logError("Fatal: " + string(e.what()));

        return EXIT_FAILURE;
    }

    logger.logInfo("Index construction completed successfully!");
    logger.logInfo("Exiting... bye!");

    return EXIT_SUCCESS;
}