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

#include "definitions.h" // for length_t, Orientation, FF, FR, RF, ALL
#include "fastq.h"       // for OutputChunk, OutputRecord, Reader, Out...
#include "indexinterface.h"
#include "bmove/bmove.h" // for BMove
#include "indexhelpers.h"               // for TextOcc, Counters
#include "logger.h"                     // for Logger, logger
#include "parameters/alignparameters.h" // for Parameters
#include "reads.h"                      // for ReadPair, Read, ReadBundle
#include "searchstrategy.h"             // for SearchStrategy

#include <algorithm>  // for max, copy, for_each, count, sort
#include <assert.h>   // for assert
#include <chrono>     // for duration, operator-, high_resolution_c...
#include <cmath>      // for sqrt, abs
#include <exception>  // for exception
#include <functional> // for ref, mem_fn, _Mem_fn
#include <iterator>   // for move_iterator, make_move_iterator
#include <memory>     // for allocator, allocator_traits<>::value_type
#include <mutex>      // for mutex, lock_guard, unique_lock
#include <numeric>    // for accumulate
#include <ostream>    // for operator<<, basic_ostream, stringstream
#include <stdint.h>   // for int32_t
#include <stdlib.h>   // for abs, EXIT_SUCCESS, EXIT_FAILURE
#include <string>     // for char_traits, operator<<, operator+
#include <thread>     // for thread
#include <utility>    // for move
#include <vector>     // for vector, _Bit_iterator, _Bit_reference

using namespace std;

//----------------------------------------------------------------------------
// Single-ended processing
//----------------------------------------------------------------------------

/**
 * Process a chunk of reads singled-ended using the given search strategy
 * @param input The input chunk of reads
 * @param output The output chunk (results will be added to this chunk)
 * @param strategy The search strategy to use
 * @param maxEDOrIdentity The maximum edit distance or minimum identity of the
 * matches
 */
void processChunk(vector<ReadBundle>& input, OutputChunk& output,
                  std::unique_ptr<SearchStrategy>& strategy,
                  const length_t& minMEMLength) {
    output.clear();
    output.reserve(input.size());
    for (auto& readBundle : input) {
        output.getMutableCounters().inc(Counters::TOTAL_READ_LENGTHS,
                                        readBundle.size());
        vector<TextOcc> matches = strategy->taggingFunction(
            readBundle, output.getMutableCounters(), minMEMLength);
        output.addSingleEndRecord(std::move(matches));
    }
}

/**
 * Entry for a worker thread for processing reads single-ended.
 * @param myReader The reader to read the input from
 * @param myWriter The writer to write the output to
 * @param strategy The search strategy to use
 * @param maxEDOrIdentity The maximum edit distance or minimum identity of
 * the matches
 */
void threadEntrySingleEnd(Reader& myReader, OutputWriter& myWriter,
                          std::unique_ptr<SearchStrategy>& strategy,
                          length_t& maxEDOrIdentity) {
    // local storage of reads
    vector<ReadBundle> input;
    OutputChunk output;

    size_t chunkID;
    while (myReader.getNextChunk(input, chunkID)) {

        auto start = chrono::high_resolution_clock::now();

        processChunk(input, output, strategy, maxEDOrIdentity);
        // mark output as ready
        output.set(chunkID);
        // add output to writer
        myWriter.commitChunk(output);

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end - start;
        myReader.addProcessingTime(elapsed.count());
    }

    myWriter.sendTermination(chunkID);
}

//----------------------------------------------------------------------------
// Logging
//----------------------------------------------------------------------------

/**
 * Log the benchmarking parameters
 * @param params The Parameters object
 * @param index The index
 * @param strategy The search strategy
 */
void logAlignmentParameters(const Parameters& params,
                            const IndexInterface& index,
                            const std::unique_ptr<SearchStrategy>& strategy) {
    stringstream ss;
    ss << "Finding maximal exact matches...\n";
    ss << "\tMinimum MEM length: " << params.minMEMLength << "\n";
#ifdef TAG_ARRAY_SUBSAMPLING
    ss << "\tMaximum number of LF steps: " << params.maxLFSteps << "\n";
#endif
    ss << "\tNumber of threads: " << params.nThreads << "\n";
    ss << "\tReads file: " << params.firstReadsFile << "\n";
    ss << "\tOutput file: " << params.outputFile << "\n";

    logger.logInfo(ss);
}

/**
 * @brief Processes command-line arguments.
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 * @param params Reference to Parameters object to store parsed arguments.
 * @return true if arguments are processed successfully.
 * @return false if there is an error in processing arguments (error will be
 * logged).
 */
bool processArguments(int argc, char* argv[], Parameters& params) {
    try {
        params = Parameters::processOptionalArguments(argc, argv);
    } catch (const std::exception& e) {
        logger.logError("Problem with processing parameters: " +
                        std::string(e.what()));
        Parameters::printHelp();
        return false;
    }

    if (!params.logFile.empty()) {
        logger.setLogFile(params.logFile);
    }

    return true;
}

namespace std {
/**
 * @brief Provides a `make_unique` function for C++11.
 *
 * @tparam T The type of the object to create.
 * @tparam Args The types of the arguments to pass to the constructor.
 * @param args The arguments to pass to the constructor.
 * @return std::unique_ptr<T> A unique pointer to the created object.
 */
template <typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}
} // namespace std

/**
 * @brief Creates a BMove object.
 *
 * @param params Parameters object containing necessary parameters.
 * @return std::unique_ptr<BMove> Pointer to the created BMove object, or
 * nullptr on failure (error will be logged).
 */
std::unique_ptr<BMove> createBMove(const Parameters& params) {
    try {
        return std::make_unique<BMove>(params.base, true, params.kmerSize
#ifdef TAG_ARRAY_SUBSAMPLING
                                       ,
                                       params.maxLFSteps
#endif
        );
    } catch (const std::exception& e) {
        logger.logError("Problem with index loading: " + std::string(e.what()));
        return nullptr;
    }
}

/**
 * @brief Creates a SearchStrategy object.
 *
 * @param params Parameters object containing necessary parameters.
 * @param index Reference to IndexInterface object.
 * @return std::unique_ptr<SearchStrategy> Pointer to the created SearchStrategy
 * object, or nullptr on failure (error will be logged).
 */
std::unique_ptr<SearchStrategy> createSearchStrategy(const Parameters& params,
                                                     IndexInterface& index) {
    try {
        return params.createStrategy(index);
    } catch (const std::exception& e) {
        logger.logError("Problem with search scheme loading: " +
                        std::string(e.what()));
        return nullptr;
    }
}

/**
 * @brief Creates a Reader object.
 *
 * @param params Parameters object containing necessary parameters.
 * @return std::unique_ptr<Reader> Pointer to the created Reader object, or
 * nullptr on failure (error will be logged).
 */
std::unique_ptr<Reader> createReader(const Parameters& params) {
    try {
        auto reader = std::make_unique<Reader>(params.firstReadsFile, "");
        size_t targetChunk = 64 * 9;
        reader->startReaderThread(targetChunk, params.nThreads);
        return reader;
    } catch (const std::exception& e) {
        logger.logError("Problem with reading the reads file(s): " +
                        std::string(e.what()));
        return nullptr;
    }
}

/**
 * @brief Creates an OutputWriter object.
 *
 * @param params Parameters object containing necessary parameters.
 * @return std::unique_ptr<OutputWriter> Pointer to the created OutputWriter
 * object, or nullptr on failure (error will be logged).
 */
std::unique_ptr<OutputWriter> createOutputWriter(const Parameters& params) {
    try {
        auto writer = std::make_unique<OutputWriter>(
            params.outputFile, params.base + ".headerSN.bin", params.command,
            SINGLE_END, params.reorder);
        writer->start(params.nThreads * 500, params.nThreads);
        return writer;
    } catch (const std::exception& e) {
        logger.logError("Problem with writing to the output file: " +
                        std::string(e.what()));
        return nullptr;
    }
}

/**
 * @brief Starts worker threads for single-end mode and waits on them
 *
 * @param workers Vector of worker threads.
 * @param reader Reference to Reader object.
 * @param writer Reference to OutputWriter object.
 * @param strategy Pointer to SearchStrategy object.
 * @param minMEMLength Reference to length_t for max or identity parameter.
 */
void startSingleEndWorkers(std::vector<std::thread>& workers, Reader& reader,
                           OutputWriter& writer,
                           std::unique_ptr<SearchStrategy>& strategy,
                           length_t& minMEMLength) {
    for (size_t i = 0; i < workers.size(); i++) {
        workers[i] = std::thread(threadEntrySingleEnd, std::ref(reader),
                                 std::ref(writer), std::ref(strategy),
                                 std::ref(minMEMLength));
    }
    for (auto& worker : workers) {
        worker.join();
    }
}

/**
 * @brief Starts worker threads for processing and waits on them as well as the
 * reader and writer to finish.
 *
 * @param params Parameters object containing necessary parameters.
 * @param reader Reference to Reader object.
 * @param writer Reference to OutputWriter object.
 * @param strategy Pointer to SearchStrategy object.
 * @return true if worker threads start and process successfully.
 * @return false if there is an error during processing (error will be logged).
 */
bool startWorkerThreads(Parameters& params, Reader& reader,
                        OutputWriter& writer,
                        std::unique_ptr<SearchStrategy>& strategy) {
    try {
        auto start = std::chrono::high_resolution_clock::now();

        std::vector<std::thread> workers(params.nThreads);
        length_t minMEMLength = params.minMEMLength;

        startSingleEndWorkers(workers, reader, writer, strategy, minMEMLength);

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;

        logger.logInfo(std::to_string(params.nThreads) +
                       " mapping threads finished mapping in " +
                       std::to_string(elapsed.count()) + "s");

        reader.joinReaderThread();
        writer.joinWriterThread();
    } catch (const std::exception& e) {
        logger.logError("Problem during mapping: " + std::string(e.what()));
        return false;
    }
    return true;
}

/**
 * @brief Cleans up resources and exits the program.
 *
 * @param indexPtr Pointer to the index interface object.
 * @param readerPtr Pointer to Reader object.
 * @param writerPtr Pointer to OutputWriter object.
 */
void cleanupAndExit(std::unique_ptr<IndexInterface>& indexPtr,
                    std::unique_ptr<Reader>& readerPtr,
                    std::unique_ptr<OutputWriter>& writerPtr) {
    indexPtr.reset();
    readerPtr.reset();
    writerPtr.reset();
    logger.logInfo("Bye...\n");
}

/**
 * @brief Entry point of the program.
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 * @return int Exit status.
 */
int main(int argc, char* argv[]) {
    logger.logDeveloper("Developer mode is on");
    Parameters params;
    if (!processArguments(argc, argv, params)) {
        return EXIT_FAILURE;
    }

    std::unique_ptr<IndexInterface> indexPtr;

    // Conditional compilation
    indexPtr = createBMove(params);

    if (!indexPtr) {
        return EXIT_FAILURE;
    }

    std::unique_ptr<SearchStrategy> strategy =
        createSearchStrategy(params, *indexPtr);
    if (!strategy) {
        return EXIT_FAILURE;
    }

    std::unique_ptr<Reader> readerPtr = createReader(params);
    if (!readerPtr) {
        return EXIT_FAILURE;
    }

    std::unique_ptr<OutputWriter> writerPtr = createOutputWriter(params);
    if (!writerPtr) {
        return EXIT_FAILURE;
    }

    logAlignmentParameters(params, *indexPtr, strategy);

    if (!startWorkerThreads(params, *readerPtr, *writerPtr, strategy)) {
        return EXIT_FAILURE;
    }

    cleanupAndExit(indexPtr, readerPtr, writerPtr);

    return EXIT_SUCCESS;
}