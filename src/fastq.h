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

#ifndef FASTQ_H
#define FASTQ_H

#include "definitions.h"  // for length_t
#include "indexhelpers.h" // for TextOcc, Counters
#include "reads.h"        // for Read
#include "seqfile.h"      // for FileType

#include <algorithm>          // for max, transform
#include <chrono>             // for minutes, seconds, steady_clock, time_p...
#include <condition_variable> // for condition_variable
#include <ctype.h>            // for toupper
#include <deque>              // for deque
#include <functional>         // for bind, function, _1
#include <memory>             // for allocator_traits<>::value_type
#include <mutex>              // for mutex, lock_guard
#include <stddef.h>           // for size_t
#include <string>             // for string, basic_string
#include <thread>             // for thread
#include <utility>            // for move
#include <vector>             // for vector

#include <parallel_hashmap/btree.h>

// ============================================================================
// SEQUENCE RECORD
// ============================================================================

/**
 * A sequence record object contains a single read that is read from file.
 */
class SequenceRecord : public ReadBundle {
  private:
    /**
     * Read record from read file
     * @param Opened read file
     * @return True upon successful read, false otherwise
     */
    bool readFromFileFASTA(SeqFile& input);
    /**
     * Read record from read file
     * @param Opened read file
     * @return True upon successful read, false otherwise
     */
    bool readFromFileFASTQ(SeqFile& input);

    std::function<bool(SeqFile&)> readFunc;

  public:
    bool readFromFile(SeqFile& input) {
        if (readFunc(input)) {
            std::transform(read.begin(), read.end(), read.begin(), ::toupper);
            return true;
        }
        return false;
    }

    SequenceRecord(bool fastq) : ReadBundle() {
        readFunc = fastq ? std::bind(&SequenceRecord::readFromFileFASTQ, this,
                                     std::placeholders::_1)
                         : std::bind(&SequenceRecord::readFromFileFASTA, this,
                                     std::placeholders::_1);
        // convert read to upper case
        std::transform(read.begin(), read.end(), read.begin(), ::toupper);
    }
};

// ============================================================================
// READ BLOCK CLASS
// ============================================================================

/**
 * A read block object contains a number of single-end reads
 * that are read from file(s). In general, a readBlock consists of
 * multiple read chunks.
 */
class ReadBlock : public std::vector<SequenceRecord> {

  private:
    size_t nextChunkOffset; // offset of the next record to be read
    size_t targetChunkSize; // the used targetChunkSize for this block

  public:
    /**
     * Construct read block from a vector of FastQ records (deep copy)
     * @param buffer A vector of FastQ records
     */
    ReadBlock(const std::vector<SequenceRecord>& buffer)
        : std::vector<SequenceRecord>(buffer), nextChunkOffset(0) {
    }

    /**
     * Default, deleted copy and default move constructor
     */
    ReadBlock() : nextChunkOffset(0) {
    }
    ReadBlock(const ReadBlock&) = delete;
    ReadBlock(ReadBlock&&) = default;

    /**
     * Deleted copy and default move assignment operator
     */
    ReadBlock& operator=(const ReadBlock&) = delete;
    ReadBlock& operator=(ReadBlock&& rhs) = default;

    /**
     * Return true of at least one chunk is available for reading
     * @return true or false
     */
    bool chunkAvailable() const {
        return nextChunkOffset < this->size();
    }

    /**
     * Get next available input record chunk
     * @param buffer Record buffer to write to (contents will be appended)
     */
    void getNextChunk(std::vector<ReadBundle>& buffer);

    /**
     * Read a block from file (single-end reads)
     * @param file1 Fastq file
     * @param targetBlockSize Desired number of nucleotides in this block
     * @param fastq True if the file is a fastq file, false if it is a fasta
     * file
     */
    void readFromFile(SeqFile& file1, size_t targetBlockSize, bool fastq);

    /**
     * Return the target chunk size that was used in this block
     * @return The target chunk size
     */
    size_t getTargetChunkSize() const {
        return targetChunkSize;
    }

    /**
     * Set the target chunk size that was used in this block
     * @param targetChunkSize The target chunk size
     */
    void setTargetChunkSize(size_t targetChunkSize) {
        this->targetChunkSize = targetChunkSize;
    }
};

// ============================================================================
// FASTQ READER
// ============================================================================

/**
 * @class Reader
 * @brief Class for reading FASTQ or FASTA files.
 *
 * The Reader class provides functionality for reading FASTQ/FASTA files,
 * single-end reads. It supports multithreaded reading,
 * allowing for efficient processing of large FASTQ/FASTA datasets.
 */
class Reader {
  private:
    std::string filename1;     // name of the input file /1
    FileType fileType1;        // file type (FASTQ, FASTQ.GZ, FASTA, FASTA.GZ)
    std::string baseFilename1; // file name w/o extension /1

    bool fastq1 = true; // true if first file is fastq, false otherwise
    bool fastq2 = true; // true if second file is fastq, false otherwise

    const int numReadBlocks = 2; // number of read blocks

    size_t targetChunkSize; // target size for a single chunk
    size_t targetBlockSize; // target size for a single block

    size_t numberOfWorkerThreads = 0; // number of worker threads

    std::thread iThread; // input thread

    // input thread variables (protect by inputMutex)
    std::mutex inputMutex;              // input thread mutex
    std::condition_variable inputReady; // input blocks ready condition
    std::vector<ReadBlock> inputBlocks; // input blocks that are empty

    // worker thread variables (protect by workMutex)
    std::mutex workMutex;              // worker thread mutex
    std::condition_variable workReady; // work ready condition
    std::deque<ReadBlock> workBlocks;  // blocks being processed
    size_t chunkID;                    // unique ID per chunk

    std::mutex processMutex; // mutex for processing time
    double lowEndProcessingTime = 0.3;
    double highEndProcessingTime = 0.5;
    double midProcessingTime =
        (highEndProcessingTime + lowEndProcessingTime) / 2;
    std::vector<double> processingTimes;

    /**
     * Entry routine for the input thread.
     * This function is responsible for reading the input files and
     * populating the input blocks.
     */
    void readerThread();

  public:
    /**
     * Default constructor.
     * @param filename1 Name of the input file /1.
     */
    Reader(const std::string& filename1);

    /**
     * Return the filename without the extension /1.
     * @return The filename without the extension.
     */
    std::string getBaseFilename1() const {
        return baseFilename1;
    }

    /**
     * Start the reader thread.
     * @param targetChunkSize Target size for a single chunk.
     * @param numWorkerThreads Number of worker threads.
     */
    void startReaderThread(size_t targetChunkSize, size_t numWorkerThreads);

    /**
     * Get the next read chunk from the input.
     * @param chunk Buffer in which to store the records (output).
     * @param chunkID Unique chunk identifier (output).
     * @return False if no more reads are available, true otherwise.
     */
    bool getNextChunk(std::vector<ReadBundle>& chunk, size_t& chunkID);

    /**
     * Join the reader thread.
     * This function blocks until the reader thread has finished.
     */
    void joinReaderThread();

    /**
     * Add processing time of a single chunk
     * @param time Processing time of a single chunk
     */
    void addProcessingTime(double time) {
        std::lock_guard<std::mutex> lock(processMutex);
        processingTimes.push_back(time);
    }
};

// ============================================================================
// OUTPUT RECORD
// ============================================================================

/**
 * An output record object contains the non-redundant occurrences for a single
 * read. 
 */
class OutputRecord {
  public:
    std::shared_ptr<std::vector<TextOcc>>
        outputOcc; // The non-redundant occurrences for this record

    /**
     * Constructor with single-end occurrences
     */
    OutputRecord(const std::vector<TextOcc>& occurrences)
        : outputOcc(std::make_shared<std::vector<TextOcc>>(occurrences)) {
    }

    // Move constructor
    OutputRecord(OutputRecord&& other) noexcept
        : outputOcc(std::move(other.outputOcc)) {};

    // Move assignment operator
    OutputRecord& operator=(OutputRecord&& other) noexcept {
        if (this != &other) {
            outputOcc = std::move(other.outputOcc);
        }
        return *this;
    }

    // Destructor (if needed)
    ~OutputRecord() = default;
};

// ============================================================================
// OUTPUT CHUNK
// ============================================================================

/**
 * An output chunk object contains a number of output records that need to be
 * written to disk. As well as the counters that were gathered during the
 * processing of the records. A special end output chunk is used to signal
 * termination.
 */
class OutputChunk {
  private:
    bool empty = true; // If the chunk is filled in
    size_t chunkID;    // The unique id of this chunk
    bool end = false;  // is this an end signal chunk
    std::vector<OutputRecord>
        records; // the records for the input of this chunk

    Counters counters; // the performance counters for this chunk

    /**
     * Set this as an end outputChunk and mark the chunk as non-empty
     */
    void setEnd(size_t id) {
        end = true;
        set(id);
    }

  public:
    /**
     * @returns, true if the chunk has not been filled in yet (= worker has
     * not finished yet), false if the chunk is filled in
     */
    bool isEmpty() const {
        return empty;
    }

    /**
     * Get the id
     * @returns the chunk id
     */
    size_t getChunkID() const {
        return chunkID;
    }
    /**
     * Default constructor. Creates an empty output chunk with no records
     */
    OutputChunk() : empty(true), chunkID(-1), end(false), records{} {
    }

    /**
     * Clear the output chunk to reset it to a default constructed output
     * chunk
     */
    void clear() {
        empty = true;
        chunkID = -1;
        end = false;
        records.clear();
        counters.resetCounters();
    }

    /**
     * Set the chunkID and mark this OutputChunk as filled in
     * @param id the required id for this chunk
     */
    void set(size_t id) {
        chunkID = id;
        empty = false;
    }

    /**
     * Static function to create an end output with the provided id. An end
     * output is *never* empty.
     * @param id the id for this chunk
     * @returns an end output chunk with the provided id
     */
    static OutputChunk makeEndOutput(size_t id) {
        OutputChunk out = OutputChunk();
        out.setEnd(id);
        return out;
    }

    void addSingleEndRecord(std::vector<TextOcc>&& occ) {
        records.emplace_back(std::move(occ));
    }

    /**
     * @returns true if this is an end OutputChunk
     */
    bool isEnd() const {
        return end;
    }

    /**
     * @returns the records for this chunk
     */
    const std::vector<OutputRecord>& getRecords() const {
        return records;
    }

    /**
     * @returns the records for this chunk
     */
    std::vector<OutputRecord>& getRecordsNonConst() {
        return records;
    }

    /**
     * Reserve space for the records
     */
    void reserve(size_t s) {
        records.reserve(s);
    }

    // Move constructor
    OutputChunk(OutputChunk&& other) noexcept
        : empty(other.empty), chunkID(other.chunkID), end(other.end),
          records(std::move(other.records)),
          counters(std::move(other.counters)) {
        other.clear();
    }

    // Move assignment operator
    OutputChunk& operator=(OutputChunk&& other) noexcept {
        if (this != &other) {
            empty = other.empty;
            chunkID = other.chunkID;
            end = other.end;
            records = std::move(other.records);
            counters = std::move(other.counters);

            other.clear();
        }
        return *this;
    }

    Counters& getMutableCounters() {
        return counters;
    }

    // Delete copy constructor
    OutputChunk(const OutputChunk&) = delete;

    ~OutputChunk() {
        clear();
    }
};

// ============================================================================
// OUTPUT WRITER
// ============================================================================

/**
 * Class that handles the writing to disk of occurrences.
 */
class OutputWriter {
  private:
    std::string writeFile; // The file to write to
    FileType fileType;     // The type of the writeFile
    bool reorder = false;  // whether we keep the original order
    size_t maxNumChunks; // the maximum no. of unprocessed chunks stored in the
                         // writer
    std::thread wThread; // thread to write to disk
    std::string headerFile; // The file to get the info about the index for
                            // the header from
    std::string commandLineParameters; // The command line parameters with
                                       // which the program was started

    // output thread variables (protect by outputMutex)
    std::mutex outputMapMutex;           // output thread mutex
    std::condition_variable outputReady; // output blocks ready condition
    std::condition_variable outputSpace; // output blocks ready condition
    phmap::btree_map<size_t, OutputChunk>
        outputChunks;   // output blocks to be written
    size_t nextChunkID; // next chunk to be written

    length_t terminationMessagesSent = 0;
    length_t terminationMessagesReceived = 0;
    length_t threads = 0;

    // variables for logging
    length_t logIntervalRecords =
        1024 * 8; // after how many processed records do we log
    std::chrono::time_point<std::chrono::steady_clock> lastLogTime;
    const std::chrono::seconds logIntervalSeconds =
        std::chrono::seconds(60); // after how many seconds must we log
    const std::chrono::seconds minimumLogSeconds =
        std::chrono::seconds(5); // minimum time between two logs

    std::chrono::time_point<std::chrono::steady_clock> startTime;

    /**
     * Writes chunks to file
     */
    void writerThread();

    /**
     * Log the processing of the records
     * @param processedRecords the number of records that have been processed
     * @param alwaysLog whether to always log, independent of how much time has
     * passed/records have been processed. (Default = false)
     */
    void logProcess(length_t processedRecords, bool alwaysLog = false);

  public:
    /**
     * Constructor
     * @param writeFile the file to write to
     * @param headerFile the file to get the info about the index for the
     * header from
     * @param commandLineParameters the command line parameters with which the
     * program was started
     * @param reorder whether to keep original order
     */
    OutputWriter(const std::string& writeFile, const std::string& headerFile,
                 const std::string& commandLineParameters, bool reorder);

    /**
     * Start writer thread.
     *
     * @param maxNumChunks the maximum number of chunks to store concurrently in
     * the writer. Must be at least the number of threads.
     * @param threads the number of threads
     */
    void start(const size_t maxNumChunks, const size_t threads);
    /**
     * Join the writer thread
     */
    void joinWriterThread();

    /**
     * Commit chunk of data to be written
     * @param chunk Buffer in which to store the records (output)
     */
    void commitChunk(OutputChunk& chunk);

    /**
     * Adds a termination message at id chunkID if no termination has been
     * set yet
     */
    void sendTermination(size_t chunkID);
};
#endif
