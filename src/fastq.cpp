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

#include "fastq.h"
#include "definitions.h" // for SINGLE_END, length_t, SUB_VERSION_NUMB...
#include "logger.h"      // for Logger, logger
#include "seqfile.h"     // for SeqFile, (anonymous), FileType, UNKNOW...
#include "util.h"        // for Util

#include <algorithm> // for transform, max, copy
#include <cassert>   // for assert
#include <fstream>   // for operator<<, basic_ostream, stringstream
#include <iomanip>   // for operator<<, setprecision
#include <numeric>   // for accumulate
#include <parallel_hashmap/btree.h> // for btree_map, btree_iterator, btree...
#include <ratio>                    // for ratio
#include <stdexcept>                // for runtime_error
#include <tuple>                    // for tie, tuple

using namespace std;

// ============================================================================
// FASTQ RECORD
// ============================================================================

bool SequenceRecord::readFromFileFASTQ(SeqFile& inputFile) {
    // read sequence identifier
    char c = inputFile.peekCharacter();

    // end of file might be reached
    if (!inputFile.good())
        return false;
    if (c != '@')
        throw ios::failure("File " + inputFile.getFileName() +
                           " doesn't appear to be in FastQ format");
    inputFile.getLine(seqID);
    if (!seqID.empty() && seqID.back() == '\n')
        seqID.pop_back();

    // read the actual read
    inputFile.getLine(read);
    if (!read.empty() && read.back() == '\n')
        read.pop_back();

    // read the + line
    string dummy;
    inputFile.getLine(dummy);
    // read the quality scores
    inputFile.getLine(qual);
    if (!qual.empty() && qual.back() == '\n')
        qual.pop_back();

    cleanUpRecord();
    ReadBundle::makeReverseComplement();

    return !read.empty();
}

bool SequenceRecord::readFromFileFASTA(SeqFile& inputFile) {
    clear();

    // end of file might be reached
    if (!inputFile.good())
        return false;
    // read sequence identifier
    char c = inputFile.peekCharacter();

    if (c != '>')
        throw ios::failure("File doesn't appear to be in Fasta format");
    inputFile.getLine(seqID);
    if (!seqID.empty() && seqID.back() == '\n')
        seqID.pop_back();

    // read the actual read line by line until a new sequence identifier is
    // found or the end of the file is reached
    string line;
    while (inputFile.good()) {
        c = inputFile.peekCharacter();
        if (c == '>') {
            break;
        }
        inputFile.getLine(line);
        if (!line.empty() && line.back() == '\n')
            line.pop_back();
        read += line;
    }

    this->qual = "*";
    cleanUpRecord();
    ReadBundle::makeReverseComplement();

    return !read.empty();
}

// ============================================================================
// READ BLOCK CLASS
// ============================================================================

void ReadBlock::getNextChunk(vector<ReadBundle>& buffer) {

    size_t remainingElements = this->size() - nextChunkOffset;
    size_t bufferSize = std::min(remainingElements, targetChunkSize);

    buffer.resize(bufferSize);
    auto sourceStart = this->begin() + nextChunkOffset;
    auto sourceEnd = sourceStart + bufferSize;
    std::move(sourceStart, sourceEnd, buffer.begin());

    nextChunkOffset += bufferSize;
}

void ReadBlock::readFromFile(SeqFile& file, size_t targetBlockSize,
                             bool fastq) {
    // clear the block
    this->clear();
    nextChunkOffset = 0;

    // read in new contents
    SequenceRecord record(fastq);
    size_t thisBlockSize = 0;

    while (thisBlockSize < targetBlockSize) {
        if (!record.readFromFile(file))
            break;
        this->push_back(record);
        thisBlockSize += 1;
    }
}

// ============================================================================
// FASTQ READER
// ============================================================================

Reader::Reader(const string& filename1)
    : filename1(filename1), fileType1(UNKNOWN_FT) {
    // try to figure out the file format based on the ext
    string ext;

    tie(fileType1, baseFilename1) = getFileType(filename1);
    if (fileType1 == UNKNOWN_FT) {
        string msg = "don't know how to open file: \"" + filename1 +
                     "\"\nExpected one of the following extensions: "
                     ".fastq, .fq, .fasta or .fa  (or .gz variants thereof)\n";
        throw runtime_error(msg);
    }

    if (!Util::fileExists(filename1))
        throw runtime_error("cannot open file " + filename1);

    fastq1 = (fileType1 == FASTQ || fileType1 == FASTQ_GZ);
}

void Reader::readerThread() {
    // open the read file(s)
    SeqFile file1(fileType1 == FASTQ_GZ || fileType1 == FASTA_GZ);
    file1.open(filename1);

    // recheck the chunksize every 8 blocks to start
    size_t checkPoint = 8;
    size_t blockCounter = 0;

    while (file1.good()) {
        // A) wait until an input block is available
        unique_lock<mutex> inputLock(inputMutex);
        inputReady.wait(inputLock, [this] { return !inputBlocks.empty(); });

        // get the input block
        ReadBlock block = std::move(inputBlocks.back());
        inputBlocks.pop_back();
        inputLock.unlock();

        block.setTargetChunkSize(targetChunkSize);

        // B) fill up the block (only this thread has access)
        // no mutexes are held by this thread at this point
        block.readFromFile(file1, targetBlockSize, fastq1);

        // empty block: end-of-file reached
        if (block.empty())
            break;

        // C) push the record block onto the worker queue
        unique_lock<mutex> workLock(workMutex);
        workBlocks.push_back(std::move(block));
        workReady.notify_all();
        workLock.unlock();

        // D) check if recalculation of target chunk size is necessary
        blockCounter++;
        if (blockCounter >= checkPoint && !processingTimes.empty()) {

            // acquire the processing time mutex
            unique_lock<mutex> processLock(processMutex);

            double averageProcessingTime =
                accumulate(processingTimes.begin(), processingTimes.end(),
                           0.0) /
                processingTimes.size();

            // check if the value should change
            bool change = false;
            auto prevSize = targetChunkSize;
            if (averageProcessingTime < lowEndProcessingTime ||
                averageProcessingTime > highEndProcessingTime) {

                int multFactor = (midProcessingTime / averageProcessingTime);
                targetChunkSize *= multFactor;
                targetChunkSize = max(targetChunkSize, 1ul);
                targetBlockSize = targetChunkSize * numberOfWorkerThreads;

                change = targetChunkSize != prevSize;
            }

            processingTimes.clear();
            processLock.unlock();

            if (change) {
                checkPoint = 8;
                logger.logDeveloper("New target chunk size: " +
                                    to_string(targetChunkSize));
            } else {
                checkPoint *= 2;
                logger.logDeveloper("Check point increased to " +
                                    to_string(checkPoint));
            }
            blockCounter = 0;
        }
    }

    file1.close();

    // send a termination message to the workers
    unique_lock<mutex> workLock(workMutex);
    workBlocks.push_back(ReadBlock()); // empty block == termination
    workReady.notify_all();
    workLock.unlock();

    logger.logDeveloper("Reader thread finished");
}

bool Reader::getNextChunk(vector<ReadBundle>& chunk, size_t& chunkID) {
    // wait until work becomes available
    unique_lock<mutex> workLock(workMutex);
    workReady.wait(workLock, [this] { return !workBlocks.empty(); });

    // check if it is a termination message (empty block)
    if (workBlocks.front().empty()) {
        chunkID = this->chunkID;
        return false;
    }

    // get a chunk of work (contents of chunk will be overwritten)
    workBlocks.front().getNextChunk(chunk);
    chunkID = this->chunkID++;

    // last chunk: move block from work queue back to the input queue
    if (!workBlocks.front().chunkAvailable()) {
        ReadBlock tmp = std::move(workBlocks.front());
        workBlocks.pop_front();
        workLock.unlock();

        unique_lock<mutex> inputLock(inputMutex);
        inputBlocks.push_back(std::move(tmp));
        inputReady.notify_one();
        inputLock.unlock();
    }

    assert(!chunk.empty());
    return true;
}

void Reader::startReaderThread(size_t targetChunkSize, size_t numWorkThreads) {
    this->targetChunkSize = targetChunkSize;
    this->numberOfWorkerThreads = numWorkThreads;
    this->targetBlockSize = targetChunkSize * numWorkThreads;

    chunkID = 0;

    // initialize (empty) blocks on the input stack
    inputBlocks.resize(numReadBlocks * 2);

    // start reader thread
    iThread = thread(&Reader::readerThread, this);
}

void Reader::joinReaderThread() {
    iThread.join();

    inputBlocks.clear();
    workBlocks.clear();
}

// ============================================================================
// OUTPUT WRITER
// ============================================================================

OutputWriter::OutputWriter(const string& writeFile, const string& headerFile,
                           const string& commandLineParameters, bool reorder)
    : writeFile(writeFile), fileType(UNKNOWN_FT), reorder(reorder),
      headerFile(headerFile), commandLineParameters(commandLineParameters),
      lastLogTime(chrono::steady_clock::now()),
      startTime(chrono::steady_clock::now()) {

    // try to figure out the file format based on the ext
    string ext;
    fileType = UNKNOWN_FT;

    if (writeFile.length() >= 4)
        ext = writeFile.substr(writeFile.length() - 4);
    transform(ext.begin(), ext.end(), ext.begin(), ::toupper);
    if (ext == ".OUT") {
        fileType = OUT;
    } 

    if (fileType == UNKNOWN_FT) {
        string msg = "don't know how to open file: \"" + writeFile +
                     "\"\nExpected the following extension: "
                     ".out\n";
        throw runtime_error(msg);
    }
}

void OutputWriter::start(const size_t size, const size_t threads) {
    maxNumChunks = size;
    terminationMessagesReceived = 0, terminationMessagesSent = 0;
    this->threads = threads;
    nextChunkID = 0; // waiting on zeroth chunk
    // start writer thread
    wThread = thread(&OutputWriter::writerThread, this);
}

void OutputWriter::joinWriterThread() {
    wThread.join();
    outputChunks.clear();
}

void OutputWriter::commitChunk(OutputChunk& chunk) {
    // wait until space becomes available to store the read block
    // there is *always* space to store the next-to-be-written block
    size_t chunkID = chunk.getChunkID();
    // acquire the mutex
    unique_lock<mutex> outputLock(outputMapMutex);
    outputSpace.wait(outputLock, [this, chunkID] {
        return (reorder) ? ((chunkID == nextChunkID) ||
                            (outputChunks.size() < maxNumChunks - 1))
                         : outputChunks.size() < maxNumChunks;
    });
    // move chunk into map
    outputChunks[chunkID] = std::move(chunk);

    bool wasEmpty = outputChunks.size() == 1;
    bool wasNext = chunkID == nextChunkID;

    outputLock.unlock();

    // signal the output thread if necessary
    if ((!reorder && wasEmpty) || wasNext) {
        outputReady.notify_one();
    }
}

void OutputWriter::logProcess(length_t processedRecords, bool alwaysLog) {
    auto now = chrono::steady_clock::now();
    auto timeElapsedSinceLastLog =
        chrono::duration_cast<chrono::seconds>(now - lastLogTime);

    if (alwaysLog || (timeElapsedSinceLastLog >= minimumLogSeconds &&
                      (processedRecords % logIntervalRecords == 0 ||
                       timeElapsedSinceLastLog >= logIntervalSeconds))) {

        auto timeElapsed =
            chrono::duration_cast<chrono::milliseconds>(now - startTime)
                .count();

        // timeElapsed cannot be zero as we exceed the minimum logTime
        double rate = (processedRecords / static_cast<double>(timeElapsed));

        stringstream ss;
        ss << "Processed " << processedRecords << " read records"
           << ". Rate: " << fixed << setprecision(0) << rate * 1000 << " "
           << "reads/second";

        logger.logInfo(ss);
        lastLogTime = now;
    }
}

void OutputWriter::writerThread() {
    // open the read file(s)
    SeqFile writeSeq(false);
    writeSeq.open(writeFile, WRITE);

    Counters counters;
    counters.resetCounters();

    writeSeq.writeLine("ReadID\tTag\n");

    length_t processedRecords = 0;

    while (true) {
        // A) wait until the next chunk is available
        OutputChunk chunk;
        {
            unique_lock<mutex> outputLock(outputMapMutex);
            outputReady.wait(outputLock, [this] {
                return !outputChunks.empty() &&
                       (!reorder ||
                        (outputChunks.begin()->first == nextChunkID));
            });

            // get the output chunk
            chunk = std::move(outputChunks.begin()->second);
            outputChunks.erase(outputChunks.begin());
            nextChunkID++; // not necessary if !reorder but does not matter
            outputSpace.notify_all();
            outputLock.unlock();
        }

        // check for termination message
        if (chunk.isEnd()) {
            terminationMessagesReceived++;
            if (terminationMessagesReceived == threads) {
                break;
            }
            continue;
        }

        // add the counters
        counters.addCounters(chunk.getMutableCounters());

        // B) write the chunk (only this thread has access)
        // no thread is held at this point
        for (const OutputRecord& record : chunk.getRecords()) {

            counters.inc(Counters::NUMBER_OF_READS);

            bool mapped = !record.outputOcc->empty() &&
                          record.outputOcc->front().isValid();
            counters.inc(Counters::MAPPED_READS, mapped);

            counters.inc(Counters::TOTAL_UNIQUE_MATCHES,
                         (mapped) ? record.outputOcc->size() : 0);

            // write the occurrences  to file
            for (const auto& match : *record.outputOcc) {
                writeSeq.writeLine(match.getOutputLine());
            }

            logProcess(++processedRecords);
        }
    }

    writeSeq.close();

    logProcess(processedRecords, true);
    logger.logInfo("Finished writing to output file: " + writeFile);

    // add dropped unique matches
    counters.inc(Counters::TOTAL_UNIQUE_MATCHES,
                 counters.get(Counters::DROPPED_UNIQUE_MATCHES));
    counters.reportStatistics();
}

void OutputWriter::sendTermination(size_t chunkID) {
    // send a termination message (empty block)
    unique_lock<mutex> outputLock(outputMapMutex);

    outputChunks[chunkID + terminationMessagesSent] =
        OutputChunk::makeEndOutput(chunkID + terminationMessagesSent);
    terminationMessagesSent++;

    // signal the output thread if necessary
    if (!reorder || chunkID == nextChunkID) {

        outputReady.notify_one();
    }

} // unique_lock goes out of scope
