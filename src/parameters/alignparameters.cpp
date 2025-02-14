#include "alignparameters.h"
#include "../indexinterface.h"
#include "../logger.h"
#include "../searchstrategy.h"
#include "parameters.h"
#include <algorithm>
#include <thread> // for thread

#include <vector>

/**
 * Option for the command line arguments considering the log file to be used.
 */
class LogFileOption : public ParameterOption {
  public:
    LogFileOption() : ParameterOption("l", "log-file", true, STRING, OUTPUT) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        params.logFile = arg;
    }

    std::string getDescription() const override {
        return "Path to the log file. Default is stdout.";
    }
};

/**
 * Option for the command line arguments considering the first reads file to be
 * used.
 */
class FirstReadsOption : public ParameterOption {
  public:
    FirstReadsOption()
        : ParameterOption("f", "reads-file", true, STRING, REQUIRED) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        params.firstReadsFile = arg;
    }

    std::string getDescription() const override {
        return "Path to the reads file.";
    }
};

/**
 * Option for the command line arguments considering the reference to be
 * used.
 */
class ReferenceOption : public ParameterOption {
  public:
    ReferenceOption()
        : ParameterOption("r", "reference", true, STRING, REQUIRED) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        params.base = arg;
    }

    std::string getDescription() const override {
        return "Path to the basename of the index.";
    }
};

/**
 * Option for the command line arguments considering the output file to be
 * used.
 */
class OutputFileOption : public ParameterOption {
  public:
    OutputFileOption()
        : ParameterOption("o", "output-file", true, STRING, OUTPUT) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        params.outputFile = arg;

        // check if the directory of the output file exists
        size_t found = params.outputFile.find_last_of("/\\");
        if (found == std::string::npos) {
            return;
        }
        std::string directory = params.outputFile.substr(0, found);
        struct stat info;
        if (stat(directory.c_str(), &info) != 0 || !(info.st_mode & S_IFDIR)) {
            std::string withoutDirectory = params.outputFile.substr(found + 1);
            logger.logWarning("Directory " + directory +
                              " does not exist. Setting output file to " +
                              withoutDirectory + " in the current directory.");
            params.outputFile = withoutDirectory;
        }
    }

    std::string getDescription() const override {
        return "Path to the output file. Should be *.out. Default is "
               "tagger.out.";
    }
};

class MinMEMLengthOption : public ParameterOption {
  public:
    MinMEMLengthOption()
        : ParameterOption("L", "min-mem-length", true, STRING, ADVANCED) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        try {
            params.minMEMLength = std::stoi(arg);
        } catch (...) {
            logger.logWarning(
                "Minimum MEM length should be a positive integer" +
                ignoreMessage());
        }
        if (params.minMEMLength < 0) {
            logger.logWarning(
                "Minimum MEM length should be a positive integer" +
                ignoreMessage());
        }
    }

    std::string getDescription() const override {
        return "The minimum MEM length to be detected. Default is " +
               std::to_string(DEFAULT_L) + ".";
    }
};

#ifdef TAG_ARRAY_SUBSAMPLING
class MaxLFStepsOption : public ParameterOption {
  public:
    MaxLFStepsOption()
        : ParameterOption("M", "max-LF-steps", true, STRING, ADVANCED) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        try {
            params.maxLFSteps = std::stoi(arg);
        } catch (...) {
            logger.logWarning(
                "Maximum number of LF steps should be a positive integer" +
                ignoreMessage());
        }
        if (params.maxLFSteps < 0) {
            logger.logWarning(
                "Maximum number of LF steps should be a positive integer" +
                ignoreMessage());
        }
    }

    std::string getDescription() const override {
        return "The maximum number of LF steps to be taken to find a sampled "
               "tag. Default is " + std::to_string(DEFAULT_MAX_LF) + ".";
    }
};
#endif // TAG_ARRAY_SUBSAMPLING

/**
 * Option for the command line arguments considering the k-mer size to be
 * used.
 */
class KmerSizeOption : public ParameterOption {
  public:
    KmerSizeOption()
        : ParameterOption("K", "kmer-size", true, INTEGER, ADVANCED) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        length_t defaultKmerSize = params.kmerSize;
        try {
            params.kmerSize = std::stoi(arg);
        } catch (...) {
            logger.logWarning("kmerSkipSize should be an integer." +
                              ignoreMessage());
        }
        if (params.kmerSize < 0 || params.kmerSize > 15) {
            logger.logWarning("kmerSkipSize should be in [0, 15]." +
                              ignoreMessage());
            params.kmerSize = defaultKmerSize;
        }
    }

    std::string getDescription() const override {
        return "The size of k-mers in the hash table. Default is 10.";
    }
};

/**
 * Option for the command line arguments considering the number of threads  to
 * be used.
 */
class ThreadsOption : public ParameterOption {
  public:
    ThreadsOption()
        : ParameterOption("t", "threads", true, INTEGER, PARALLELIZATION) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        try {
            params.nThreads = std::stoi(arg);
            if (params.nThreads <= 0) {
                logger.logWarning("Only positive values are allowed for the "
                                  "number of threads. Using 1 thread instead.");
                params.nThreads = 1;
            }

        } catch (...) {
            logger.logWarning(
                "Number of threads should be an integer. Set to 1 by default.");
            params.nThreads = 1;
        }
    }

    std::string getDescription() const override {
        return "The number of threads to be used. Default is 1.";
    }
};

/**
 * Option for the command line arguments to trigger the help.
 */
class HelpOption : public ParameterOption {
  public:
    HelpOption() : ParameterOption("h", "help", false, NONE, HELP) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        ParametersInterface::printHelp();
        exit(0);
    }

    std::string getDescription() const override {
        return "Print this help message.";
    }
};

/**
 * Option for the command line arguments to ensure that the output is in the
 * same order as the input.
 */
class ReorderOption : public ParameterOption {
  public:
    ReorderOption() : ParameterOption("R", "reorder", false, NONE, OUTPUT) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        params.reorder = true;
    }

    std::string getDescription() const override {
        return "Guarantees that output records are printed in the "
               "order corresponding to the order of reads in the original "
               "file. Setting this will cause tagger to be somewhat slower "
               "and use somewhat more memory.";
    }
};

const std::vector<std::shared_ptr<Option>> ParametersInterface::options = {
    std::make_shared<LogFileOption>(),
    std::make_shared<FirstReadsOption>(),
    std::make_shared<ReferenceOption>(),
    std::make_shared<OutputFileOption>(),
    std::make_shared<KmerSizeOption>(),
    std::make_shared<MinMEMLengthOption>(),
#ifdef TAG_ARRAY_SUBSAMPLING
    std::make_shared<MaxLFStepsOption>(),
#endif // TAG_ARRAY_SUBSAMPLING
    std::make_shared<ThreadsOption>(),
    std::make_shared<ReorderOption>(),
    std::make_shared<HelpOption>()
};

void ParametersInterface::printHeader() {
    std::cout << "Usage: tagger [OPTIONS]\n\n";
}

Parameters Parameters::processOptionalArguments(int argc, char** argv) {

    Parameters params;
    params.command = getCommand(argc, argv);

    params.parse(argc, argv);

    // sanity checks

    // check if threads does not exceed hardware options
    if (params.nThreads > std::thread::hardware_concurrency()) {
        std::stringstream ss;
        ss << "The entered number of threads: " << params.nThreads
           << " is higher than the number of threads "
              "available: "
           << std::thread::hardware_concurrency()
           << ". Setting number of threads to "
           << std::thread::hardware_concurrency();
        logger.logWarning(ss);
        params.nThreads = std::thread::hardware_concurrency();
    }

    if(params.minMEMLength < params.kmerSize){
        params.kmerSize = params.minMEMLength;
        logger.logWarning("Minimum MEM length is smaller than k-mer size. "
                          "Setting k-mer size to minimum MEM length of " +
                          std::to_string(params.minMEMLength));
    }

    return params;
}

std::unique_ptr<SearchStrategy>
Parameters::createStrategy(IndexInterface& index) const {
    std::unique_ptr<SearchStrategy> strategy;

    strategy.reset(new SearchStrategy(index));

    return strategy;
}
