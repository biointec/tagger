#ifndef ALIGN_PARAMETERS_H
#define ALIGN_PARAMETERS_H

#include "../definitions.h"
#include "parameters.h"

// forward declaration
class SearchStrategy;
class IndexInterface;

/**
 * Struct Parameters contains all the parameters that can be set by the user
 * Provides a method to process the optional arguments.
 * Provides a method to create a search strategy based on the parameters.
 * Provides a method to get the maximum distance or identity based on the mode.
 */
struct Parameters : public ParametersInterface {

    // Input output parameters
    std::string firstReadsFile = "";  // path to first reads file
    std::string base = "";            // path to the basename of the index
    std::string outputFile = "tagger.out"; // path to the output file
    std::string logFile = "";                     // Path to the log file
    std::string command = ""; // The command that was used to run the program
    bool reorder = false;          // Ensure order is same as original file

    // Numerical parameters
    length_t nThreads = 1;    // the number of threads to be used
    length_t kmerSize = 10;   // The size of k-mers in the hash table (used as
                              // seeds during partitioning)
    length_t minMEMLength = DEFAULT_L; // The minimum MEM length for the MEM finder
#ifdef TAG_ARRAY_SUBSAMPLING
    length_t maxLFSteps = DEFAULT_MAX_LF; // The maximum number of LF steps to find a tag
#endif
    /**
     * Get the command that was used to run the program
     */
    static std::string getCommand(int argc, char** argv) {
        std::string command;
        for (int i = 0; i < argc; ++i) {
            command += argv[i];
            if (i < argc - 1) {
                command += " ";
            }
        }
        return command;
    }

    /**
     * Process the optional arguments given by the user
     * @param argc The number of arguments
     * @param argv The arguments
     * @return The parameters struct
     */
    static Parameters processOptionalArguments(int argc, char** argv);

    /**
     * Create a search strategy based on the parameters
     * @param index Pointer to he index to search in
     * @return The search strategy
     */
    std::unique_ptr<SearchStrategy> createStrategy(IndexInterface& index) const;
};

class ParameterOption : public Option {
  public:
    ParameterOption(const std::string& shortOpt, const std::string& longOpt,
                    bool hasArgument, ArgumentType argumentType,
                    OptionType optionType)
        : Option(shortOpt, longOpt, hasArgument, argumentType, optionType) {
    }

    void process(const std::string& arg,
                 ParametersInterface& params) const override {
        Parameters& p = dynamic_cast<Parameters&>(params);
        process(arg, p);
    };

    virtual void process(const std::string& arg, Parameters& params) const = 0;
};
#endif