#include "buildparameters.h"
#include "../logger.h"
#include "parameters.h"
#include <array>
#include <cstdlib>
#include <fstream>
#include <string>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/stat.h>
#endif

bool directoryExists(const std::string& path) {
#ifdef _WIN32
    DWORD ftyp = GetFileAttributesA(path.c_str());
    if (ftyp == INVALID_FILE_ATTRIBUTES)
        return false; // something is wrong with your path!

    if (ftyp & FILE_ATTRIBUTE_DIRECTORY)
        return true; // this is a directory!

    return false; // this is not a directory!
#else
    struct stat info;
    return (stat(path.c_str(), &info) == 0 && (info.st_mode & S_IFDIR));
#endif
}

std::array<std::string, 6> BuildParameters::allowedExtensionsFasta = {
    ".fasta", ".fa", ".FASTA", ".FA", ".fna", ".FNA"};

/**
 * @brief Remove the extension of a fasta file.
 * @param baseFN The filename.
 */
void BuildParameters::removeFastaExtension(std::string& baseFN) {
    for (const auto& ext : allowedExtensionsFasta) {
        if (baseFN.size() >= ext.size() &&
            std::equal(ext.rbegin(), ext.rend(), baseFN.rbegin())) {
            baseFN = baseFN.substr(0, baseFN.size() - ext.size());
            break;
        }

        std::string gzExt = ext + ".gz";
        if (baseFN.size() >= gzExt.size() &&
            std::equal(gzExt.rbegin(), gzExt.rend(), baseFN.rbegin())) {
            baseFN = baseFN.substr(0, baseFN.size() - gzExt.size());
            break;
        }
    }
}

/**
 * @brief Check if the extension of the fasta file is valid. If the file has no
 * extension it is checked whether the file with one of the valid extensions
 * exists.
 * @param fastaFile The path to the fasta file.
 * @return A boolean indicating whether the extension is valid.
 */
bool BuildParameters::validFastaExtension(std::string& fastaFile) {
    // Check if the extension of the fasta file is valid

    size_t lastDot = fastaFile.find_last_of('.');
    std::string fastaExtension =
        (lastDot == std::string::npos) ? "" : fastaFile.substr(lastDot);

#ifdef HAVE_ZLIB
    if (fastaExtension == ".gz") {
        // find the extension of the file before the .gz
        size_t secondToLastDot = fastaFile.find_last_of('.', lastDot - 1);
        fastaExtension =
            (secondToLastDot == std::string::npos)
                ? ""
                : fastaFile.substr(secondToLastDot, lastDot - secondToLastDot);
    }
#endif

    if (!fastaExtension.empty()) {
        for (const auto& ext : allowedExtensionsFasta) {
            if (fastaExtension == ext) {
                return true;
            }
        }
    }

    // file has no extension, check if fastaFile  + ext exists
    for (const auto& ext : allowedExtensionsFasta) {
#ifdef HAVE_ZLIB
        // Check if .gz variation exists
        std::ifstream ifsGz(fastaFile + ext + ".gz");
        if (ifsGz) {
            logger.logInfo("Detected fasta extension: " + ext + ".gz for " +
                           fastaFile);
            fastaFile += ext + ".gz";
            return true;
        }
#endif
        std::ifstream ifs(fastaFile + ext);
        if (ifs) {
            logger.logInfo("Detected fasta extension: " + ext + " for " +
                           fastaFile);
            fastaFile += ext;
            return true;
        }
    }
    return false;
}

//==============================================================================
// BuildParameters options
//==============================================================================
class PFPOption : public BuildParameterOption {
  public:
    PFPOption() : BuildParameterOption("p", "pfp", false, NONE, ADVANCED) {
#ifndef BIG_BWT_USABLE
        noPrint = true;
#endif
    }

    void process(const std::string& arg,
                 BuildParameters& params) const override {
#ifdef BIG_BWT_USABLE
        params.pfp = true;
#else
        logger.logError("Prefix-free parsing (option" + getOptionTag() +
                        ") is not supported due to missing "
                        "dependencies (Python 3.8 or greater and psutil).");
        exit(1);
#endif
    }

    std::string getDescription() const override {
        return "Start the index construction after the prefix-free parsing "
               "step.";
    }
};

class PreprocessOnlyOption : public BuildParameterOption {
  public:
    PreprocessOnlyOption()
        : BuildParameterOption("P", "preprocess", false, NONE, ADVANCED) {
#ifndef BIG_BWT_USABLE
        noPrint = true;
#endif
    }

    void process(const std::string& arg,
                 BuildParameters& params) const override {
#ifdef BIG_BWT_USABLE
        params.preprocessOnly = true;
#else
        logger.logError("Preprocessing only (option" + getOptionTag() +
                        ") is not supported due to missing "
                        "dependencies (Python 3.8 or greater and psutil).");
        exit(1);
#endif
    }

    std::string getDescription() const override {
        return "Preprocess the fasta file only. Do not build the index.";
    }
};

class TaggingCategoriesOption : public BuildParameterOption {
  public:
    TaggingCategoriesOption()
        : BuildParameterOption("t", "taggingcategories", true, STRING,
                               TAGGING_OPTION) {
    }

    void process(const std::string& arg,
                 BuildParameters& params) const override {

        params.taggingCategories = arg;
    }

    std::string getDescription() const override {
        return "File containing the tagging categories for classification.";
    }
};

class MaxLFOption : public BuildParameterOption {
  public:
    MaxLFOption()
        : BuildParameterOption("M", "max-LF", true, INTEGER, TAGGING_OPTION) {
    }

    void process(const std::string& arg,
                 BuildParameters& params) const override {

        // Split the argument by spaces to extract multiple maxLF values
        std::istringstream iss(arg);
        std::string value;
        std::vector<int> maxLFValues;

        while (iss >> value) {
            try {
                int maxLF = std::stoi(value); // Convert string to integer
                if (maxLF < 0) {
                    logger.logWarning("Maximum number of LF steps should be a "
                                      "positive integer: " +
                                      value + " ignored.");
                } else {
                    maxLFValues.push_back(maxLF); // Add valid value
                }
            } catch (...) {
                logger.logWarning(
                    "Invalid value for maximum number of LF steps: " + value +
                    " ignored.");
            }
        }

        if (maxLFValues.empty()) {
            logger.logWarning("No valid maximum LF values provided.");
        } else {
            params.maxLFValues =
                maxLFValues; // Assign the vector to the parameter
        }
    }

    std::string getDescription() const override {
        return "Space-separated list of maximum LF steps allowed to find a "
               "sampled tag. Each value should be a positive integer. "
               "Default is " +
               std::to_string(DEFAULT_MAX_LF) +
               ". When a maximum LF value is set to 0, no subsampling is "
               "performed.";
    }
};

/**
 * Option for the command line arguments to trigger the help.
 */
class HelpOption : public BuildParameterOption {
  public:
    HelpOption() : BuildParameterOption("h", "help", false, NONE, HELP) {
    }

    void process(const std::string& arg,
                 BuildParameters& params) const override {
        ParametersInterface::printHelp();
        exit(0);
    }

    std::string getDescription() const override {
        return "Print this help message.";
    }
};

class ReferenceBaseNameOption : public BuildParameterOption {
  public:
    ReferenceBaseNameOption()
        : BuildParameterOption("r", "reference-base-name", true, STRING,
                               BUILD_OUTPUT_FILE) {
    }

    void process(const std::string& arg,
                 BuildParameters& params) const override {

        params.baseFN = arg;
    }

    std::string getDescription() const override {
        return "Name/location of the index to be created. The index will be "
               "stored in files with this name as base name. This option is "
               "required!";
    }
};

class FastaFilesOption : public BuildParameterOption {
  public:
    FastaFilesOption()
        : BuildParameterOption("f", "fasta-files", true, STRING,
                               BUILD_INPUT_FILE) {
    }

    void process(const std::string& arg,
                 BuildParameters& params) const override {

        // split the argument by spaces
        std::istringstream iss(arg);
        std::string fastaFile;
        while (iss >> fastaFile) {
            if (!BuildParameters::validFastaExtension(fastaFile)) {
                logger.logWarning("Invalid fasta file: " + fastaFile +
                                  " ignoring this file.");
            }
            params.fastaFiles.push_back(fastaFile);
        }
    }

    std::string getDescription() const override {
        return "Space separated list to the FASTA files containing the "
               "reference sequences. Allowed extensions are .fa, .fna, .fasta "
               "or capitalized versions. Can be a /path/to/directory/*.ext for "
               "any of the allowed extensions. This option can be omitted if "
               "the -F option is used.";
    }
};

class TextFileWithFastaOption : public BuildParameterOption {
  public:
    TextFileWithFastaOption()
        : BuildParameterOption("F", "text-file-with-fasta", true, STRING,
                               BUILD_INPUT_FILE) {
    }

    void process(const std::string& arg,
                 BuildParameters& params) const override {

        std::string fileName = arg;
        std::ifstream ifs(fileName);
        if (!ifs) {
            logger.logWarning("Cannot open file: " + fileName + ". " +
                              ignoreMessage());
            return;
        }

        std::string fastaFile;
        while (std::getline(ifs, fastaFile)) {
            if (!BuildParameters::validFastaExtension(fastaFile)) {
                logger.logWarning("Invalid fasta file: " + fastaFile +
                                  " ignoring this file.");
            }
            params.fastaFiles.push_back(fastaFile);
        }
    }

    std::string getDescription() const override {
        return "A text file containing the paths to the FASTA files with the "
               "reference sequences. Allowed extensions are .fa, .fna, .fasta "
               "or capitalized versions. Each file path should be on a "
               "separate line. This option can be omitted if the -f option is "
               "used.";
    }
};

const std::vector<std::shared_ptr<Option>> ParametersInterface::options = {
    std::make_shared<PFPOption>(),
    std::make_shared<PreprocessOnlyOption>(),
    std::make_shared<TaggingCategoriesOption>(),
    std::make_shared<MaxLFOption>(),
    std::make_shared<ReferenceBaseNameOption>(),
    std::make_shared<FastaFilesOption>(),
    std::make_shared<TextFileWithFastaOption>(),
    std::make_shared<HelpOption>()};

BuildParameters BuildParameters::parse(int argc, char** argv) {
    // set the default values
    BuildParameters params;
    params.ParametersInterface::parse(argc, argv);

    // sanity checks
    // if not BIG_BWT_USABLE, give an error when the
    // params.pfp is true
#ifndef BIG_BWT_USABLE
    if (params.pfp || params.preprocessOnly) {
        // this is not valid so return, an error message should already be
        // printed
        return params;
    }
#endif

    bool notPFP = !params.pfp;

    if (notPFP && params.fastaFiles.empty()) {
        logger.logError("No reference files provided. Please provide reference "
                        "fasta files with the -f or -F flag!");
        return params;
    }

    if (params.baseFN.empty()) {
        if (params.fastaFiles.size() == 1) {
            logger.logWarning(
                "No reference index name provided with -r, using the sole "
                "reference file name as index base name.");
            params.baseFN = params.fastaFiles[0];
            removeFastaExtension(params.baseFN);
        } else {
            logger.logError("No reference index name provided, provide one "
                            "with the -r flag!");
            return params;
        }
    }

    if(params.taggingCategories.empty()){
        logger.logError("No tagging categories file provided. Please provide a "
                        "tagging categories file with the -t flag!");
        exit(1);
    }

    logger.logInfo("Index base filename: " + params.baseFN);

    // check if the base filename is in a directory that exists
    std::string baseDir;
    std::size_t lastSlashPos =
        params.baseFN.find_last_of("/\\"); // Check for both '/' and '\'

    if (lastSlashPos != std::string::npos) {
        baseDir = params.baseFN.substr(0, lastSlashPos);
    }

    if (!baseDir.empty()) {
        if (!directoryExists(baseDir)) {
            logger.logError(
                "The directory of the base filename does not exist: " +
                baseDir);
            return params;
        }
    }

    if(params.maxLFValues.empty()){
        params.maxLFValues.push_back(0);
        params.maxLFValues.push_back(DEFAULT_MAX_LF);
    }

    // end of logic reached, parameters are valid!
    params.valid = true;

    return params;
}

void ParametersInterface::printHeader() {
    std::cout << "Usage: tagger_build [OPTIONS]\n\n";
}
