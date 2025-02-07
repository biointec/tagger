# tagger  

**tagger** is a high-performance, run-length compressed metagenomic read classifier that utilizes SMEM-finding and a tag array for efficient classification. It supports multi-threaded execution and works with pre-built indexes for rapid read classification.  

## Features  

- **Run-length compression** using move tables for efficient storage and lookup  
- **SMEM-based classification** using SMEMs of length at least $L$  
- **Tag array toehold** to enable classification without locating  
- **Multi-threaded execution** for scalable performance  
- **Support for long and short reads**  

## Installation  

### Prerequisites  

Before you can build **tagger**, ensure you have the following installed:

- A C++ compiler supporting **C++11** or higher
- **CMake 3.14** or higher
- The **SDSL-lite library**: Follow [these installation instructions](https://github.com/simongog/sdsl-lite) to install it on your system.

### Install tagger  

Clone and build **tagger**:

```bash
git clone https://github.com/biointec/tagger.git
cd tagger
mkdir build
cd build
cmake ..
make -j$(nproc)  # Use all available cores for faster compilation
```

This will generate the `tagger` and `tagger_build` binaries.

### Dependencies  

tagger requires the [SDSL-lite library](https://github.com/simongog/sdsl-lite).  
You can follow their installation instructions to install it on your system.  

If you have installed SDSL to a non-standard location, you can point CMake to the installation location by adding `-DSDSL_INCLUDE_DIR=<path-to-sdsl>/` and `-DSDSL_LIBRARY=<path-to-sdsl-lib>` to the CMake command.  

## Usage  

### Index Construction  

Before classifying reads, you need to build an index from reference genomes.  

#### Build Modes  

tagger has two build modes:

1. **Standard Mode**: Uses `n*8` more memory than the second mode, where `n` is the total length of the reference sequences. This mode requires a higher memory usage but can be faster.
2. **Prefix-Free Parsing (PFP) Mode**: Uses [Big-BWT](https://gitlab.com/manzai/Big-BWT) for a more memory-efficient index construction. This mode requires Python 3.8 (or greater) and the `psutil` package. PFP Mode may be beneficial for larger genomes or when memory constraints are a concern.

**Standard Mode:**  
```bash
./tagger_build -f reference.fa -r index_name -t taggingcategories.txt
```  

Key options:
- `-f` / `--fasta-files`: Space-separated list of FASTA files containing reference sequences. Supports `.fa`, `.fna`, `.fasta` extensions. Can use `/path/to/directory/*.ext` to include all matching files.  
- `-F` / `--text-file-with-fasta`: File containing paths to FASTA files (one per line). Supports `.fa`, `.fna`, `.fasta` extensions.  
- `-r` / `--reference-base-name`: Name for the generated index (required).  
- `-t` / `--taggingcategories`: File containing the tagging categories (required).  
- `-l` / `--seed-length`: Seed length for replacing non-ACGT characters (default: 100, 0 means no seed).  
- `-M` / `--max-LF`: Space-separated list of maximum LF steps allowed to find a sampled tag (default: 0, meaning no subsampling).  

**PFP Mode:**  
```bash
bash tagger_build_pfp.sh -r <index_name> [-f <fasta_files>] [-F <fasta_file_list>] [-l <seedLength>] [-w <ws>] [-p <mod>] [-t <taggingCategoriesFile>] [-M <maxLFSteps>]
```

or  

```bash
bash columba_build_pfp.sh -r <index_name> [-F <fasta_file_list>] [-l <seedLength>] [-w <ws>] [-p <mod>] [-t <taggingCategoriesFile>] [-M <maxLFSteps>]
```

Required argument:
- `-r <index_name>`: Name/location of the index to be created.
- `-t <taggingCategoriesFile>`: File with tagging categories for the input file.

Optional arguments:
- `-f <fasta_files>`: Space-separated list of FASTA files.
- `-F <fasta_file_list>`: Path to a file containing a list of FASTA files, one per line.
- `-l <seedLength>`: Seed length for replacing non-ACGT characters (default: 100). 0 means no seed is used.
- `-w <ws>`: Window size for Big-BWT. If unset, Big-BWT will use its default window size.
- `-p <mod>`: Mod value for Big-BWT. If unset, Big-BWT will use its default mod value.
- `-M <maxLFSteps>`: Space-separated list of maximum LF steps for subsampling the tag samples.

**Note**: If you encounter hash collisions during Big-BWT processing, try increasing the window size (`-w`) and/or the mod value (`-p`) to resolve the issue.  

#### Seeded Replacement  

To improve memory usage when working with references that contain large blocks of non-ACGT characters, tagger uses seeded replacement. This ensures that non-ACGT characters are replaced in a predictable way to avoid the inefficiencies caused by randomness.

Use the `-l <seed_length>` option to specify the seed length (default is 100).

### Read Classification  

Once the index is built, classify reads using the `tagger` binary.  

**Basic usage:**  
```bash
./tagger -f reads.fq -r index_name
```  

Key options:
- `-f` / `--first-reads-file`: Path to the reads file.  
- `-r` / `--reference`: Path to the basename of the index.  
- `-t` / `--threads`: Number of threads (default: 1).  
- `-l` / `--log-file`: Path to the log file (default: stdout).  
- `-o` / `--output-file`: Path to the output file (default: `tagger.out`).  
- `-R` / `--reorder`: Ensures output records match input order (slower but uses more memory).  
- `-K` / `--kmer-size`: K-mer size for partitioning (default: 10).  
- `-L` / `--min-mem-length`: Minimum MEM length to detect (default: 25).  
- `-M` / `--max-LF-steps`: Maximum number of LF steps to find a sampled tag (only with the subsampling compiler option, default: 0).  

### Example Output  
When running read classification, the output file will contain the classification for each read. Example output might look like this:

```
ReadID  Tag
read_id_1   tag_1
read_id_2   tag_2
read_id_3   tag_1
...
```

## Compiler Options  

tagger provides several compiler options that you can enable or disable during the build process:

- **Use 32-bit types**  
  To use 32-bit types instead of the default 64-bit types, enable the `THIRTY_TWO` option:
  ```bash
  cmake -DTHIRTY_TWO=ON ..
  ```
  This is only possible if the index size is less than 2^32.

- **Enable tag array subsampling**  
  To enable tag array subsampling, use the `TAG_ARRAY_SUBSAMPLING` option:
  ```bash
  cmake -DTAG_ARRAY_SUBSAMPLING=ON ..
  ```
  Enabling this option introduces a memory-runtime trade-off where memory usage is reduced, but the classification process may take slightly longer due to the subsampling technique.

These options can be set during the `cmake` command to customize the build according to your system and needs.

## License and Dependencies  

This project is licensed under the **AGPL-3.0 License**. See the [LICENSE](./LICENSE) file for more details.  

### Dependencies  

In addition to the aforementioned dependencies, SDSL and Big-BWT, **tagger** also builds upon the foundation of the [b-move](https://github.com/biointec/b-move) implementation.  

tagger further relies on the following libraries:  
- **[{fmt} library](https://github.com/fmtlib/fmt)**: Falls under the exception of its license.  
- **[libsais](https://github.com/IlyaGrebnov/libsais)**: Licensed under the [Apache-2.0 License](./licenses_dependencies/Apache-2.0_LICENSE).  
- **[divsufsort](https://github.com/y-256/libdivsufsort)**: Licensed under the MIT license. See the [license](./licenses_dependencies/divsufsort_MIT_LICENSE) file in the repository.  
- **[parallel-hashmap](https://github.com/greg7mdp/parallel-hashmap)**: Licensed under the [Apache-2.0 License](./licenses_dependencies/Apache-2.0_LICENSE).

We acknowledge and appreciate the contributions made by these projects.

## Contact  

For any questions, feedback, or issues, please feel free to reach out via one of the following channels:

- **Email**: [Lore.Depuydt@UGent.be](mailto:Lore.Depuydt@UGent.be)  
- **GitHub Issues**: [github.com/biointec/tagger/issues](https://github.com/biointec/tagger/issues)  

We welcome contributions and feedback from the community!