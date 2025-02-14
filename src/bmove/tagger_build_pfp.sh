#!/bin/bash

# Script: tagger_build_pfp.sh
# Description: This script performs the tagger build process for PFP.
# Author: Lore Depuydt - lore.depuydt@ugent.be

# Capture start time
start_time=$(date +%s)

# We assume that this script is run from the build folder
tagger_build_exe="./tagger_build"
big_bwt_exe="./../external/Big-BWT/bigbwt"

taggingCategories=""

# Optional Big-BWT parameters
ws=0  # Default window size (unset means 0)
mod=0 # Default mod value (unset means 0)

# Array to store fasta files
fasta_files=()

# Array to store max LF values
maxLFValues=()

# Function to show usage
showUsage() {
	echo "Usage: $0 [-w <ws>] [-p <mod>] -r <index_name> [-f <fasta_files>] [-F <fasta_file_list>] [-t <taggingCategoriesFile>] [-M <maxLFSteps>]"
	echo
	echo "Required arguments:"
	echo "  -r <index_name>   Name/location of the index to be created."
	echo
	echo "Optional arguments:"
	echo "  -f <fasta_files>  Space-separated list of FASTA files."
	echo "  -F <fasta_file_list>  Path to a file containing a list of FASTA files, one per line."
	echo "  -w <ws>           Window size for Big-BWT. If unset, Big-BWT will use its default window size."
	echo "  -p <mod>          Mod value for Big-BWT. If unset, Big-BWT will use its default mod value."
	echo "  -t <taggingCategoriesFile>  File with tagging categories for the input file."
	echo "  -M <maxLFSteps>   Space-separated list of maximum LF steps for subsampling the tag samples."
}

# Function to run a command with /usr/bin/time -v and extract time and memory usage
# Usage: runCommand <command> [<args>...]
runCommand() {
	local command="$1"
	shift
	# Run the command and capture output, while measuring time and memory usage
	("$command" "$@") || {
		local status=$?
		echo "Error: Command '$command $@' failed with exit status $status." >&2
		exit $status
	}
}

# Function to parse command-line options
parseOptions() {
	# Parse command-line options
	while getopts ":r:f:F:w:p:t:M:" opt; do
		case $opt in
		t)
			taggingCategories=$OPTARG
			;;
		r)
			index_name=$OPTARG
			;;
		f)
			# Collect all subsequent arguments as the list of FASTA files
			fasta_files+=("$OPTARG")
			while [[ $OPTIND -le $# && ! ${!OPTIND} =~ ^- ]]; do
				fasta_files+=("${!OPTIND}")
				OPTIND=$((OPTIND + 1))
			done
			;;
		F)
			# Read each line in the specified file and add it to the fasta_files array
			if [[ -f $OPTARG ]]; then
				while IFS= read -r line; do
					fasta_files+=("$line")
				done <"$OPTARG"
			else
				echo "Error: File '$OPTARG' not found." >&2
				exit 1
			fi
			;;
		w)
			ws=$OPTARG
			;;
		p)
			mod=$OPTARG
			;;
		M)
			maxLFValues+=("$OPTARG")
			while [[ $OPTIND -le $# && ! ${!OPTIND} =~ ^- ]]; do
				maxLFValues+=("${!OPTIND}")
				OPTIND=$((OPTIND + 1))
			done
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			showUsage
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			showUsage
			exit 1
			;;
		esac
	done
	# Shift off the options and optional --
	shift $((OPTIND - 1))

	# Ensure required arguments are provided
	if [ -z "$index_name" ] || [ "${#fasta_files[@]}" -eq 0 ]; then
		showUsage
		exit 1
	fi
}

# Main script logic

# Parse command-line options
parseOptions "$@"

echo "Welcome to the tagger build process with prefix-free parsing!"
echo "-------------------------------------------------------------"
echo "Index name: $index_name"
echo "Input FASTA files:"
count=0
total_files=${#fasta_files[@]}
for file in "${fasta_files[@]}"; do
  if [ $count -ge 3 ]; then
    echo "  ... $((total_files - 3)) more files not shown"
    break
  fi
  echo "  $file"
  ((count++))
done
echo "Tagging categories file: $taggingCategories"
echo "Big-BWT window size: ${ws:-not set}"
echo "Big-BWT mod value: ${mod:-not set}"
echo "Max LF steps: ${maxLFValues[*]}"
echo "-------------------------------------------------------------"

# Start the preprocessing
echo "Start preprocessing the fasta file(s) with tagger..."
runCommand "$tagger_build_exe" --preprocess -r "$index_name" -f "${fasta_files[@]}" --taggingcategories "$taggingCategories" 
echo "Preprocessing done!"
echo "-------------------------------------------------------------"

base="${index_name}"

# Build Big-BWT command with optional -w and -p arguments
big_bwt_args=("$big_bwt_exe" -e -s -v "$base")
if [ "$ws" -gt 0 ]; then
	big_bwt_args+=(-w "$ws")
fi
if [ "$mod" -gt 0 ]; then
	big_bwt_args+=(-p "$mod")
fi

# Start the prefix-free parsing
echo "Start prefix-free parsing for the original string..."
runCommand "${big_bwt_args[@]}"
echo "Prefix-free parsing done!"
echo "-------------------------------------------------------------"

# Adjust the arguments for the reverse string
big_bwt_args=("$big_bwt_exe" -v "${base}.rev")
if [ "$ws" -gt 0 ]; then
	big_bwt_args+=(-w "$ws")
fi
if [ "$mod" -gt 0 ]; then
	big_bwt_args+=(-p "$mod")
fi

echo "Start prefix-free parsing for the reverse string..."
runCommand "${big_bwt_args[@]}"
echo "Prefix-free parsing done!"
echo "-------------------------------------------------------------"

# Start building the tagger index
echo "Start building the tagger index..."
pfp_args=("$tagger_build_exe" --pfp -r "$index_name")
if [ "${#maxLFValues[@]}" -gt 0 ]; then
	pfp_args+=(-M "${maxLFValues[@]}")
fi
runCommand "${pfp_args[@]}"
echo "tagger index built!"
echo "-------------------------------------------------------------"

# Remove the temporary files
echo "Remove temporary files..."
rm "${base}.bwt"
rm "${base}.rev.bwt"
rm "${base}.ssa"
rm "${base}.esa"
rm "${base}.log"
rm "${base}.rev.log"
rm "${base}"
rm "${base}.rev"
echo "Temporary files removed!"
echo "-------------------------------------------------------------"

# Capture end time
end_time=$(date +%s)

# Calculate total elapsed time
total_time=$((end_time - start_time))

echo "The tagger build process with prefix-free parsing is finished!"
echo "Total time elapsed: $total_time seconds."
