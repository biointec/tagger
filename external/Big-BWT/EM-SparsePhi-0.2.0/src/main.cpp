/**
 * @file    src/main.cpp
 * @section LICENCE
 *
 * This file is part of EM-SparsePhi v0.2.0
 * See: http://www.cs.helsinki.fi/group/pads/
 *
 * Copyright (C) 2016-2017
 *   Juha Karkkainen <juha.karkkainen (at) cs.helsinki.fi>
 *   Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <ctime>
#include <string>
#include <sstream>
#include <getopt.h>
#include <unistd.h>

#include "em_sparse_phi_src/em_sparse_phi.hpp"
#include "uint24.hpp"
#include "uint40.hpp"
#ifndef MAX40BITS
#include "uint48.hpp"
#include "uint56.hpp"
#endif

char *program_name;

void usage(int status) {
  printf(

"Usage: %s [OPTION]... FILE\n"
"Construct the LCP array for text stored in FILE.\n"
"\n"
"Mandatory arguments to long options are mandatory for short options too.\n"
"  -c, --charsize=SIZE     use symbols of SIZE bytes. Default: 1. Currently\n"
"                          supported values are: 1-8\n"
"  -f, --fast              enable faster computation at the price of higher disk\n"
"                          use (the algorithm in this mode is no longer inplace)\n"
"  -h, --help              display this help and exit\n"
"  -i, --intsize=SIZE      use integers of SIZE bytes. Int type needs to be wide\n"
"                          enough to encode any position in FILE. Default: 5.\n"
"                          Currently supported values are: 4-8\n"
"  -m, --mem=MEM           use MEM bytes of RAM for computation. Metric and IEC\n"
"                          suffixes are recognized, e.g., -l 10k, -l 1Mi, -l 3G\n"
"                          gives MEM = 10^4, 2^20, 3*10^6. Default: 3584Mi\n"
"  -o, --output=OUTFILE    specify output filename. Default: FILE.lcpX, where\n"
"                          X = integer size, see the -i flag\n"
"  -w, --overwrite         overwrite without asking an existing output file\n"
"  -s, --sa=SUFARRAY       specify the location of the suffix array of FILE.\n"
"                          Default: FILE.saX, X = integer size, see -i flag\n",
    program_name);

  std::exit(status);
}

bool file_exists(std::string filename) {
  std::FILE *f = std::fopen(filename.c_str(), "r");
  bool ret = (f != NULL);
  if (f != NULL) std::fclose(f);

  return ret;
}

template<typename int_type>
std::string intToStr(int_type x) {
  std::stringstream ss;
  ss << x;
  return ss.str();
}

template<typename int_type>
bool parse_number(char *str, int_type *ret) {
  *ret = 0;
  std::uint64_t n_digits = 0;
  std::uint64_t str_len = std::strlen(str);
  while (n_digits < str_len && std::isdigit(str[n_digits])) {
    std::uint64_t digit = str[n_digits] - '0';
    *ret = (*ret) * 10 + digit;
    ++n_digits;
  }

  if (n_digits == 0)
    return false;

  std::uint64_t suffix_length = str_len - n_digits;
  if (suffix_length > 0) {
    if (suffix_length > 2)
      return false;

    for (std::uint64_t j = 0; j < suffix_length; ++j)
      str[n_digits + j] = std::tolower(str[n_digits + j]);
    if (suffix_length == 2 && str[n_digits + 1] != 'i')
      return false;

    switch(str[n_digits]) {
      case 'k':
        if (suffix_length == 1)
          *ret *= 1000;
        else
          *ret <<= 10;
        break;
      case 'm':
        if (suffix_length == 1)
          *ret *= 1000000;
        else
          *ret <<= 20;
        break;
      case 'g':
        if (suffix_length == 1)
          *ret *= 1000000000;
        else
          *ret <<= 30;
        break;
      case 't':
        if (suffix_length == 1)
          *ret *= 1000000000000;
        else
          *ret <<= 40;
        break;
      default:
        return false;
    }
  }

  return true;
}

template<typename char_type>
void compute_lcp(
    std::string text_filename,
    std::string sa_filename,
    std::string output_filename,
    std::uint64_t ram_use,
    std::uint64_t text_offset_size,
    bool inplace_mode) {
  if (text_offset_size == 4)
    em_sparse_phi<char_type, std::uint32_t>(text_filename,
        sa_filename, output_filename, ram_use, inplace_mode);
  else if (text_offset_size == 5)
    em_sparse_phi<char_type, uint40>(text_filename,
        sa_filename, output_filename, ram_use, inplace_mode);
#ifndef MAX40BITS
  else if (text_offset_size == 6)
    em_sparse_phi<char_type, uint48>(text_filename,
        sa_filename, output_filename, ram_use, inplace_mode);
  else if (text_offset_size == 7)
    em_sparse_phi<char_type, uint56>(text_filename,
        sa_filename, output_filename, ram_use, inplace_mode);
  else if (text_offset_size == 8)
    em_sparse_phi<char_type, std::uint64_t>(text_filename,
        sa_filename, output_filename, ram_use, inplace_mode);
#endif
  else {
    fprintf(stderr, "\nError: compute_lcp: unsupported int type!\n");
    std::exit(EXIT_FAILURE);
  }
}

void compute_lcp(
    std::string text_filename,
    std::string sa_filename,
    std::string output_filename,
    std::uint64_t ram_use,
    std::uint64_t char_size,
    std::uint64_t text_offset_size,
    bool inplace_mode) {

  if (char_size == 1)
    compute_lcp<std::uint8_t>(text_filename, sa_filename,
        output_filename, ram_use, text_offset_size, inplace_mode);
  else if (char_size == 2)
    compute_lcp<std::uint16_t>(text_filename, sa_filename,
        output_filename, ram_use, text_offset_size, inplace_mode);
  else if (char_size == 3)
    compute_lcp<uint24>(text_filename, sa_filename,
        output_filename, ram_use, text_offset_size, inplace_mode);
  else if (char_size == 4)
    compute_lcp<std::uint32_t>(text_filename, sa_filename,
        output_filename, ram_use, text_offset_size, inplace_mode);
  else if (char_size == 5)
    compute_lcp<uint40>(text_filename, sa_filename,
        output_filename, ram_use, text_offset_size, inplace_mode);
#ifndef MAX40BITS
  else if (char_size == 6)
    compute_lcp<uint48>(text_filename, sa_filename,
        output_filename, ram_use, text_offset_size, inplace_mode);
  else if (char_size == 7)
    compute_lcp<uint56>(text_filename, sa_filename,
        output_filename, ram_use, text_offset_size, inplace_mode);
  else if (char_size == 8)
    compute_lcp<std::uint64_t>(text_filename, sa_filename,
        output_filename, ram_use, text_offset_size, inplace_mode);
#endif
  else {
    fprintf(stderr, "\nError: compute_lcp: unsupported char type!\n");
    std::exit(EXIT_FAILURE);
  }
}

int main(int argc, char **argv) {
  srand(time(0) + getpid());
  program_name = argv[0];

  static struct option long_options[] = {
    {"charsize", required_argument, NULL, 'c'},
    {"fast",     no_argument,       NULL, 'f'},
    {"help",     no_argument,       NULL, 'h'},
    {"intsize",  required_argument, NULL, 'i'},
    {"mem",      required_argument, NULL, 'm'},
    {"output",   required_argument, NULL, 'o'},
    {"overwrite",no_argument,       NULL, 'w'},
    {"sa",       required_argument, NULL, 's'},
    {NULL,       0,                 NULL, 0}
  };

  std::uint64_t ram_use = 3584UL << 20;
  std::uint64_t char_size = 1;
  std::uint64_t int_size = 5;
  std::string sa_filename("");
  std::string output_filename("");
  bool inplace_mode = true;
  bool check_overwrite = true;

  // Parse command-line options.
  int c;
  while ((c = getopt_long(argc, argv, "c:fwhi:m:o:s:",
          long_options, NULL)) != -1) {
    switch(c) {
      case 'c':
        char_size = std::atol(optarg);
        if (char_size < 1 || char_size > 8) {
          fprintf(stderr, "Error: invalid char size (%lu)\n\n", char_size);
          usage(EXIT_FAILURE);
        }
        break;
      case 'f':
        inplace_mode = false;
        break;
      case 'w':
        check_overwrite = false;
        break;
      case 'h':
        usage(EXIT_FAILURE);
        break;
      case 'i':
        int_size = std::atol(optarg);
        if (int_size < 4 || int_size > 8) {
          fprintf(stderr, "Error: invalid int size (%lu)\n\n", int_size);
          usage(EXIT_FAILURE);
        }
        break;
      case 'm':
        {
          bool ok = parse_number(optarg, &ram_use);
          if (!ok) {
            fprintf(stderr, "Error: parsing phrase length "
                "limit (%s) failed\n\n", optarg);
            usage(EXIT_FAILURE);
          }
          if (ram_use == 0) {
            fprintf(stderr, "Error: invalid RAM limit (%lu)\n\n", ram_use);
            usage(EXIT_FAILURE);
          }
          break;
        }
      case 'o':
        output_filename = std::string(optarg);
        break;
      case 's':
        sa_filename = std::string(optarg);
        break;
      default:
        usage(EXIT_FAILURE);
    }
  }

  if (optind >= argc) {
    fprintf(stderr, "Error: FILE not provided\n\n");
    usage(EXIT_FAILURE);
  }

  // Parse the text filename.
  std::string text_filename = std::string(argv[optind++]);
  if (optind < argc) {
    fprintf(stderr, "Warning: multiple input files provided. "
    "Only the first will be processed.\n");
  }

  // Set default SA and output filenames (if not provided).
  if (sa_filename.empty())
    sa_filename = text_filename + ".sa" + intToStr(int_size);
  if (output_filename.empty())
    output_filename = text_filename + ".lcp" + intToStr(int_size);

  // Check for the existence of text and suffix array.
  if (!file_exists(text_filename)) {
    fprintf(stderr, "Error: input file (%s) does not exist\n\n",
        text_filename.c_str());
    usage(EXIT_FAILURE);
  }
  if (!file_exists(sa_filename)) {
    fprintf(stderr, "Error: suffix array (%s) does not exist\n\n",
        sa_filename.c_str());
    usage(EXIT_FAILURE);
  }

  if (check_overwrite and file_exists(output_filename)) {

    // Output file exists, should we proceed?
    char *line = NULL;
    std::uint64_t buflen = 0;
    std::int64_t len = 0L;

    do {
      printf("Output file (%s) exists. Overwrite? [y/n]: ",
          output_filename.c_str());
      if ((len = getline(&line, &buflen, stdin)) == -1) {
        printf("\nError: failed to read answer\n\n");
        std::fflush(stdout);
        usage(EXIT_FAILURE);
      }
    } while (len != 2 || (line[0] != 'y' && line[0] != 'n'));

    if (line[0] == 'n') {
      free(line);
      std::exit(EXIT_FAILURE);
    }
    free(line);
  }
  
  // Run the algorithm.
  compute_lcp(text_filename, sa_filename, output_filename,
      ram_use, char_size, int_size, inplace_mode);
}
