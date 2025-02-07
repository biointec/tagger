/**
 * @file    src/em_sparse_phi_src/em_sparse_phi.hpp
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

#ifndef __SRC_EM_SPARSE_PHI_SRC_EM_SPARSE_PHI_HPP_INCLUDED
#define __SRC_EM_SPARSE_PHI_SRC_EM_SPARSE_PHI_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <string>
#include <algorithm>
#include <parallel/algorithm>
#include <omp.h>
#include <unistd.h>

//#define PRODUCTION_VER

#include "compute_sparse_plcp.hpp"
#include "compute_final_lcp.hpp"
#include "compute_lcp_lower_bounds.hpp"
#include "compute_lcp_delta.hpp"
#include "preprocess_lcp_delta_values.hpp"
#include "utils.hpp"


namespace em_sparse_phi_private {

template<typename text_offset_type>
void write_sparse_plcp_to_file(
    text_offset_type *sparse_plcp,
    std::uint64_t text_length,
    std::uint64_t sampling_rate,
    std::string output_filename,
    std::uint64_t &total_io_volume) {

  // Compute space PLCP array size.
  std::uint64_t size = (text_length + sampling_rate - 1) / sampling_rate;

  // Depending on sampling_rate we
  // chose smaller of two encodings.
  if (sampling_rate >= 40) {

    // Plain array encoding using sizeof(text_offset_type) *
    // text_length / sampling_rate bytes.
    // Print initial message and start the timer.
    fprintf(stderr, "  Write sparse PLCP (array) to disk: ");
    long double start = utils::wclock();

    // Write data to disk.
    utils::write_to_file(sparse_plcp, size, output_filename);

    // Update I/O volume.
    std::uint64_t io_vol = sizeof(text_offset_type) * size;
    total_io_volume += io_vol;

    // Print summary.
    long double write_time = utils::wclock() - start;
    fprintf(stderr, "time = %.2Lfs, I/O = %.2LfMiB/s, "
        "total I/O vol = %.2Lfbytes/input symbol\n",
        write_time, ((1.L * io_vol) / (1L << 20)) / write_time,
        (1.L * total_io_volume) / text_length);
  } else {

    // Space-efficient encoding using at most
    // text_length * (1 + 1 / sampling_rate) bits.
    // Print initial message and start the timer.
    fprintf(stderr, "  Write sparse PLCP (bitvector) to disk: ");
    long double start = utils::wclock();

    // Initialize the writer.
    typedef async_stream_writer<std::uint64_t> writer_type;
    writer_type *writer =
      new writer_type(output_filename, (1UL << 20), 2);  // XXX

    // Initialize counters.
    std::uint64_t buffer = 0;
    std::uint64_t filled = 0;
    std::uint64_t prev = 0;

    // Write data to disk.
    for (std::uint64_t j = 0; j < size; ++j) {
      std::uint64_t cur = (std::uint64_t)sparse_plcp[j];
      std::uint64_t diff = cur - std::max(0L,
          (std::int64_t)prev - (std::int64_t)sampling_rate);
      prev = cur;

      while (diff > 0) {
        ++filled;

        if (filled == 64) {
          writer->write(buffer);
          buffer = 0;
          filled = 0;
        }
        --diff;
      }

      buffer |= (1UL << filled);
      ++filled;

      if (filled == 64) {
        writer->write(buffer);
        buffer = 0;
        filled = 0;
      }
    }

    if (filled > 0)
      writer->write(buffer);

    // Update I/O volume.
    std::uint64_t io_vol = writer->bytes_written();
    total_io_volume += io_vol;

    // Print summary.
    long double write_time = utils::wclock() - start;
    fprintf(stderr, "time = %.2Lfs, I/O = %.2LfMiB/s, "
        "total I/O vol = %.2Lfbytes/input symbol\n",
        write_time, ((1.L * io_vol) / (1L << 20)) / write_time,
        (1.L * total_io_volume) / text_length);

    // Clean up.
    delete writer;
  }
}

template<typename text_offset_type>
void read_sparse_plcp_from_file(
    text_offset_type *sparse_plcp,
    std::uint64_t text_length,
    std::uint64_t sampling_rate,
    std::string input_filename,
    std::uint64_t &total_io_volume) {

  // Compute space PLCP array size.
  std::uint64_t size = (text_length + sampling_rate - 1) / sampling_rate;

  // Depending on sampling_rate we
  // chose smaller of two encodings.
  if (sampling_rate >= 40) {

    // Plain array encoding using sizeof(text_offset_type) *
    // text_length / sampling_rate bytes.
    // Print initial message and start the timer.
    fprintf(stderr, "  Read sparse PLCP (array) from disk: ");
    long double start = utils::wclock();

    // Read data from disk.
    utils::read_from_file(sparse_plcp, size, input_filename);

    // Update I/O volume.
    std::uint64_t io_vol = sizeof(text_offset_type) * size;
    total_io_volume += io_vol;

    // Print summary.
    long double read_time = utils::wclock() - start;
    fprintf(stderr, "time = %.2Lfs, I/O = %.2LfMiB/s, "
        "total I/O vol = %.2Lfbytes/input symbol\n",
        read_time, ((1.L * io_vol) / (1L << 20)) / read_time,
        (1.L * total_io_volume) / text_length);
  } else {

    // Space-efficient encoding using at most
    // text_length * (1 + 1 / sampling_rate) bits.
    // Print initial message and start the timer.
    fprintf(stderr, "  Read sparse PLCP (bitvector) from disk: ");
    long double start = utils::wclock();
    std::uint64_t io_vol = 0;

    // Initialize the reader.
    typedef async_stream_reader<std::uint64_t> reader_type;
    reader_type *reader =
      new reader_type(input_filename, (1UL << 20), 2);  // XXX

    // Initialize counters.
    std::uint64_t buffer = 0;
    std::uint64_t filled = 0;
    std::uint64_t prev = 0;

    // Read data from disk.
    for (std::uint64_t j = 0; j < size; ++j) {
      std::uint64_t diff = 0;
      while (true) {
        if (filled == 0) {
          buffer = reader->read();
          filled = 64;
        }

        if (buffer & 1) {
          --filled;
          buffer >>= 1;
          break;
        } else {
          --filled;
          buffer >>= 1;
          ++diff;
        }
      }

      std::uint64_t cur = std::max(0L,
          (std::int64_t)prev - (std::int64_t)sampling_rate) + diff;
      prev = cur;
      sparse_plcp[j] = (text_offset_type)cur;
    }

    // Stop I/O threads.
    reader->stop_reading();

    // Update I/O volume.
    io_vol +=
      reader->bytes_read();
    total_io_volume += io_vol;

    // Print summary.
    long double read_time = utils::wclock() - start;
    fprintf(stderr, "time = %.2Lfs, I/O = %.2LfMiB/s, "
        "total I/O vol = %.2Lfbytes/input symbol\n",
        read_time, ((1.L * io_vol) / (1L << 20)) / read_time,
        (1.L * total_io_volume) / text_length);

    // Clean up.
    delete reader;
  }
}

//-----------------------------------------------------------------------------
// NOTE: it is quite important to only use this algorithm in the case
//   when ram_use <= text_length (otherwise there is a better algorithm),
//   because only then the disk space usage is really we that promise.
//   For example, when plcp_sampling_rate = 5 (which is the smallest possible
//   value), the pairs (i, phi[i]) take up 2n bytes (when computing sparse
//   PLCP). They are then replaced (though in practice they are only deleted
//   later) by pairs (i, PLCP[i]), which take up another 2n, 4n in total.
//   This is OK, since at this point there is no LCP array so we have 5n
//   bytes of free disk space available. But the point here is that if
//   sparse_plcp_sampling_rate is smaller than 5, the whole argument above
//   does not work any and in fact we could use a lot more than 12n of disk
//   space (for the whole algorithm).
//
//   maybe we should even add a check if ram_use > text_length
//-----------------------------------------------------------------------------
template<typename char_type,
  typename text_offset_type>
void em_sparse_phi(
    std::string text_filename,
    std::string sa_filename,
    std::string output_filename,
    std::uint64_t ram_use,
    bool inplace_mode) {

  utils::initialize_stats();

  // Initialize pseudo-random number generator.
  srand(time(0) + getpid());

  // Empty page cache.
  utils::empty_page_cache(text_filename);
  utils::empty_page_cache(sa_filename);

  // Start the timer.
  long double start = utils::wclock();

  // Initialize the I/O volume.
  std::uint64_t total_io_volume = 0;

  // Compute some basic parameters of computation.
  std::uint64_t text_file_size = utils::file_size(text_filename);
  std::uint64_t text_length = text_file_size / sizeof(char_type);

  // Check if the text is non-empty.
  if (text_length == 0) {
    fprintf(stderr, "Error: the input file is empty!\n");
    std::exit(EXIT_FAILURE);
  }

  // Sanity check.
  if (text_file_size % sizeof(char_type)) {
    fprintf(stderr, "\nError: size of text file is not a "
        "multiple of sizeof(char_type)!\n");
    std::exit(EXIT_FAILURE);
  }

  // XXX should there be minimal RAM, 1MiB?

  // Lex partitioning needs
  //   text_length * sizeof(text_offset_type) +
  //   (3 * n * n * sizeof(char_type)) / ram_use bytes of I/O.
  // Text partitioning needs
  //   3 * text_length * sizeof(text_offset_type) +
  //   (n * n * sizeof(char_type)) / ram_use bytes of I/O.
  // Thus, text partitioning uses less I/O iff
  //   ram_use * sizeof(text_offset_type) <
  //   text_length * sizeof(char_type).
  enum partitioning_t { text_partitioning, lex_partitioning };
  partitioning_t partitioning_type = lex_partitioning;
  if (inplace_mode == true &&
      ram_use * sizeof(text_offset_type) <
      text_length * sizeof(char_type))
    partitioning_type = text_partitioning;

  // Compute the number and sizes of parts.
  std::uint64_t n_parts = 0;
  std::vector<std::uint64_t> part_sizes;
  if (inplace_mode == true) {

    // Hardcoded values determined using
    // the analysis of vbyte encoding size.
    n_parts = 3;
    part_sizes.push_back(0.395L * text_length);
    part_sizes.push_back(0.330L * text_length);
    part_sizes.push_back(text_length - part_sizes[0] - part_sizes[1]);
  } else {
    n_parts = 1;
    part_sizes.push_back(text_length);
  }

  // Compute max halfsegment size and buffer
  // sizes used when computing lcp deltas.
  static const std::uint64_t opt_overflow_ram_compute_lcp_delta = (1UL << 20);
  static const std::uint64_t opt_in_buf_ram_compute_lcp_delta = (32UL << 20);
  static const std::uint64_t opt_out_buf_ram_compute_lcp_delta = (4UL << 20);
  static const std::uint64_t opt_local_buf_ram_compute_lcp_delta = (64UL << 20);

  std::uint64_t max_halfsegment_ram =
    (std::uint64_t)((long double)ram_use * 0.45L);  // min acceptable value
  std::uint64_t overflow_ram_compute_lcp_delta =
    opt_overflow_ram_compute_lcp_delta;  // XXX passed to more than one function
  std::uint64_t in_buf_ram_compute_lcp_delta =
    opt_in_buf_ram_compute_lcp_delta;
  std::uint64_t out_buf_ram_compute_lcp_delta =
    opt_out_buf_ram_compute_lcp_delta;
  std::uint64_t local_buf_ram_compute_lcp_delta =
    opt_local_buf_ram_compute_lcp_delta;

  // If necessary, shrink buffers, so that they fit within the given RAM budget.
  {
    std::uint64_t total_buf_size_compute_lcp_delta =
      2 * overflow_ram_compute_lcp_delta +
      in_buf_ram_compute_lcp_delta +
      out_buf_ram_compute_lcp_delta +
      local_buf_ram_compute_lcp_delta;
    if (2 * max_halfsegment_ram +
        total_buf_size_compute_lcp_delta > ram_use) {
      std::uint64_t ram_budget = ram_use - 2 * max_halfsegment_ram;
      long double shrink_factor =
        (long double)ram_budget /
        (long double)total_buf_size_compute_lcp_delta;
      overflow_ram_compute_lcp_delta =
        (std::uint64_t)((long double)overflow_ram_compute_lcp_delta *
            shrink_factor);
      in_buf_ram_compute_lcp_delta =
        (std::uint64_t)((long double)in_buf_ram_compute_lcp_delta *
            shrink_factor);
      out_buf_ram_compute_lcp_delta =
        (std::uint64_t)((long double)out_buf_ram_compute_lcp_delta *
            shrink_factor);
      local_buf_ram_compute_lcp_delta =
        (std::uint64_t)((long double)local_buf_ram_compute_lcp_delta *
            shrink_factor);
    } else max_halfsegment_ram =
      (ram_use - total_buf_size_compute_lcp_delta) / 2;
  }

  std::uint64_t max_halfsegment_size =
    std::max(1UL, max_halfsegment_ram / sizeof(char_type));
  std::uint64_t n_halfsegments =
    (text_length + max_halfsegment_size - 1) / max_halfsegment_size;
  std::uint64_t n_different_halfseg_pairs =
    (n_halfsegments * (n_halfsegments + 1)) / 2;
  std::uint64_t overflow_size_compute_lcp_delta =
    std::max(1UL, overflow_ram_compute_lcp_delta / sizeof(char_type));

  // Compute PLCP sampling rate and size
  // of buffer used for each halfsegment pair.
  static const std::uint64_t opt_in_buf_ram = (32UL << 20);
  static const std::uint64_t opt_out_buf_ram = (4UL << 20);
  std::uint64_t opt_local_buf_ram = 0;
#ifdef _OPENMP
  if (partitioning_type == lex_partitioning)
    opt_local_buf_ram = (100UL << 20);
  else opt_local_buf_ram = (140UL << 20);
#else
  opt_local_buf_ram = (64UL << 20);
#endif

  std::uint64_t sparse_plcp_and_halfseg_buffers_ram =
    (std::uint64_t)((long double)ram_use * 0.9L);  // min acceptable value
  std::uint64_t sparse_plcp_ram =
    (std::uint64_t)((long double)ram_use * 0.2L);  // min acceptable value
  std::uint64_t halfseg_buffers_ram = n_different_halfseg_pairs * (4UL << 20);
  std::uint64_t in_buf_ram = opt_in_buf_ram;
  std::uint64_t out_buf_ram = opt_out_buf_ram;
  std::uint64_t local_buf_ram = opt_local_buf_ram;

  // Compute the actual size of sparse PLCP
  // array and update sparse_plcp_ram.
  std::uint64_t sparse_plcp_size = 0;
  std::uint64_t plcp_sampling_rate = 0;
  {
    std::uint64_t max_possible_plcp_size =
      std::max(1UL, sparse_plcp_ram / sizeof(text_offset_type));
    plcp_sampling_rate = std::max(4UL,
        (text_length + max_possible_plcp_size - 1) / max_possible_plcp_size);
    sparse_plcp_size =
      (text_length + plcp_sampling_rate - 1) / plcp_sampling_rate;
    sparse_plcp_ram = sparse_plcp_size * sizeof(text_offset_type);
  }

  // Try to utilize best the available RAM
  // and make sure evertything is as large
  // as possible.
  {
    std::uint64_t total_buf_ram = in_buf_ram + out_buf_ram + local_buf_ram;
    if (sparse_plcp_and_halfseg_buffers_ram + total_buf_ram > ram_use) {

      // Shrink reader/writer/local buffers.
      {
        std::uint64_t ram_budget = ram_use -
          sparse_plcp_and_halfseg_buffers_ram;
        long double shrink_factor =
          (long double)ram_budget / (long double)total_buf_ram;
        in_buf_ram =
          (std::uint64_t)((long double)in_buf_ram * shrink_factor);
        out_buf_ram =
          (std::uint64_t)((long double)out_buf_ram * shrink_factor);
        local_buf_ram =
          (std::uint64_t)((long double)local_buf_ram * shrink_factor);
      }

      if (sparse_plcp_ram +
          halfseg_buffers_ram >
          sparse_plcp_and_halfseg_buffers_ram) {

        // Decrease RAM for halfsegment buffers.
        halfseg_buffers_ram =
          sparse_plcp_and_halfseg_buffers_ram - sparse_plcp_ram;
      } else {

        // Try to increase the space effectively used by the sparse
        // PLCP array. Then left halfsegmetn buffers use the rest.
        sparse_plcp_ram =
          sparse_plcp_and_halfseg_buffers_ram - halfseg_buffers_ram;
        std::uint64_t max_possible_plcp_size =
          std::max(1UL, sparse_plcp_ram / sizeof(text_offset_type));
        plcp_sampling_rate = std::max(4UL,
            (text_length + max_possible_plcp_size - 1) /
            max_possible_plcp_size);
        sparse_plcp_size =
          (text_length + plcp_sampling_rate - 1) / plcp_sampling_rate;
        sparse_plcp_ram = sparse_plcp_size * sizeof(text_offset_type);
        halfseg_buffers_ram =
          sparse_plcp_and_halfseg_buffers_ram - sparse_plcp_ram;
      }
    } else {

      // Increase the space for sparse PLCP and halfseg buffers.
      sparse_plcp_and_halfseg_buffers_ram = ram_use - total_buf_ram;
      if (sparse_plcp_ram + halfseg_buffers_ram >
          sparse_plcp_and_halfseg_buffers_ram) {

        // Decrese RAM for halfseg buffers.
        halfseg_buffers_ram =
          sparse_plcp_and_halfseg_buffers_ram - sparse_plcp_ram;
      } else {

        // Try to increase the space effectively used by the sparse
        // PLCP array. Then left halfsegmetn buffers use the rest.
        sparse_plcp_ram =
          sparse_plcp_and_halfseg_buffers_ram - halfseg_buffers_ram;
        std::uint64_t max_possible_plcp_size =
          std::max(1UL, sparse_plcp_ram / sizeof(text_offset_type));
        plcp_sampling_rate = std::max(4UL,
            (text_length + max_possible_plcp_size - 1) /
            max_possible_plcp_size);
        sparse_plcp_size =
          (text_length + plcp_sampling_rate - 1) / plcp_sampling_rate;
        sparse_plcp_ram = sparse_plcp_size * sizeof(text_offset_type);
        halfseg_buffers_ram =
          sparse_plcp_and_halfseg_buffers_ram - sparse_plcp_ram;
      }
    }
  }

  // Turn paths absolute.
  text_filename = utils::absolute_path(text_filename);
  sa_filename = utils::absolute_path(sa_filename);
  output_filename = utils::absolute_path(output_filename);

  // Print summary of basic parameters.
  fprintf(stderr, "Running EM-SparsePhi v0.2.0\n");
  fprintf(stderr, "Timestamp = %s", utils::get_timestamp().c_str());
  fprintf(stderr, "Text filename = %s\n", text_filename.c_str());
  fprintf(stderr, "SA filename = %s\n", sa_filename.c_str());
  fprintf(stderr, "Output (LCP) filename = %s\n", output_filename.c_str());
  fprintf(stderr, "Text length = %lu (%.2LfMiB)\n", text_length,
      (1.L * text_length * sizeof(char_type)) / (1 << 20));
  fprintf(stderr, "RAM use = %lu bytes (%.2LfMiB)\n",
      ram_use, ram_use / (1024.L * 1024));
  fprintf(stderr, "sizeof(char_type) = %lu\n", sizeof(char_type));
  fprintf(stderr, "sizeof(text_offset_type) = %lu\n",
      sizeof(text_offset_type));
  fprintf(stderr, "Halfsegment size = %lu symbols (%.2LfMiB)\n",
      max_halfsegment_size,
      (1.L * max_halfsegment_size * sizeof(char_type)) / (1L << 20));
  fprintf(stderr, "Number of halfsegments = %lu\n", n_halfsegments);
  fprintf(stderr, "Number of halfsegment pairs = %lu\n",
      n_different_halfseg_pairs);
  fprintf(stderr, "PLCP sampling rate = %lu\n", plcp_sampling_rate);
  fprintf(stderr, "RAM for sparse PLCP = %lu bytes (%.2LfMiB)\n",
      sparse_plcp_ram,
      (1.L * sparse_plcp_ram) / (1L << 20));
  fprintf(stderr, "RAM for halfsegment buffers = %lu bytes (%.2LfMiB)\n",
      halfseg_buffers_ram,
      (1.L * halfseg_buffers_ram) / (1L << 20));

  fprintf(stderr, "Parallel mode = ");
#ifdef _OPENMP
  fprintf(stderr, "ON\n");
  fprintf(stderr, "Number of threads = %d\n", omp_get_max_threads());
#else
  fprintf(stderr, "OFF\n");
#endif
  fprintf(stderr, "Inplace mode = %s\n", inplace_mode ? "ON" : "OFF");
  if (inplace_mode == true) {
    fprintf(stderr, "Partitioning type = %s\n",
        partitioning_type == text_partitioning ? "text-order" : "lex-order");

    fprintf(stderr, "Part sizes:\n");
    for (std::uint64_t part_id = 0; part_id < n_parts; ++part_id)
      fprintf(stderr, "  #%lu: %lu (%.2Lf%%)\n", part_id + 1,
          part_sizes[part_id], (100.L * part_sizes[part_id]) / text_length);
  }
  fprintf(stderr, "\n\n");

  // Allocate info about the number of
  // items in each pair of halfsegments.
  std::uint64_t **pairs_count = NULL;
  if (partitioning_type == text_partitioning) {
    pairs_count = new std::uint64_t*[n_halfsegments];
    for (std::uint64_t j = 0; j < n_halfsegments; ++j) {
      pairs_count[j] = new std::uint64_t[n_halfsegments];
      std::fill(pairs_count[j], pairs_count[j] + n_halfsegments, 0UL);
    }
  }

  // Compute sparse PLCP array.
  text_offset_type *sparse_plcp = NULL;
  sparse_plcp = compute_sparse_plcp<char_type, text_offset_type>(
      text_filename, sa_filename, output_filename, text_length,
      max_halfsegment_size, pairs_count, plcp_sampling_rate,
      ram_use, total_io_volume);

  std::uint64_t lcp_sum = 0;
  std::uint64_t max_lcp = 0;
  std::string sparse_plcp_filename = output_filename +
    ".sparse_plcp" + utils::random_string_hash();

  if (partitioning_type == lex_partitioning) {
    std::string **delta_filenames = new std::string*[n_halfsegments];
    for (std::uint64_t i = 0; i < n_halfsegments; ++i) {
      delta_filenames[i] = new std::string[n_halfsegments];
      for (std::uint64_t j = i; j < n_halfsegments; ++j)
        delta_filenames[i][j] = output_filename +
          "." + utils::random_string_hash();
    }

    // Compute lcp delta values.
    for (std::uint64_t part_id = 0; part_id < n_parts; ++part_id) {

      // Compute part boundaries.
      std::uint64_t part_beg = 0;
      for (std::uint64_t i = 0; i < part_id; ++i)
        part_beg += part_sizes[i];
      std::uint64_t part_end = part_beg + part_sizes[part_id];

      // Start the timer.
      long double process_part_start = utils::wclock();
      if (inplace_mode == true)
        fprintf(stderr, "Compute LCP delta values, part %lu/%lu:\n",
            part_id + 1, n_parts);
      else
        fprintf(stderr, "Compute LCP delta values:\n");

      // Compute LCP lower bounds
      compute_lcp_lower_bounds_lex_partitioning(sa_filename, output_filename,
          text_length, ram_use, sparse_plcp, plcp_sampling_rate,
          max_halfsegment_size, part_beg, part_end, halfseg_buffers_ram,
          in_buf_ram, local_buf_ram, overflow_size_compute_lcp_delta,
          total_io_volume);

      // Write sparse PLCP to disk.
      if (part_id == 0)
        write_sparse_plcp_to_file(sparse_plcp, text_length,
            plcp_sampling_rate, sparse_plcp_filename, total_io_volume);

      // Deallocate sparse PLCP.
      utils::deallocate(sparse_plcp);

      // Compute LCP delta-values.
      compute_lcp_delta_lex_partitioning<char_type, text_offset_type>(
          text_filename, output_filename, text_length,
          overflow_size_compute_lcp_delta, max_halfsegment_size,
          (part_id == 0), in_buf_ram_compute_lcp_delta,
          out_buf_ram_compute_lcp_delta, local_buf_ram_compute_lcp_delta,
          delta_filenames, total_io_volume);

      // Load sparse PLCP from disk.
      sparse_plcp = utils::allocate_array<text_offset_type>(sparse_plcp_size);
      read_sparse_plcp_from_file(sparse_plcp, text_length,
          plcp_sampling_rate, sparse_plcp_filename, total_io_volume);

      // Print summary.
      fprintf(stderr, "  Summary: time = %.2Lfs, "
          "total I/O vol = %.2Lfbytes/input symbol\n",
          utils::wclock() - process_part_start,
          (1.L * total_io_volume) / text_length);
    }

    // Delete sparse PLCP file.
    utils::file_delete(sparse_plcp_filename);

    if (inplace_mode == true) {

      // Preprocess LCP delta values.
      std::string lcp_delta_filename =
        output_filename + "." + utils::random_string_hash();
      preprocess_lcp_delta_values_lex_partitioning(sa_filename,
          lcp_delta_filename, text_length, sparse_plcp, plcp_sampling_rate,
          max_halfsegment_size, halfseg_buffers_ram, in_buf_ram, out_buf_ram,
          local_buf_ram, overflow_size_compute_lcp_delta, delta_filenames,
          total_io_volume);

      // Compute the final LCP array.
      compute_final_lcp_inplace(sa_filename, lcp_delta_filename,
          output_filename, text_length, sparse_plcp, plcp_sampling_rate,
          in_buf_ram, out_buf_ram, local_buf_ram, lcp_sum, max_lcp,
          total_io_volume);
    } else {

      // Compute final LCP array using
      // simple (not inplace) method.
      compute_final_lcp(sa_filename, output_filename, text_length,
          sparse_plcp, plcp_sampling_rate, max_halfsegment_size,
          halfseg_buffers_ram, in_buf_ram, out_buf_ram, local_buf_ram,
          overflow_size_compute_lcp_delta, delta_filenames, lcp_sum,
          max_lcp, total_io_volume);
    }

    // Deallocate sparse PLCP.
    utils::deallocate(sparse_plcp);

  } else {

    // Compute the number of pairs to be processed in
    // each of the halfsegments pairs for all text parts.
    std::uint64_t ***items_per_halfseg_pair = new std::uint64_t**[n_parts];
    for (std::uint64_t part_id = 0; part_id < n_parts; ++part_id) {
      items_per_halfseg_pair[part_id] = new std::uint64_t*[n_halfsegments];

      for (std::uint64_t i = 0; i < n_halfsegments; ++i) {
        items_per_halfseg_pair[part_id][i] =
          new std::uint64_t[n_halfsegments];
        std::fill(items_per_halfseg_pair[part_id][i],
            items_per_halfseg_pair[part_id][i] + n_halfsegments, 0UL);
      }

      std::uint64_t items_left = part_sizes[part_id];
      for (std::uint64_t diff = 0; diff < n_halfsegments; ++diff) {
        for (std::uint64_t j = n_halfsegments; j > diff; --j) {
          std::uint64_t i = (j - 1) - diff;
          std::uint64_t count = std::min(items_left, pairs_count[i][j - 1]);
          items_per_halfseg_pair[part_id][i][j - 1] = count;
          pairs_count[i][j - 1] -= count;
          items_left -= count;

          if (items_left == 0)
            break;
        }

        if (items_left == 0)
          break;
      }
    }

    // Initialize filenames for files
    // containing LCP delta values.
    std::string ***delta_filenames = new std::string**[n_parts];
    for (std::uint64_t part_id = 0; part_id < n_parts; ++part_id) {
      delta_filenames[part_id] = new std::string*[n_halfsegments];
      for (std::uint64_t i = 0; i < n_halfsegments; ++i) {
        delta_filenames[part_id][i] = new std::string[n_halfsegments];

        for (std::uint64_t j = i; j < n_halfsegments; ++j)
          delta_filenames[part_id][i][j] =
            output_filename + "." + utils::random_string_hash();
      }
    }

    // Process all text parts.
    for (std::uint64_t part_id = 0; part_id < n_parts; ++part_id) {

      // Print initial message and start the timer.
      long double process_part_start = utils::wclock();
      fprintf(stderr, "Compute LCP delta values, part %lu/%lu:\n",
          part_id + 1, n_parts);

      // Compute LCP lower bounds.
      compute_lcp_lower_bounds_text_partitioning(sa_filename,
          output_filename, text_length, ram_use, sparse_plcp,
          plcp_sampling_rate, max_halfsegment_size, halfseg_buffers_ram,
          in_buf_ram, local_buf_ram, overflow_size_compute_lcp_delta,
          items_per_halfseg_pair, part_id, total_io_volume);

      // write sparse PLCP to disk.
      if (part_id == 0)
        write_sparse_plcp_to_file(sparse_plcp, text_length,
            plcp_sampling_rate, sparse_plcp_filename, total_io_volume);

      // Deallocate sparse PLCP.
      utils::deallocate(sparse_plcp);

      // Compute LCP delta-values.
      compute_lcp_delta_text_partitioning<char_type, text_offset_type>(
          text_filename, output_filename, text_length,
          overflow_size_compute_lcp_delta, max_halfsegment_size,
          in_buf_ram_compute_lcp_delta, out_buf_ram_compute_lcp_delta,
          local_buf_ram_compute_lcp_delta, delta_filenames, part_id,
          total_io_volume);

      // Load sparse PLCP from disk.
      sparse_plcp = utils::allocate_array<text_offset_type>(sparse_plcp_size);
      read_sparse_plcp_from_file(sparse_plcp, text_length,
          plcp_sampling_rate, sparse_plcp_filename, total_io_volume);

      // Print summary.
      fprintf(stderr, "  Summary: time = %.2Lfs, "
          "total I/O vol = %.2Lfbytes/input symbol\n",
          utils::wclock() - process_part_start,
          (1.L * total_io_volume) / text_length);
    }

    // Delete sparse PLCP file.
    utils::file_delete(sparse_plcp_filename);

    // Preprocess LCP delta values.
    std::string lcp_delta_filename =
      output_filename + "." + utils::random_string_hash();
    preprocess_lcp_delta_values_text_partitioning(sa_filename,
        lcp_delta_filename, text_length, sparse_plcp, plcp_sampling_rate,
        max_halfsegment_size, halfseg_buffers_ram, in_buf_ram, out_buf_ram,
        local_buf_ram, overflow_size_compute_lcp_delta, n_parts,
        items_per_halfseg_pair, delta_filenames, total_io_volume);

    // Compute the final LCP array.
    compute_final_lcp_inplace(sa_filename, lcp_delta_filename,
        output_filename, text_length, sparse_plcp, plcp_sampling_rate,
        in_buf_ram, out_buf_ram, local_buf_ram, lcp_sum, max_lcp,
        total_io_volume);

    // Clean up.
    utils::deallocate(sparse_plcp);
    for (std::uint64_t part_id = 0; part_id < n_parts; ++part_id) {
      for (std::uint64_t i = 0; i < n_halfsegments; ++i)
        delete[] items_per_halfseg_pair[part_id][i];
      delete[] items_per_halfseg_pair[part_id];
    }

    delete[] items_per_halfseg_pair;

    for (std::uint64_t part_id = 0; part_id < n_parts; ++part_id) {
      for (std::uint64_t i = 0; i < n_halfsegments; ++i)
        delete[] delta_filenames[part_id][i];
      delete[] delta_filenames[part_id];
    }

    delete[] delta_filenames;
  }

  // Clean up.
  if (pairs_count != NULL) {
    for (std::uint64_t j = 0; j < n_halfsegments; ++j)
      delete[] pairs_count[j];
    delete[] pairs_count;
  }

  // Print summary.
  long double total_time = utils::wclock() - start;
  long double avg_lcp = (long double)lcp_sum / text_length;
  fprintf(stderr, "\n\nComputation finished. Summary:\n");
  fprintf(stderr, "  Total time = %.2Lfs\n", total_time);
  fprintf(stderr, "  Relative time = %.2Lfus/input symbol\n",
      (1000000.0 * total_time) / text_length);
  fprintf(stderr, "  I/O volume = %lu bytes (%.2Lfbytes/input symbol)\n",
      total_io_volume, (1.L * total_io_volume) / text_length);

#ifdef MONITOR_DISK_USAGE
  fprintf(stderr, "  Internal I/O volume counter = %lu\n",
      utils::get_current_io_volume());
#endif

  fprintf(stderr, "  RAM allocation: cur = %lu bytes, peak = %.2LfMiB\n",
      utils::get_current_ram_allocation(),
      (1.L * utils::get_peak_ram_allocation()) / (1UL << 20));

#ifdef MONITOR_DISK_USAGE
  fprintf(stderr, "  Disk allocation: cur = %.2LfGiB, peak = %.2LfGiB\n",
      (1.L * utils::get_current_disk_allocation()) / (1UL << 30),
      (1.L * utils::get_peak_disk_allocation()) / (1UL << 30));
#endif

  fprintf(stderr, "  Sum of LCP = %lu\n", lcp_sum);  // XXX 128-bit integers?
  fprintf(stderr, "  Average LCP = %.2Lf\n", avg_lcp);
  fprintf(stderr, "  Maximal LCP = %lu\n", max_lcp);
}

}  // namespace em_sparse_phi_private

template<typename char_type,
  typename text_offset_type>
void em_sparse_phi(
    std::string text_filename,
    std::string sa_filename,
    std::string output_filename,
    std::uint64_t ram_use,
    bool inplace_mode) {
  em_sparse_phi_private::em_sparse_phi<char_type, text_offset_type>(
      text_filename, sa_filename, output_filename, ram_use, inplace_mode);
}

#endif  // __SRC_EM_SPARSE_PHI_SRC_EM_SPARSE_PHI_HPP_INCLUDED
