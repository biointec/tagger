/**
 * @file    src/em_sparse_phi_src/compute_lcp_delta.hpp
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

#ifndef __SRC_EM_SPARSE_PHI_SRC_COMPUTE_LCP_DELTA_HPP_INCLUDED
#define __SRC_EM_SPARSE_PHI_SRC_COMPUTE_LCP_DELTA_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <algorithm>
#include <parallel/algorithm>
#include <omp.h>

#include "io/async_stream_reader.hpp"
#include "io/async_stream_writer.hpp"
#include "utils.hpp"
#include "convert_to_vbyte_slab.hpp"


namespace em_sparse_phi_private {

struct buf_item_3 {
  buf_item_3() {}
  buf_item_3(
      std::uint64_t left_pos,
      std::uint64_t right_pos,
      std::uint64_t max_lcp_delta,
      std::uint64_t orig_id)
    : m_left_pos(left_pos),
    m_right_pos(right_pos),
    m_max_lcp_delta(max_lcp_delta),
    m_orig_id(orig_id) {}

  bool operator < (const buf_item_3 &x) const {
    std::int64_t my_diff =
      (std::int64_t)m_right_pos -
      (std::int64_t)m_left_pos;
    std::int64_t his_diff =
      (std::int64_t)x.m_right_pos -
      (std::int64_t)x.m_left_pos;
    return (my_diff < his_diff) ||
      (my_diff == his_diff && m_left_pos < x.m_left_pos);
  }

  std::uint64_t m_left_pos;
  std::uint64_t m_right_pos;
  std::uint64_t m_max_lcp_delta;
  std::uint64_t m_orig_id;
  std::uint64_t m_lcp_delta;
};

// Compute lcp delta-values.
template<typename char_type,
  typename text_offset_type>
void compute_lcp_delta_lex_partitioning(
    std::string text_filename,
    std::string output_filename,
    std::uint64_t text_length,
    std::uint64_t max_overflow_size,
    std::uint64_t max_halfsegment_size,
    bool first_sa_part,
    std::uint64_t in_buf_ram,
    std::uint64_t out_buf_ram,
    std::uint64_t local_buf_ram,
    std::string **delta_filenames,
    std::uint64_t &total_io_volume) {

  // Print initial message and start the timer.
  fprintf(stderr, "  Process halfsegments:\n");
  long double start = utils::wclock();

  // Compute the number of halfsegments.
  std::uint64_t n_halfsegments =
    (text_length + max_halfsegment_size - 1) / max_halfsegment_size;

  // Load every possible pair of half-segments and compute
  // all plcp values in that half-segments by brute force.
  std::uint64_t max_ext_halfsegment_size =
    max_halfsegment_size + max_overflow_size;
  char_type *left_halfsegment =
    utils::allocate_array<char_type>(max_ext_halfsegment_size);
  char_type *right_halfsegment =
    utils::allocate_array<char_type>(max_ext_halfsegment_size);

  // Scan all halfsegments left-to-right.
  for (std::uint64_t left_halfsegment_id = 0;
      left_halfsegment_id < n_halfsegments; ++left_halfsegment_id) {
    std::uint64_t left_halfsegment_beg =
      left_halfsegment_id * max_halfsegment_size;
    std::uint64_t left_halfsegment_end =
      std::min(left_halfsegment_beg + max_halfsegment_size, text_length);
    std::uint64_t left_halfsegment_ext_end =
      std::min(left_halfsegment_end + max_overflow_size, text_length);
    std::uint64_t left_halfsegment_ext_size =
      left_halfsegment_ext_end - left_halfsegment_beg;
    bool left_halfsegment_loaded = false;

    // Scan all halfsegments to the right of left_halfsegment_id.
    for (std::uint64_t right_halfsegment_id = left_halfsegment_id;
        right_halfsegment_id < n_halfsegments; right_halfsegment_id++) {
      std::uint64_t right_halfsegment_beg =
        right_halfsegment_id * max_halfsegment_size;
      std::uint64_t right_halfsegment_end =
        std::min(right_halfsegment_beg + max_halfsegment_size, text_length);
      std::uint64_t right_halfsegment_ext_end =
        std::min(right_halfsegment_end + max_overflow_size, text_length);
      std::uint64_t right_halfsegment_ext_size =
        right_halfsegment_ext_end - right_halfsegment_beg;

      // Check if that pairs of halfsegments has any associated pairs.
      std::string pairs_filename = output_filename + ".pairs." +
        utils::intToStr(left_halfsegment_id) + "_" +
        utils::intToStr(right_halfsegment_id);
      if (utils::file_exists(pairs_filename) == false ||
          utils::file_size(pairs_filename) == 0) {

        // Delete empty file.
        if (utils::file_exists(pairs_filename))
          utils::file_delete(pairs_filename);
        continue;
      }

      // Print initial message.
      fprintf(stderr, "    Process halfsegments %lu and %lu: ",
          left_halfsegment_id, right_halfsegment_id);

      // Initialize the timer.
      long double halfsegment_process_start = utils::wclock();

      // Initialize stats.
      std::uint64_t lcp_delta_sum = 0;
      std::uint64_t extra_io = 0;

      // Initialize reading from file associated
      // with current pair of halfsegments.
      typedef async_stream_reader<text_offset_type> pair_reader_type;
      std::uint64_t n_pairs =
        utils::file_size(pairs_filename) / (2 * sizeof(text_offset_type));
      pair_reader_type *pair_reader = new pair_reader_type(pairs_filename,
          in_buf_ram, std::max(4UL, in_buf_ram / (2UL << 20)));

      // Read left halfsegment from disk (if it wasn't already)
      if (left_halfsegment_loaded == false) {
        std::uint64_t offset = left_halfsegment_beg * sizeof(char_type);
        utils::read_at_offset(left_halfsegment, offset,
            left_halfsegment_ext_size, text_filename);
        left_halfsegment_loaded = true;
        extra_io += left_halfsegment_ext_size * sizeof(char_type);
      }

      // Read right halfsegment from disk.
      char_type *right_halfsegment_ptr = right_halfsegment;
      if (right_halfsegment_id != left_halfsegment_id) {
        std::uint64_t offset = right_halfsegment_beg * sizeof(char_type);
        utils::read_at_offset(right_halfsegment, offset,
            right_halfsegment_ext_size, text_filename);
        extra_io += right_halfsegment_ext_size * sizeof(char_type);
      } else right_halfsegment_ptr = left_halfsegment;

      // Initialize writing of the lcp values
      // computed when processing current halfsegments.
      typedef async_stream_writer<std::uint8_t> lcp_writer_type;
      std::string lcp_values_filename =
        delta_filenames[left_halfsegment_id][right_halfsegment_id];
      lcp_writer_type *lcp_writer = new lcp_writer_type(lcp_values_filename,
          out_buf_ram, std::max(4UL, out_buf_ram / (2UL << 20)),
          first_sa_part ? "w" : "a");

      static const std::uint64_t long_lcp_threshold = 50000;
      static const std::uint64_t very_long_lcp_threshold = 200000;
      static const std::uint64_t long_lcp_count_threshold = 20;

#ifdef _OPENMP
      std::uint64_t local_buf_size = local_buf_ram /
        (2 * sizeof(text_offset_type) +
         9 * sizeof(std::uint8_t) +
         1 * sizeof(std::uint64_t) +
         2 * sizeof(buf_item_3));

      std::uint64_t local_buf_filled = 0;
      text_offset_type *pair_buffer =
        utils::allocate_array<text_offset_type>(local_buf_size * 2);
      std::uint64_t *lcp_deltas =
        utils::allocate_array<std::uint64_t>(local_buf_size);
      std::uint8_t *vbyte_slab =
        utils::allocate_array<std::uint8_t>(local_buf_size * 9);
      buf_item_3 *item_buffer =
        utils::allocate_array<buf_item_3>(local_buf_size);

      std::uint64_t pairs_read = 0;
      while (pairs_read < n_pairs) {
        local_buf_filled = std::min(n_pairs - pairs_read, local_buf_size);
        pair_reader->read(pair_buffer, local_buf_filled * 2);

        bool long_lcp_mode = false;
        std::uint64_t buf_lcp_delta_sum = 0;
        std::uint64_t long_lcps_spotted = 0;
        std::uint64_t max_threads = omp_get_max_threads();
        std::uint64_t max_block_size =
          (local_buf_filled + max_threads - 1) / max_threads;
        std::uint64_t n_blocks =
          (local_buf_filled + max_block_size - 1) / max_block_size;

        #pragma omp parallel num_threads(n_blocks)
        {

          // Initialize private variables.
          std::uint64_t local_lcp_delta_sum = 0;
          std::uint64_t block_id = omp_get_thread_num();
          std::uint64_t block_beg = block_id * max_block_size;
          std::uint64_t block_end =
            std::min(block_beg + max_block_size, local_buf_filled);

          for (std::uint64_t j = block_beg; j < block_end; ++j) {
            std::uint64_t left_pos = pair_buffer[2 * j];
            std::uint64_t right_pos = pair_buffer[2 * j + 1];
            std::uint64_t lcp_delta = 0;
            std::uint64_t max_lcp_delta = text_length;

            if (long_lcp_mode == true)
              break;

            // Check if the current pair was the last one in the
            // multipair chain spawned by a single pair (i, Phi[i])).
            // If yes, one of the integers encodes the upper bound
            // on max_lcp_delta.
            std::uint64_t msb_bit =
              (1UL << (8 * sizeof(text_offset_type) - 1));
            if (left_pos & msb_bit) {
              if (right_pos & msb_bit) {
                max_lcp_delta = left_pos - msb_bit;
                right_pos  = right_pos - msb_bit;
                left_pos = max_overflow_size;
              } else {
                max_lcp_delta = right_pos;
                left_pos = left_pos - msb_bit;
                right_pos = max_overflow_size;
              }
            }

            // Compute the LCP value naively.
            while (lcp_delta < max_lcp_delta &&
                left_pos + lcp_delta < left_halfsegment_ext_size &&
                right_pos + lcp_delta < right_halfsegment_ext_size &&
                left_halfsegment[left_pos + lcp_delta] ==
                right_halfsegment_ptr[right_pos + lcp_delta])
              ++lcp_delta;

            lcp_deltas[j] = lcp_delta;
            local_lcp_delta_sum += lcp_delta;

            if (lcp_delta >= long_lcp_threshold) ++long_lcps_spotted;
            if (long_lcps_spotted >= long_lcp_count_threshold ||
                lcp_delta >= very_long_lcp_threshold)
              long_lcp_mode = true;
          }

          // Reduction of the private variables.
          #pragma omp critical
          {
            buf_lcp_delta_sum += local_lcp_delta_sum;
          }
        }

        if (long_lcp_mode == true) {
          buf_lcp_delta_sum = 0;

          // Long LCPs mode.
          #pragma omp parallel for
          for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
            std::uint64_t left_pos = pair_buffer[2 * j];
            std::uint64_t right_pos = pair_buffer[2 * j + 1];
            std::uint64_t max_lcp_delta = text_length;

            // Check if the current pair was the last one in
            // the chain spawned by a single pair (i, Phi[i])).
            // If yes, one of the integers encodes the upper
            // bound on max_lcp_delta.
            std::uint64_t msb_bit =
              (1UL << (8 * sizeof(text_offset_type) - 1));
            if (left_pos & msb_bit) {
              if (right_pos & msb_bit) {
                max_lcp_delta = left_pos - msb_bit;
                right_pos  = right_pos - msb_bit;
                left_pos = max_overflow_size;
              } else {
                max_lcp_delta = right_pos;
                left_pos = left_pos - msb_bit;
                right_pos = max_overflow_size;
              }
            }

            item_buffer[j] =
              buf_item_3(left_pos, right_pos, max_lcp_delta, j);
          }

          // Sort item buffer.
          __gnu_parallel::sort(item_buffer, item_buffer + local_buf_filled);

          // Process item buffer.
          #pragma omp parallel num_threads(n_blocks)
          {
            std::uint64_t block_id = omp_get_thread_num();
            std::uint64_t block_beg = block_id * max_block_size;
            std::uint64_t block_end =
              std::min(local_buf_filled, block_beg + max_block_size);
            std::uint64_t local_lcp_delta_sum = 0;

            std::int64_t prev_diff = 0;
            for (std::uint64_t j = block_beg; j < block_end; ++j) {
              std::uint64_t left_pos = item_buffer[j].m_left_pos;
              std::uint64_t right_pos = item_buffer[j].m_right_pos;
              std::int64_t cur_diff =
                (std::int64_t)right_pos -
                (std::int64_t)left_pos;
              std::uint64_t max_lcp_delta = item_buffer[j].m_max_lcp_delta;
              std::uint64_t lcp_delta = 0;

              if (j > block_beg && prev_diff == cur_diff &&
                  item_buffer[j - 1].m_left_pos +
                  item_buffer[j - 1].m_lcp_delta >= left_pos) {

                // Main trick of the long-lcp mode. Obtain the lcp from
                // the previously processed pair in the item buffer.
                lcp_delta = std::min(max_lcp_delta,
                    (item_buffer[j - 1].m_left_pos +
                     item_buffer[j - 1].m_lcp_delta) - left_pos);
              } else {

                // Compute the LCP value naively.
                while (lcp_delta < max_lcp_delta &&
                    left_pos + lcp_delta < left_halfsegment_ext_size &&
                    right_pos + lcp_delta < right_halfsegment_ext_size &&
                    left_halfsegment[left_pos + lcp_delta] ==
                    right_halfsegment_ptr[right_pos + lcp_delta])
                  ++lcp_delta;
              }

              item_buffer[j].m_lcp_delta = lcp_delta;
              local_lcp_delta_sum += lcp_delta;
              prev_diff = cur_diff;
            }

            #pragma omp critical
            {
              buf_lcp_delta_sum += local_lcp_delta_sum;
            }
          }

          // Permute LCP deltas back to the correct order.
          #pragma omp parallel for
          for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
            std::uint64_t addr = item_buffer[j].m_orig_id;
            std::uint64_t lcp_delta = item_buffer[j].m_lcp_delta;
            lcp_deltas[addr] = lcp_delta;
          }
        }

        // Update delta sum.
        lcp_delta_sum += buf_lcp_delta_sum;

        // Convert LCP values to vbyte encoding.
        std::uint64_t vbyte_slab_length =
          convert_to_vbyte_slab(lcp_deltas, local_buf_filled, vbyte_slab);

        // Write LCP values to disk.
        lcp_writer->write(vbyte_slab, vbyte_slab_length);
        pairs_read += local_buf_filled;
      }

      utils::deallocate(item_buffer);
      utils::deallocate(vbyte_slab);
      utils::deallocate(lcp_deltas);
      utils::deallocate(pair_buffer);

#else
      {
        std::uint64_t local_buf_size = local_buf_ram /
          (2 * sizeof(text_offset_type) +
           1 * sizeof(std::uint64_t) +
           1 * sizeof(buf_item_3));
        std::uint64_t local_buf_filled = 0;

        text_offset_type *pair_buffer =
          utils::allocate_array<text_offset_type>(local_buf_size * 2);
        std::uint64_t *lcp_deltas =
          utils::allocate_array<std::uint64_t>(local_buf_size);
        buf_item_3 *item_buffer =
          utils::allocate_array<buf_item_3>(local_buf_size);

        std::uint64_t pairs_read = 0;
        while (pairs_read < n_pairs) {
          local_buf_filled = std::min(n_pairs - pairs_read, local_buf_size);
          pair_reader->read(pair_buffer, local_buf_filled * 2);

          // Process the buffer using naive method.
          bool long_lcp_mode = false;
          std::uint64_t long_lcps_spotted = 0;
          for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
            std::uint64_t left_pos = pair_buffer[2 * j];
            std::uint64_t right_pos = pair_buffer[2 * j + 1];
            std::uint64_t lcp_delta = 0;
            std::uint64_t max_lcp_delta = text_length;

            // Check if the current pair was the last one in
            // the chain spawned by a single pair (i, Phi[i])).
            // If yes, one of the integers encodes the upper
            // bound on max_lcp_delta.
            std::uint64_t msb_bit =
              (1UL << (8 * sizeof(text_offset_type) - 1));
            if (left_pos & msb_bit) {
              if (right_pos & msb_bit) {
                max_lcp_delta = left_pos - msb_bit;
                right_pos  = right_pos - msb_bit;
                left_pos = max_overflow_size;
              } else {
                max_lcp_delta = right_pos;
                left_pos = left_pos - msb_bit;
                right_pos = max_overflow_size;
              }
            }

            // Compute the LCP value.
            while (lcp_delta < max_lcp_delta &&
                left_pos + lcp_delta < left_halfsegment_ext_size &&
                right_pos + lcp_delta < right_halfsegment_ext_size &&
                left_halfsegment[left_pos + lcp_delta] ==
                right_halfsegment_ptr[right_pos + lcp_delta])
              ++lcp_delta;

            lcp_deltas[j] = lcp_delta;
            if (lcp_delta >= long_lcp_threshold)
              ++long_lcps_spotted;

            if (long_lcps_spotted >= long_lcp_count_threshold ||
                lcp_delta >= very_long_lcp_threshold) {
              long_lcp_mode = true;
              break;
            }
          }

          if (long_lcp_mode == true) {

            // Long LCPs mode.
            for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
              std::uint64_t left_pos = pair_buffer[2 * j];
              std::uint64_t right_pos = pair_buffer[2 * j + 1];
              std::uint64_t max_lcp_delta = text_length;

              // Check if the current pair was the last one in
              // the chain spawned by a single pair (i, Phi[i])).
              // If yes, one of the integers encodes the upper
              // bound on max_lcp_delta.
              std::uint64_t msb_bit =
                (1UL << (8 * sizeof(text_offset_type) - 1));
              if (left_pos & msb_bit) {
                if (right_pos & msb_bit) {
                  max_lcp_delta = left_pos - msb_bit;
                  right_pos  = right_pos - msb_bit;
                  left_pos = max_overflow_size;
                } else {
                  max_lcp_delta = right_pos;
                  left_pos = left_pos - msb_bit;
                  right_pos = max_overflow_size;
                }
              }

              item_buffer[j] =
                buf_item_3(left_pos, right_pos, max_lcp_delta, j);
            }

            // Sort item buffer.
            std::sort(item_buffer, item_buffer + local_buf_filled);

            // Process item buffer.
            std::int64_t prev_diff = 0;
            for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
              std::uint64_t left_pos = item_buffer[j].m_left_pos;
              std::uint64_t right_pos = item_buffer[j].m_right_pos;
              std::int64_t cur_diff =
                (std::int64_t)right_pos -
                (std::int64_t)left_pos;
              std::uint64_t max_lcp_delta = item_buffer[j].m_max_lcp_delta;
              std::uint64_t lcp_delta = 0;

              if (j > 0 && prev_diff == cur_diff &&
                  item_buffer[j - 1].m_left_pos +
                  item_buffer[j - 1].m_lcp_delta >= left_pos) {

                // Main trick of the long-lcp mode. Obtain the lcp from
                // the previously processed pair in the item buffer.
                lcp_delta = std::min(max_lcp_delta,
                    (item_buffer[j - 1].m_left_pos +
                     item_buffer[j - 1].m_lcp_delta) - left_pos);
              } else {

                // Compute the LCP value naively.
                while (lcp_delta < max_lcp_delta &&
                    left_pos + lcp_delta < left_halfsegment_ext_size &&
                    right_pos + lcp_delta < right_halfsegment_ext_size &&
                    left_halfsegment[left_pos + lcp_delta] ==
                    right_halfsegment_ptr[right_pos + lcp_delta])
                  ++lcp_delta;
              }

              item_buffer[j].m_lcp_delta = lcp_delta;
              prev_diff = cur_diff;
            }

            // Permute LCP deltas back to the correct order.
            for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
              std::uint64_t addr = item_buffer[j].m_orig_id;
              std::uint64_t lcp_delta = item_buffer[j].m_lcp_delta;
              lcp_deltas[addr] = lcp_delta;
            }
          }

          // Write LCP deltas to file using v-byte encoding.
          for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
            std::uint64_t lcp_delta = lcp_deltas[j];
            lcp_delta_sum += lcp_delta;
            while (lcp_delta > 127) {
              std::uint64_t val = (lcp_delta & 0x7f) | 0x80;
              lcp_writer->write((std::uint8_t)val);
              lcp_delta >>= 7;
            }
            lcp_writer->write((std::uint8_t)lcp_delta);
          }

          pairs_read += local_buf_filled;
        }

        utils::deallocate(item_buffer);
        utils::deallocate(lcp_deltas);
        utils::deallocate(pair_buffer);
      }
#endif

      // Stop I/O threads.
      pair_reader->stop_reading();

      // Update I/O volume.
      std::uint64_t io_vol =
        lcp_writer->bytes_written() +
        pair_reader->bytes_read() +
        extra_io;
      total_io_volume += io_vol;

      // Print summary.
      long double avg_lcp_delta = (long double)lcp_delta_sum /
        (long double)std::max(1UL, n_pairs);
      long double elapsed = utils::wclock() - halfsegment_process_start;
      fprintf(stderr, "time = %.1Lfs, I/O = %.1LfMiB/s, "
          "avg. LCP delta = %.2Lf\n",
          elapsed, (1.L * io_vol / (1L << 20)) / elapsed, avg_lcp_delta);

      // Clean up.
      delete lcp_writer;
      delete pair_reader;
      utils::file_delete(pairs_filename);
    }
  }

  // Clean up.
  utils::deallocate(right_halfsegment);
  utils::deallocate(left_halfsegment);

  // Print summary.
  long double total_time = utils::wclock() - start;
  fprintf(stderr, "    Summary: time = %.2Lfs, "
      "total I/O vol = %.2Lfbytes/input symbol\n",
      total_time, (1.L * total_io_volume) / text_length);
}

// Compute LCP delta values.
template<typename char_type,
  typename text_offset_type>
void compute_lcp_delta_text_partitioning(
    std::string text_filename,
    std::string output_filename,
    std::uint64_t text_length,
    std::uint64_t max_overflow_size,
    std::uint64_t max_halfsegment_size,
    std::uint64_t in_buf_ram,
    std::uint64_t out_buf_ram,
    std::uint64_t local_buf_ram,
    std::string ***delta_filenames,
    std::uint64_t part_id,
    std::uint64_t &total_io_volume) {

  // Print initial message.
  fprintf(stderr, "  Process halfsegments:\n");

  // Initialize the timer.
  long double start = utils::wclock();

  // Open text file.
  std::FILE *f_text = utils::file_open_nobuf(text_filename, "r");

  // Compute basic stats.
  std::uint64_t n_halfsegments =
    (text_length + max_halfsegment_size - 1) / max_halfsegment_size;

  // Allocate halfsegments.
  std::uint64_t max_ext_halfsegment_size =
    max_halfsegment_size + max_overflow_size;
  char_type *left_halfsegment =
    utils::allocate_array<char_type>(max_ext_halfsegment_size);
  char_type *right_halfsegment =
    utils::allocate_array<char_type>(max_ext_halfsegment_size);

  // Load every possible pair of halfsegments and compute
  // all PLCP values in that halfsegments by brute force.
  for (std::uint64_t left_halfsegment_id = 0;
      left_halfsegment_id < n_halfsegments; ++left_halfsegment_id) {
    std::uint64_t left_halfsegment_beg =
      left_halfsegment_id * max_halfsegment_size;
    std::uint64_t left_halfsegment_end =
      std::min(left_halfsegment_beg + max_halfsegment_size, text_length);
    std::uint64_t left_halfsegment_ext_end =
      std::min(left_halfsegment_end + max_overflow_size, text_length);
    std::uint64_t left_halfsegment_ext_size =
      left_halfsegment_ext_end - left_halfsegment_beg;
    bool left_halfsegment_loaded = false;

    // Scan all halfsegments to the right of left_halfsegment_id.
    for (std::uint64_t right_halfsegment_id = left_halfsegment_id;
        right_halfsegment_id < n_halfsegments; right_halfsegment_id++) {
      std::uint64_t right_halfsegment_beg =
        right_halfsegment_id * max_halfsegment_size;
      std::uint64_t right_halfsegment_end =
        std::min(right_halfsegment_beg + max_halfsegment_size, text_length);
      std::uint64_t right_halfsegment_ext_end =
        std::min(right_halfsegment_end + max_overflow_size, text_length);
      std::uint64_t right_halfsegment_ext_size =
        right_halfsegment_ext_end - right_halfsegment_beg;

      // Check if that pairs of halfsegments has any associated pairs.
      std::string pairs_filename = output_filename + ".pairs." +
        utils::intToStr(left_halfsegment_id) + "_" +
        utils::intToStr(right_halfsegment_id);
      if (utils::file_exists(pairs_filename) == false ||
          utils::file_size(pairs_filename) == 0) {
        if (utils::file_exists(pairs_filename))
          utils::file_delete(pairs_filename);
        continue;
      }

      // Print initial progress message.
      fprintf(stderr, "    Process halfsegments %lu and %lu: ",
          left_halfsegment_id, right_halfsegment_id);

      // Start the timer.
      long double halfsegment_process_start = utils::wclock();

      // Initialize stats.
      std::uint64_t lcp_delta_sum = 0;
      std::uint64_t extra_io = 0;

      // Initialize reading from file associated
      // with current pair of halfsegments.
      typedef async_stream_reader<text_offset_type> pair_reader_type;
      std::uint64_t n_pairs =
        utils::file_size(pairs_filename) / (2 * sizeof(text_offset_type));
      pair_reader_type *pair_reader = new pair_reader_type(pairs_filename,
          in_buf_ram, std::max(4UL, in_buf_ram / (2UL << 20)));

      // Read left halfsegment from disk (if it wasn't already)
      if (left_halfsegment_loaded == false) {
        std::uint64_t offset = left_halfsegment_beg * sizeof(char_type);
        utils::read_at_offset(left_halfsegment, offset,
            left_halfsegment_ext_size, text_filename);
        left_halfsegment_loaded = true;
        extra_io += left_halfsegment_ext_size * sizeof(char_type);
      }

      // Read right halfsegment from disk.
      char_type *right_halfsegment_ptr = right_halfsegment;
      if (right_halfsegment_id != left_halfsegment_id) {
        std::uint64_t offset = right_halfsegment_beg * sizeof(char_type);
        utils::read_at_offset(right_halfsegment, offset,
            right_halfsegment_ext_size, text_filename);
        extra_io += right_halfsegment_ext_size * sizeof(char_type);
      } else right_halfsegment_ptr = left_halfsegment;

      // Initialize writing of the LCP values computed
      // during the processing of current halfsegments.
      typedef async_stream_writer<std::uint8_t> lcp_writer_type;
      lcp_writer_type *lcp_writer = new lcp_writer_type(
          delta_filenames[part_id][left_halfsegment_id][right_halfsegment_id],
          out_buf_ram, std::max(4UL, out_buf_ram / (2UL << 20)), "w");

      static const std::uint64_t long_lcp_threshold = 50000;
      static const std::uint64_t very_long_lcp_threshold = 200000;
      static const std::uint64_t long_lcp_count_threshold = 20;

#ifdef _OPENMP
      std::uint64_t local_buf_size = local_buf_ram /
        (2 * sizeof(text_offset_type) +
         9 * sizeof(std::uint8_t) +
         1 * sizeof(std::uint64_t) +
         2 * sizeof(buf_item_3));  // because it is not sorted in place

      std::uint64_t local_buf_filled = 0;
      text_offset_type *pair_buffer =
        utils::allocate_array<text_offset_type>(local_buf_size * 2);
      std::uint64_t *lcp_deltas =
        utils::allocate_array<std::uint64_t>(local_buf_size);
      std::uint8_t *vbyte_slab =
        utils::allocate_array<std::uint8_t>(local_buf_size * 9);
      buf_item_3 *item_buffer =
        utils::allocate_array<buf_item_3>(local_buf_size);

      std::uint64_t pairs_read = 0;
      while (pairs_read < n_pairs) {
        local_buf_filled = std::min(n_pairs - pairs_read, local_buf_size);
        pair_reader->read(pair_buffer, local_buf_filled * 2);

        bool long_lcp_mode = false;
        std::uint64_t buf_lcp_delta_sum = 0;
        std::uint64_t long_lcps_spotted = 0;
        std::uint64_t max_threads = omp_get_max_threads();
        std::uint64_t max_block_size =
          (local_buf_filled + max_threads - 1) / max_threads;
        std::uint64_t n_blocks =
          (local_buf_filled + max_block_size - 1) / max_block_size;

        #pragma omp parallel num_threads(n_blocks)
        {

          // Initialize private variables.
          std::uint64_t local_lcp_delta_sum = 0;
          std::uint64_t block_id = omp_get_thread_num();
          std::uint64_t block_beg = block_id * max_block_size;
          std::uint64_t block_end =
            std::min(block_beg + max_block_size, local_buf_filled);

          for (std::uint64_t j = block_beg; j < block_end; ++j) {
            std::uint64_t left_pos = pair_buffer[2 * j];
            std::uint64_t right_pos = pair_buffer[2 * j + 1];
            std::uint64_t lcp_delta = 0;
            std::uint64_t max_lcp_delta = text_length;

            if (long_lcp_mode == true)
              break;

            // Check if the current pair was the last one in
            // the chain spawned by a single pair (i, Phi[i])).
            // If yes, one of the integers encodes the upper
            // bound on max_lcp_delta.
            std::uint64_t msb_bit =
              (1UL << (8 * sizeof(text_offset_type) - 1));
            if (left_pos & msb_bit) {
              if (right_pos & msb_bit) {
                max_lcp_delta = left_pos - msb_bit;
                right_pos  = right_pos - msb_bit;
                left_pos = max_overflow_size;
              } else {
                max_lcp_delta = right_pos;
                left_pos = left_pos - msb_bit;
                right_pos = max_overflow_size;
              }
            }

            // Compute the LCP value naively.
            while (lcp_delta < max_lcp_delta &&
                left_pos + lcp_delta < left_halfsegment_ext_size &&
                right_pos + lcp_delta < right_halfsegment_ext_size &&
                left_halfsegment[left_pos + lcp_delta] ==
                right_halfsegment_ptr[right_pos + lcp_delta])
              ++lcp_delta;

            lcp_deltas[j] = lcp_delta;
            local_lcp_delta_sum += lcp_delta;

            if (lcp_delta >= long_lcp_threshold) ++long_lcps_spotted;
            if (long_lcps_spotted >= long_lcp_count_threshold ||
                lcp_delta >= very_long_lcp_threshold)
              long_lcp_mode = true;
          }

          // Reduction of the private variables.
          #pragma omp critical
          {
            buf_lcp_delta_sum += local_lcp_delta_sum;
          }
        }

        if (long_lcp_mode == true) {
          buf_lcp_delta_sum = 0;

          // Long LCPs mode.
          #pragma omp parallel for
          for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
            std::uint64_t left_pos = pair_buffer[2 * j];
            std::uint64_t right_pos = pair_buffer[2 * j + 1];
            std::uint64_t max_lcp_delta = text_length;

            // Check if the current pair was the last one in
            // the chain spawned by a single pair (i, Phi[i])).
            // If yes, one of the integers encodes the upper
            // bound on max_lcp_delta.
            std::uint64_t msb_bit =
              (1UL << (8 * sizeof(text_offset_type) - 1));
            if (left_pos & msb_bit) {
              if (right_pos & msb_bit) {
                max_lcp_delta = left_pos - msb_bit;
                right_pos  = right_pos - msb_bit;
                left_pos = max_overflow_size;
              } else {
                max_lcp_delta = right_pos;
                left_pos = left_pos - msb_bit;
                right_pos = max_overflow_size;
              }
            }

            item_buffer[j] =
              buf_item_3(left_pos, right_pos, max_lcp_delta, j);
          }

          // Sort item buffer.
          __gnu_parallel::sort(item_buffer, item_buffer + local_buf_filled);

          // Process item buffer.
          #pragma omp parallel num_threads(n_blocks)
          {
            std::uint64_t block_id = omp_get_thread_num();
            std::uint64_t block_beg = block_id * max_block_size;
            std::uint64_t block_end =
              std::min(local_buf_filled, block_beg + max_block_size);
            std::uint64_t local_lcp_delta_sum = 0;

            std::int64_t prev_diff = 0;
            for (std::uint64_t j = block_beg; j < block_end; ++j) {
              std::uint64_t left_pos = item_buffer[j].m_left_pos;
              std::uint64_t right_pos = item_buffer[j].m_right_pos;
              std::int64_t cur_diff =
                (std::int64_t)right_pos -
                (std::int64_t)left_pos;
              std::uint64_t max_lcp_delta = item_buffer[j].m_max_lcp_delta;
              std::uint64_t lcp_delta = 0;

              if (j > block_beg && prev_diff == cur_diff &&
                  item_buffer[j - 1].m_left_pos +
                  item_buffer[j - 1].m_lcp_delta >= left_pos) {

                // Main trick of the long LCP mode. Obtain the LCP from
                // the previously processed pair in the item buffer.
                lcp_delta = std::min(max_lcp_delta,
                    (item_buffer[j - 1].m_left_pos +
                     item_buffer[j - 1].m_lcp_delta) - left_pos);
              } else {

                // Compute the LCP value naively.
                while (lcp_delta < max_lcp_delta &&
                    left_pos + lcp_delta < left_halfsegment_ext_size &&
                    right_pos + lcp_delta < right_halfsegment_ext_size &&
                    left_halfsegment[left_pos + lcp_delta] ==
                    right_halfsegment_ptr[right_pos + lcp_delta])
                  ++lcp_delta;
              }

              item_buffer[j].m_lcp_delta = lcp_delta;
              local_lcp_delta_sum += lcp_delta;
              prev_diff = cur_diff;
            }

            #pragma omp critical
            {
              buf_lcp_delta_sum += local_lcp_delta_sum;
            }
          }

          // Permute LCP deltas back to the correct order.
          #pragma omp parallel for
          for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
            std::uint64_t addr = item_buffer[j].m_orig_id;
            std::uint64_t lcp_delta = item_buffer[j].m_lcp_delta;
            lcp_deltas[addr] = lcp_delta;
          }
        }

        // Update delta sum.
        lcp_delta_sum += buf_lcp_delta_sum;

        // Convert LCP values to vbyte encoding.
        std::uint64_t vbyte_slab_length =
          convert_to_vbyte_slab(lcp_deltas, local_buf_filled, vbyte_slab);

        // Write LCP values to disk.
        lcp_writer->write(vbyte_slab, vbyte_slab_length);
        pairs_read += local_buf_filled;
      }

      // Clean up.
      utils::deallocate(item_buffer);
      utils::deallocate(vbyte_slab);
      utils::deallocate(lcp_deltas);
      utils::deallocate(pair_buffer);

#else
      {
        std::uint64_t local_buf_size = local_buf_ram /
          (2 * sizeof(text_offset_type) +
           1 * sizeof(std::uint64_t) +
           1 * sizeof(buf_item_3));

        std::uint64_t local_buf_filled = 0;
        text_offset_type *pair_buffer =
          utils::allocate_array<text_offset_type>(local_buf_size * 2);
        std::uint64_t *lcp_deltas =
          utils::allocate_array<std::uint64_t>(local_buf_size);
        buf_item_3 *item_buffer =
          utils::allocate_array<buf_item_3>(local_buf_size);

        std::uint64_t pairs_read = 0;
        while (pairs_read < n_pairs) {
          local_buf_filled = std::min(n_pairs - pairs_read, local_buf_size);
          pair_reader->read(pair_buffer, local_buf_filled * 2);

          // Process the buffer using naive method.
          bool long_lcp_mode = false;
          std::uint64_t long_lcps_spotted = 0;
          for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
            std::uint64_t left_pos = pair_buffer[2 * j];
            std::uint64_t right_pos = pair_buffer[2 * j + 1];
            std::uint64_t lcp_delta = 0;
            std::uint64_t max_lcp_delta = text_length;

            // Check if the current pair was the last one in
            // the chain spawned by a single pair (i, Phi[i])).
            // If yes, one of the integers encodes the upper
            // bound on max_lcp_delta.
            std::uint64_t msb_bit =
              (1UL << (8 * sizeof(text_offset_type) - 1));
            if (left_pos & msb_bit) {
              if (right_pos & msb_bit) {
                max_lcp_delta = left_pos - msb_bit;
                right_pos  = right_pos - msb_bit;
                left_pos = max_overflow_size;
              } else {
                max_lcp_delta = right_pos;
                left_pos = left_pos - msb_bit;
                right_pos = max_overflow_size;
              }
            }

            // Compute the LCP value.
            while (lcp_delta < max_lcp_delta &&
                left_pos + lcp_delta < left_halfsegment_ext_size &&
                right_pos + lcp_delta < right_halfsegment_ext_size &&
                left_halfsegment[left_pos + lcp_delta] ==
                right_halfsegment_ptr[right_pos + lcp_delta])
              ++lcp_delta;

            lcp_deltas[j] = lcp_delta;
            if (lcp_delta >= long_lcp_threshold)
              ++long_lcps_spotted;

            if (long_lcps_spotted >= long_lcp_count_threshold ||
                lcp_delta >= very_long_lcp_threshold) {
              long_lcp_mode = true;
              break;
            }
          }

          if (long_lcp_mode == true) {

            // Long LCP mode.
            for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
              std::uint64_t left_pos = pair_buffer[2 * j];
              std::uint64_t right_pos = pair_buffer[2 * j + 1];
              std::uint64_t max_lcp_delta = text_length;

              // Check if the current pair was the last one in
              // the chain spawned by a single pair (i, Phi[i])).
              // If yes, one of the integers encodes the upper
              // bound on max_lcp_delta.
              std::uint64_t msb_bit =
                (1UL << (8 * sizeof(text_offset_type) - 1));
              if (left_pos & msb_bit) {
                if (right_pos & msb_bit) {
                  max_lcp_delta = left_pos - msb_bit;
                  right_pos  = right_pos - msb_bit;
                  left_pos = max_overflow_size;
                } else {
                  max_lcp_delta = right_pos;
                  left_pos = left_pos - msb_bit;
                  right_pos = max_overflow_size;
                }
              }

              item_buffer[j] =
                buf_item_3(left_pos, right_pos, max_lcp_delta, j);
            }

            // Sort item buffer.
            std::sort(item_buffer, item_buffer + local_buf_filled);

            // Process item buffer.
            std::int64_t prev_diff = 0;
            for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
              std::uint64_t left_pos = item_buffer[j].m_left_pos;
              std::uint64_t right_pos = item_buffer[j].m_right_pos;
              std::int64_t cur_diff =
                (std::int64_t)right_pos -
                (std::int64_t)left_pos;
              std::uint64_t max_lcp_delta = item_buffer[j].m_max_lcp_delta;
              std::uint64_t lcp_delta = 0;

              if (j > 0 && prev_diff == cur_diff &&
                  item_buffer[j - 1].m_left_pos +
                  item_buffer[j - 1].m_lcp_delta >= left_pos) {

                // Main trick of the long LCP mode. Obtain the LCP from
                // the previously processed pair in the item buffer.
                lcp_delta = std::min(max_lcp_delta,
                    (item_buffer[j - 1].m_left_pos +
                     item_buffer[j - 1].m_lcp_delta) - left_pos);
              } else {

                // Compute the LCP value naively.
                while (lcp_delta < max_lcp_delta &&
                    left_pos + lcp_delta < left_halfsegment_ext_size &&
                    right_pos + lcp_delta < right_halfsegment_ext_size &&
                    left_halfsegment[left_pos + lcp_delta] ==
                    right_halfsegment_ptr[right_pos + lcp_delta])
                  ++lcp_delta;
              }

              item_buffer[j].m_lcp_delta = lcp_delta;
              prev_diff = cur_diff;
            }

            // Permute LCP deltas back to the correct order.
            for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
              std::uint64_t addr = item_buffer[j].m_orig_id;
              std::uint64_t lcp_delta = item_buffer[j].m_lcp_delta;
              lcp_deltas[addr] = lcp_delta;
            }
          }

          // Write LCP deltas to file using v-byte encoding.
          for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
            std::uint64_t lcp_delta = lcp_deltas[j];
            lcp_delta_sum += lcp_delta;
            while (lcp_delta > 127) {
              std::uint64_t val = (lcp_delta & 0x7f) | 0x80;
              lcp_writer->write((std::uint8_t)val);
              lcp_delta >>= 7;
            }
            lcp_writer->write((std::uint8_t)lcp_delta);
          }

          pairs_read += local_buf_filled;
        }

        // Clean up.
        utils::deallocate(item_buffer);
        utils::deallocate(lcp_deltas);
        utils::deallocate(pair_buffer);
      }
#endif

      // Stop I/O threads.
      pair_reader->stop_reading();

      // Update I/O volume.
      std::uint64_t io_vol =
        lcp_writer->bytes_written() +
        pair_reader->bytes_read() +
        extra_io;
      total_io_volume += io_vol;

      // Print summary.
      long double avg_lcp_delta = (long double)lcp_delta_sum /
        (long double)std::max(1UL, n_pairs);
      long double elapsed = utils::wclock() - halfsegment_process_start;
      fprintf(stderr, "time = %.1Lfs, I/O = %.1LfMiB/s, "
          "avg. LCP delta = %.2Lf\n", elapsed,
          (1.L * io_vol / (1L << 20)) / elapsed, avg_lcp_delta);

      // Clean up.
      delete lcp_writer;
      delete pair_reader;
      utils::file_delete(pairs_filename);
    }
  }

  // Clean up.
  utils::deallocate(right_halfsegment);
  utils::deallocate(left_halfsegment);
  std::fclose(f_text);

  // Print summary.
  long double total_time = utils::wclock() - start;
  fprintf(stderr, "    Summary: time = %.2Lfs, "
      "total I/O vol = %.2Lfbytes/input symbol\n", total_time,
      (1.L * total_io_volume) / text_length);
}

}  // namespace em_sparse_phi_private

#endif  // __SRC_EM_SPARSE_PHI_SRC_COMPUTE_LCP_DELTA_HPP_INCLUDED
