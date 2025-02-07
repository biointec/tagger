/**
 * @file    src/em_sparse_phi_src/compute_lcp_lower_bounds.hpp
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

#ifndef __SRC_EM_SPARSE_PHI_SRC_COMPUTE_LCP_LOWER_BOUNDS_HPP_INCLUDED
#define __SRC_EM_SPARSE_PHI_SRC_COMPUTE_LCP_LOWER_BOUNDS_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <algorithm>
#include <omp.h>

#include "io/async_stream_reader.hpp"
#include "io/async_multi_stream_writer.hpp"
#include "utils.hpp"


namespace em_sparse_phi_private {

// Bijectively map a pair of ints (i, j), where (0 <= i <= j < n)
// into a set of integers [0..k), where k = (n * (n + 1)) / 2.
inline std::uint64_t pair_to_int(std::uint64_t i,
    std::uint64_t j, std::uint64_t n) {
  return (j - i) + i * (n - i + 1) + (i * (i - 1)) / 2;
}

template<typename T, typename S>
struct buf_item {
  buf_item() {}
  buf_item(T ad, S v1, S v2, std::uint8_t sp)
    : addr(ad), val1(v1), val2(v2), special(sp) {}
  T addr;
  S val1;
  S val2;
  std::uint8_t special;
} __attribute__ ((packed));

// Extended buf item. Used only for text partitioning.
template<typename T, typename S, typename U>
struct ext_buf_item {
  ext_buf_item() {}
  ext_buf_item(T ad, U i1, U i2, S v1, S v2, std::uint8_t sp)
    : addr(ad), id1(i1), id2(i2), val1(v1), val2(v2), special(sp) {}
  T addr;
  U id1;
  U id2;
  S val1;
  S val2;
  std::uint8_t special;
} __attribute__ ((packed));

// Scan SA and for every pair (i, Phi[i]) use PLCP array to find
// the lower bound `lb' on PLCP[i]. Then add pair (i + lb, Phi[i] + lb)
// to the file corresponding to the pairs of segments containing positions
// i + lb and Phi[i] + lb. Large values of `lb' may spawn many pairs.
// A version of the function using lex-partitioning.
template<typename text_offset_type>
void compute_lcp_lower_bounds_lex_partitioning(
    std::string sa_filename,
    std::string output_filename,
    std::uint64_t text_length,
    std::uint64_t ram_use,
    text_offset_type *sparse_plcp,
    std::uint64_t plcp_sampling_rate,
    std::uint64_t max_halfsegment_size,
    std::uint64_t sa_part_beg,
    std::uint64_t sa_part_end,
    std::uint64_t halfseg_buffers_ram,
    std::uint64_t in_buf_ram,
    std::uint64_t local_buf_ram,
    std::uint64_t max_overflow_size,
    std::uint64_t &total_io_volume) {

  std::uint64_t n_halfsegments =
    (text_length + max_halfsegment_size - 1) / max_halfsegment_size;
  std::uint64_t sparse_plcp_size =
    (text_length + plcp_sampling_rate - 1) / plcp_sampling_rate;
  std::uint64_t sparse_plcp_ram = sparse_plcp_size * sizeof(text_offset_type);
  std::uint64_t sa_part_size = sa_part_end - sa_part_beg;

  // Print initial message.
  fprintf(stderr, "  Compute and store LCP lower bounds: ");

  // Start the timer.
  long double start = utils::wclock();

  // Initialize I/O volume.
  std::uint64_t io_volume = 0;

  // Proportionally enlarge the buffers to fit the available RAM.
  {
    std::uint64_t ram_budget = ram_use -
      sparse_plcp_ram - halfseg_buffers_ram;
    std::uint64_t total_buf_ram = in_buf_ram + local_buf_ram;
    long double scale_factor =
      (long double)ram_budget / (long double)total_buf_ram;
    in_buf_ram = (std::uint64_t)((long double)in_buf_ram * scale_factor);
    local_buf_ram =
      (std::uint64_t)((long double)local_buf_ram * scale_factor);
  }

  // Create SA reader.
  typedef async_stream_reader<text_offset_type> sa_reader_type;
  sa_reader_type *sa_reader = new sa_reader_type(sa_filename,
      in_buf_ram, std::max(4UL, in_buf_ram / (2UL << 20)), sa_part_beg);

  // Create a writer for every possible pair of halfsegments.
  static const std::uint64_t n_free_buffers = 4;
  std::uint64_t n_different_halfseg_pairs =
    (n_halfsegments * (n_halfsegments + 1)) / 2;
  std::uint64_t buffer_size = halfseg_buffers_ram /
    (n_different_halfseg_pairs + n_free_buffers);
  typedef async_multi_stream_writer<text_offset_type> pair_multiwriter_type;
  pair_multiwriter_type *pair_multiwriter = new pair_multiwriter_type(
      n_different_halfseg_pairs, buffer_size, n_free_buffers);
  for (std::uint64_t i = 0; i < n_halfsegments; ++i)
    for (std::uint64_t j = i; j < n_halfsegments; ++j)
      pair_multiwriter->add_file(output_filename + ".pairs." +
          utils::intToStr(i) + "_" + utils::intToStr(j));

  // Allocate buffers.
#ifdef _OPENMP
  typedef buf_item<std::uint64_t, std::uint64_t> buf_item_type;
  std::uint64_t local_buf_size = local_buf_ram /
    (1 * sizeof(buf_item_type) +
     1 * sizeof(text_offset_type) +
     3 * sizeof(std::uint64_t));
#else
  std::uint64_t local_buf_size = local_buf_ram /
    (1 * sizeof(text_offset_type) +
     3 * sizeof(std::uint64_t));
#endif

#ifdef _OPENMP
  buf_item_type *item_buffer =
    utils::allocate_array<buf_item_type>(local_buf_size);
#endif

  text_offset_type *sa_buffer =
    utils::allocate_array<text_offset_type>(local_buf_size);
  std::uint64_t *addr_buffer =
    utils::allocate_array<std::uint64_t>(local_buf_size);
  std::uint64_t *plcp_buffer =
    utils::allocate_array<std::uint64_t>(2 * local_buf_size);

  // Read SA values, one buffer at a time.
  std::uint64_t sa_items_read = 0;
//  std::uint64_t nontrivial_pairs = 0;
  std::uint64_t step_counter = 0;
  text_offset_type prev_sa = text_length;
  if (sa_part_beg > 0) {
    std::uint64_t offset = (sa_part_beg - 1) * sizeof(text_offset_type);
    utils::read_at_offset(&prev_sa, offset, 1, sa_filename);
    io_volume += sizeof(text_offset_type);
  }

  while (sa_items_read < sa_part_size) {
    ++step_counter;
    if (step_counter == (1UL << 25) / local_buf_size) {
      step_counter = 0;
      long double elapsed = utils::wclock() - start;
      std::uint64_t io_vol =
        sa_reader->bytes_read() +
        pair_multiwriter->bytes_written();
      fprintf(stderr, "\r  Compute and store LCP lower bounds: %.1Lf%%, "
          "time = %.1Lfs, I/O = %.2LfMiB/s",
          (100.L * sa_items_read) / sa_part_size, elapsed,
          (1.L * io_vol / (1L << 20)) / elapsed);
    }

    // Fill in the buffer with SA values.
    std::uint64_t local_buf_filled =
      std::min(sa_part_size - sa_items_read, local_buf_size);
    sa_reader->read(sa_buffer, local_buf_filled);

    // Process items in the buffer.
    // First, fetch PLCP values.
#ifdef _OPENMP
#ifdef PRODUCTION_VER
    #pragma omp parallel for
    for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
      std::uint64_t addr = sa_buffer[j] / plcp_sampling_rate;
      addr_buffer[j] = addr;
    }

    #pragma omp parallel for
    for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
      std::uint64_t addr = addr_buffer[j];
      plcp_buffer[2 * j] = sparse_plcp[addr];
      if (addr + 1 < sparse_plcp_size)
        plcp_buffer[2 * j + 1] = sparse_plcp[addr + 1];
    }
#else
    #pragma omp parallel for
    for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
      std::uint64_t addr = sa_buffer[j] / plcp_sampling_rate;
      addr_buffer[j] = addr;
      plcp_buffer[2 * j] = sparse_plcp[addr];
      if (addr + 1 < sparse_plcp_size)
        plcp_buffer[2 * j + 1] = sparse_plcp[addr + 1];
    }
#endif
#else
    for (std::uint64_t j = 0; j < local_buf_filled; ++j)
      addr_buffer[j] = sa_buffer[j] / plcp_sampling_rate;
    for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
      std::uint64_t addr = addr_buffer[j];
      plcp_buffer[2 * j] = sparse_plcp[addr];
      if (addr + 1 < sparse_plcp_size)
        plcp_buffer[2 * j + 1] = sparse_plcp[addr + 1];
    }
#endif

#ifdef _OPENMP
    {
      std::uint64_t different_halfsegments_pairs =
        (n_halfsegments * (n_halfsegments + 1)) / 2;
      buf_item_type dummy_buf_item(different_halfsegments_pairs,
          text_offset_type(0), text_offset_type(0), false);

      #pragma omp parallel for
      for (std::uint64_t i = 0; i < local_buf_filled; ++i) {
        std::uint64_t prevsa = (i == 0) ?
          (std::uint64_t)prev_sa :
          (std::uint64_t)sa_buffer[i - 1];
        std::uint64_t currsa = sa_buffer[i];

        if (prevsa == text_length) {
          item_buffer[i] = dummy_buf_item;
          continue;
        }

        // Compute lower bound for PLCP[i].
        std::uint64_t sampled_plcp_ptr = addr_buffer[i];
        std::uint64_t sampled_plcp_idx =
          sampled_plcp_ptr * plcp_sampling_rate;
        std::uint64_t prev_sample_dist = currsa - sampled_plcp_idx;
        std::uint64_t next_sample_dist =
          plcp_sampling_rate - prev_sample_dist;
        std::uint64_t plcp_val = plcp_buffer[2 * i];
        std::uint64_t plcp_lower_bound = std::max(0L,
            (std::int64_t)plcp_val -
            (std::int64_t)prev_sample_dist);
        std::uint64_t pos1 = currsa;
        std::uint64_t pos2 = prevsa;
        if (pos1 > pos2) std::swap(pos1, pos2);
        pos1 += plcp_lower_bound;
        pos2 += plcp_lower_bound;
        if (pos2 >= text_length) {
          item_buffer[i] = dummy_buf_item;
          continue;
        }

        // Compute the maximal LCP delta in
        // case PLCP[currsa] was not sampled.
        std::uint64_t max_lcp_delta = 0;
        if (sampled_plcp_idx != currsa) {
          std::uint64_t max_lcp_delta2 = text_length -
            (std::max(currsa, prevsa) + plcp_lower_bound);
          if (sampled_plcp_idx + plcp_sampling_rate < text_length) {
            std::uint64_t plcp_upper_bound =
              plcp_buffer[2 * i + 1] + next_sample_dist;
            max_lcp_delta = std::min(max_lcp_delta2,
                plcp_upper_bound -
                plcp_lower_bound);
          } else max_lcp_delta = max_lcp_delta2;
        }

        if (max_lcp_delta > 0) {

          // Compute the segments containing the positions to compare and
          // determine if the pair needs to be split into multiple pairs.
          std::uint64_t seg_1_idx = pos1 / max_halfsegment_size;
          std::uint64_t seg_2_idx = pos2 / max_halfsegment_size;
          std::uint64_t seg_1_beg = seg_1_idx * max_halfsegment_size;
          std::uint64_t seg_2_beg = seg_2_idx * max_halfsegment_size;
          std::uint64_t seg_1_ext_end = std::min(text_length,
              seg_1_beg + max_halfsegment_size + max_overflow_size);
          std::uint64_t seg_2_ext_end = std::min(text_length,
              seg_2_beg + max_halfsegment_size + max_overflow_size);
          std::uint64_t seg_1_maxlcp = seg_1_ext_end - pos1;
          std::uint64_t seg_2_maxlcp = seg_2_ext_end - pos2;
          std::uint64_t seg_maxlcp = std::min(seg_1_maxlcp, seg_2_maxlcp);

          // Add the appropriate item to the buffer.
          if (max_lcp_delta > seg_maxlcp) {

            // Special item. Will be split into chain of pairs when writing.
            item_buffer[i] = buf_item_type(0,
                text_offset_type(0), text_offset_type(0), true);
          } else {

            // Normal pair. Will be written as a single pair.
            std::uint64_t file_id =
              pair_to_int(seg_1_idx, seg_2_idx, n_halfsegments);
            item_buffer[i] = buf_item_type(file_id,
                text_offset_type(pos1 - seg_1_beg),
                text_offset_type(pos2 - seg_2_beg),
                false);
          }
        } else item_buffer[i] = dummy_buf_item;
      }

      // Distribute the items from the pair buffer to files.
      for (std::uint64_t i = 0; i < local_buf_filled; ++i) {
        std::uint64_t addr = item_buffer[i].addr;
        if (addr != different_halfsegments_pairs) {
          if (item_buffer[i].special == false) {

            // Write a single lcp pair.
            pair_multiwriter->write_to_ith_file(addr, item_buffer[i].val1);
            pair_multiwriter->write_to_ith_file(addr, item_buffer[i].val2);
          } else {

            // Recompute and write the chain of pairs.
            std::uint64_t prevsa = (i == 0) ?
              (std::uint64_t)prev_sa :
              (std::uint64_t)sa_buffer[i - 1];
            std::uint64_t currsa = sa_buffer[i];

            // Compute lower bound for PLCP[i].
            std::uint64_t sampled_plcp_ptr = addr_buffer[i];
            std::uint64_t sampled_plcp_idx =
              sampled_plcp_ptr * plcp_sampling_rate;
            std::uint64_t prev_sample_dist = currsa - sampled_plcp_idx;
            std::uint64_t next_sample_dist =
              plcp_sampling_rate - prev_sample_dist;
            std::uint64_t plcp_val = plcp_buffer[2 * i];
            std::uint64_t plcp_lower_bound = std::max(0L,
                (std::int64_t)plcp_val -
                (std::int64_t)prev_sample_dist);
            std::uint64_t pos1 = currsa + plcp_lower_bound;
            std::uint64_t pos2 = prevsa + plcp_lower_bound;
            if (pos1 > pos2) std::swap(pos1, pos2);

            // Compute the maximal LCP delta in
            // case PLCP[currsa] was not sampled.
            std::uint64_t max_lcp_delta = 0;
            if (sampled_plcp_idx != currsa) {
              std::uint64_t max_lcp_delta2 = text_length -
                (std::max(currsa, prevsa) + plcp_lower_bound);
              if (sampled_plcp_idx + plcp_sampling_rate < text_length) {
                std::uint64_t plcp_upper_bound =
                  plcp_buffer[2 * i + 1] + next_sample_dist;
                max_lcp_delta = std::min(max_lcp_delta2,
                    plcp_upper_bound -
                    plcp_lower_bound);
              } else max_lcp_delta = max_lcp_delta2;
            }

            // Compute the segments containing the positions to compare and
            // determine if the pair needs to be split into multiple pairs.
            std::uint64_t seg_1_idx = pos1 / max_halfsegment_size;
            std::uint64_t seg_2_idx = pos2 / max_halfsegment_size;
            std::uint64_t seg_1_beg = seg_1_idx * max_halfsegment_size;
            std::uint64_t seg_2_beg = seg_2_idx * max_halfsegment_size;
            std::uint64_t seg_1_ext_end = std::min(text_length,
                seg_1_beg + max_halfsegment_size + max_overflow_size);
            std::uint64_t seg_2_ext_end = std::min(text_length,
                seg_2_beg + max_halfsegment_size + max_overflow_size);
            std::uint64_t seg_1_maxlcp = seg_1_ext_end - pos1;
            std::uint64_t seg_2_maxlcp = seg_2_ext_end - pos2;
            std::uint64_t seg_maxlcp = std::min(seg_1_maxlcp, seg_2_maxlcp);

            // Split the pair into multiple pairs.
            while (max_lcp_delta > seg_maxlcp) {
              std::uint64_t file_id =
                pair_to_int(seg_1_idx, seg_2_idx, n_halfsegments);
              pair_multiwriter->write_to_ith_file(file_id, pos1 - seg_1_beg);
              pair_multiwriter->write_to_ith_file(file_id, pos2 - seg_2_beg);
              pos1 += seg_maxlcp;
              pos2 += seg_maxlcp;
              max_lcp_delta -= seg_maxlcp;
              seg_1_idx = pos1 / max_halfsegment_size;
              seg_2_idx = pos2 / max_halfsegment_size;
              seg_1_beg = seg_1_idx * max_halfsegment_size;
              seg_2_beg = seg_2_idx * max_halfsegment_size;
              seg_1_ext_end = std::min(text_length,
                  seg_1_beg + max_halfsegment_size + max_overflow_size);
              seg_2_ext_end = std::min(text_length,
                  seg_2_beg + max_halfsegment_size + max_overflow_size);
              seg_1_maxlcp = seg_1_ext_end - pos1;
              seg_2_maxlcp = seg_2_ext_end - pos2;
              seg_maxlcp = std::min(seg_1_maxlcp, seg_2_maxlcp);
            }

            // Write the last pair in the chain. Use one of the
            // integers in the pair to encode the upper bound for LCP.
            std::uint64_t file_id =
              pair_to_int(seg_1_idx, seg_2_idx, n_halfsegments);
            std::uint64_t msb_bit =
              (1UL << (8 * sizeof(text_offset_type) - 1));
            if (pos1 == seg_1_beg + max_overflow_size) {

              // We don't have to encode pos1 - seg_1_beg. Instead, we
              // encode max_lcp_delta. Of course we also have to set
              // the MSB bit in both integers to recognize this case.
              pair_multiwriter->write_to_ith_file(
                  file_id, max_lcp_delta + msb_bit);
              pair_multiwriter->write_to_ith_file(
                  file_id, (pos2 - seg_2_beg) + msb_bit);
            } else {

              // Analogously, but we don't encode pos2 - seg_2_beg.
              pair_multiwriter->write_to_ith_file(
                  file_id, (pos1 - seg_1_beg) + msb_bit);
              pair_multiwriter->write_to_ith_file(
                  file_id, max_lcp_delta);
            }
          }
        }
      }
    }
#else
    {
      for (std::uint64_t i = 0; i < local_buf_filled; ++i) {
        std::uint64_t prevsa = (i == 0) ?
          (std::uint64_t)prev_sa :
          (std::uint64_t)sa_buffer[i - 1];
        std::uint64_t currsa = sa_buffer[i];

        if (prevsa == text_length)
          continue;

        // Compute lower bound for PLCP[i].
        std::uint64_t sampled_plcp_ptr = addr_buffer[i];
        std::uint64_t sampled_plcp_idx =
          sampled_plcp_ptr * plcp_sampling_rate;
        std::uint64_t prev_sample_dist = currsa - sampled_plcp_idx;
        std::uint64_t next_sample_dist =
          plcp_sampling_rate - prev_sample_dist;
        std::uint64_t plcp_val = plcp_buffer[2 * i];
        std::uint64_t plcp_lower_bound = std::max(0L,
            (std::int64_t)plcp_val -
            (std::int64_t)prev_sample_dist);
        std::uint64_t pos1 = currsa;
        std::uint64_t pos2 = prevsa;
        if (pos1 > pos2) std::swap(pos1, pos2);
        pos1 += plcp_lower_bound;
        pos2 += plcp_lower_bound;
        if (pos2 >= text_length)
          continue;

        // Compute the maximal LCP delta in
        // case PLCP[currsa] was not sampled.
        std::uint64_t max_lcp_delta = 0;
        if (sampled_plcp_idx != currsa) {
          std::uint64_t max_lcp_delta2 = text_length -
            (std::max(currsa, prevsa) + plcp_lower_bound);
          if (sampled_plcp_idx + plcp_sampling_rate < text_length) {
            std::uint64_t plcp_upper_bound =
              plcp_buffer[2 * i + 1] + next_sample_dist;
            max_lcp_delta = std::min(max_lcp_delta2,
               plcp_upper_bound -
               plcp_lower_bound);
          } else max_lcp_delta = max_lcp_delta2;
        }

        if (max_lcp_delta > 0) {
//          ++nontrivial_pairs;

          // Compute the segments containing the positions to compare and
          // determine if the pair needs to be split into multiple pairs.
          std::uint64_t seg_1_idx = pos1 / max_halfsegment_size;
          std::uint64_t seg_2_idx = pos2 / max_halfsegment_size;
          std::uint64_t seg_1_beg = seg_1_idx * max_halfsegment_size;
          std::uint64_t seg_2_beg = seg_2_idx * max_halfsegment_size;
          std::uint64_t seg_1_ext_end = std::min(text_length,
              seg_1_beg + max_halfsegment_size + max_overflow_size);
          std::uint64_t seg_2_ext_end = std::min(text_length,
              seg_2_beg + max_halfsegment_size + max_overflow_size);
          std::uint64_t seg_1_maxlcp = seg_1_ext_end - pos1;
          std::uint64_t seg_2_maxlcp = seg_2_ext_end - pos2;
          std::uint64_t seg_maxlcp = std::min(seg_1_maxlcp, seg_2_maxlcp);

          if (max_lcp_delta > seg_maxlcp) {

            // Split the pair into multiple pairs.
            while (max_lcp_delta > seg_maxlcp) {
              std::uint64_t file_id =
                pair_to_int(seg_1_idx, seg_2_idx, n_halfsegments);
              pair_multiwriter->write_to_ith_file(file_id, pos1 - seg_1_beg);
              pair_multiwriter->write_to_ith_file(file_id, pos2 - seg_2_beg);
              pos1 += seg_maxlcp;
              pos2 += seg_maxlcp;
              max_lcp_delta -= seg_maxlcp;
              seg_1_idx = pos1 / max_halfsegment_size;
              seg_2_idx = pos2 / max_halfsegment_size;
              seg_1_beg = seg_1_idx * max_halfsegment_size;
              seg_2_beg = seg_2_idx * max_halfsegment_size;
              seg_1_ext_end = std::min(text_length,
                  seg_1_beg + max_halfsegment_size + max_overflow_size);
              seg_2_ext_end = std::min(text_length,
                  seg_2_beg + max_halfsegment_size + max_overflow_size);
              seg_1_maxlcp = seg_1_ext_end - pos1;
              seg_2_maxlcp = seg_2_ext_end - pos2;
              seg_maxlcp = std::min(seg_1_maxlcp, seg_2_maxlcp);
            }

            // Write the last pair in the chain. Use one of the
            // integers in the pair to endode the upper bound for LCP.
            std::uint64_t file_id =
              pair_to_int(seg_1_idx, seg_2_idx, n_halfsegments);
            std::uint64_t msb_bit =
              (1UL << (8 * sizeof(text_offset_type) - 1));
            if (pos1 == seg_1_beg + max_overflow_size) {

              // We don't have to encode pos1 - seg_1_beg. Instead, we
              // encode max_lcp_delta. Of course we also have to set
              // the MSB bit in both integers to recognize this case.
              pair_multiwriter->write_to_ith_file(
                  file_id, max_lcp_delta + msb_bit);
              pair_multiwriter->write_to_ith_file(
                  file_id, (pos2 - seg_2_beg) + msb_bit);
            } else {

              // Analogously, but we don't encode pos2 - seg_2_beg.
              pair_multiwriter->write_to_ith_file(
                  file_id, (pos1 - seg_1_beg) + msb_bit);
              pair_multiwriter->write_to_ith_file(
                  file_id, max_lcp_delta);
            }
          } else {
            std::uint64_t file_id =
              pair_to_int(seg_1_idx, seg_2_idx, n_halfsegments);
            pair_multiwriter->write_to_ith_file(file_id, pos1 - seg_1_beg);
            pair_multiwriter->write_to_ith_file(file_id, pos2 - seg_2_beg);
          }
        }
      }
    }
#endif

    // Update the number of read SA items.
    prev_sa = sa_buffer[local_buf_filled - 1];
    sa_items_read += local_buf_filled;
  }

  // Stop I/O threads.
  sa_reader->stop_reading();

  // Update I/O volume.
  io_volume +=
    pair_multiwriter->bytes_written() +
    sa_reader->bytes_read();
  total_io_volume += io_volume;

  // Print summary.
  long double total_time = utils::wclock() - start;
  fprintf(stderr, "\r  Compute and store LCP lower bounds: time = %.1Lfs, "
      "I/O = %.2LfMiB/s, total I/O vol = %.2Lfbytes/input symbol\n",
      total_time, ((1.L * io_volume) / (1L << 20)) / total_time,
      (1.L * total_io_volume) / text_length);

//  fprintf(stderr, "nontrivial_pairs = %lu/%lu (%.2Lf%%)\n", nontrivial_pairs,
//      text_length, (100.L * nontrivial_pairs) / text_length);

  // Clean up.
  utils::deallocate(plcp_buffer);
  utils::deallocate(addr_buffer);
  utils::deallocate(sa_buffer);

#ifdef _OPENMP
  utils::deallocate(item_buffer);
#endif

  delete pair_multiwriter;
  delete sa_reader;
}

// Scan SA and for every pair (i, Phi[i]) use PLCP array to find
// the lower bound `lb' on PLCP[i]. Then add pair (i + lb, Phi[i] + lb)
// to the file corresponding to the pairs of segments containing positions
// i + lb and Phi[i] + lb. Large values of `lb' may spawn many pairs.
// A version of the function using text-partitioning.
template<typename text_offset_type>
void compute_lcp_lower_bounds_text_partitioning(
    std::string sa_filename,
    std::string output_filename,
    std::uint64_t text_length,
    std::uint64_t ram_use,
    text_offset_type *sparse_plcp,
    std::uint64_t plcp_sampling_rate,
    std::uint64_t max_halfsegment_size,
    std::uint64_t halfseg_buffers_ram,
    std::uint64_t in_buf_ram,
    std::uint64_t local_buf_ram,
    std::uint64_t max_overflow_size,
    std::uint64_t ***items_per_halfseg_pair,
    std::uint64_t part_id,
    std::uint64_t &total_io_volume) {

  // Compute basic parameters.
  std::uint64_t n_halfsegments =
    (text_length + max_halfsegment_size - 1) / max_halfsegment_size;
  std::uint64_t sparse_plcp_size =
    (text_length + plcp_sampling_rate - 1) / plcp_sampling_rate;
  std::uint64_t sparse_plcp_ram = sparse_plcp_size * sizeof(text_offset_type);

  // Compute how many items belonging for each pairs of halfsegments
  // were processed. Negative values allows skipping items.
  std::int64_t **pairs_processed = new std::int64_t*[n_halfsegments];
  for (std::uint64_t i = 0; i < n_halfsegments; ++i) {
    pairs_processed[i] = new std::int64_t[n_halfsegments];
    std::fill(pairs_processed[i], pairs_processed[i] + n_halfsegments, 0UL);
    for (std::uint64_t j = i; j < n_halfsegments; ++j)
      for (std::uint64_t prev_part_id = 0;
          prev_part_id < part_id; ++prev_part_id)
        pairs_processed[i][j] -=
          (std::int64_t)items_per_halfseg_pair[prev_part_id][i][j];
  }

  // Proportionally enlarge the buffers to fit the available RAM.
  {
    std::uint64_t ram_budget = ram_use -
      sparse_plcp_ram - halfseg_buffers_ram;
    std::uint64_t total_buf_ram = in_buf_ram + local_buf_ram;
    long double scale_factor =
      (long double)ram_budget / (long double)total_buf_ram;
    in_buf_ram = (std::uint64_t)((long double)in_buf_ram * scale_factor);
    local_buf_ram =
      (std::uint64_t)((long double)local_buf_ram * scale_factor);
  }

  // Print the initial message.
  fprintf(stderr, "  Compute and store LCP lower bounds: ");

  // Start the timer.
  long double start = utils::wclock();

  // Create SA reader.
  typedef async_stream_reader<text_offset_type> sa_reader_type;
  sa_reader_type *sa_reader = new sa_reader_type(sa_filename,
      in_buf_ram, std::max(4UL, in_buf_ram / (2UL << 20)));

  // Create a writer for every pairs of halfsegments possibly used.
  std::uint64_t cur_part_min_halfsegment_id_diff = n_halfsegments;
  std::uint64_t cur_part_max_halfsegment_id_diff = 0UL;
  for (std::uint64_t i = 0; i < n_halfsegments; ++i) {
    for (std::uint64_t j = i; j < n_halfsegments; ++j) {
      if (items_per_halfseg_pair[part_id][i][j] > 0) {
        cur_part_max_halfsegment_id_diff =
          std::max(cur_part_max_halfsegment_id_diff, j - i);
        cur_part_min_halfsegment_id_diff =
          std::min(cur_part_min_halfsegment_id_diff, j - i);
      }
    }
  }
  std::uint64_t different_possibly_used_pairs_of_halfsegments = 0;
  for (std::uint64_t i = 0; i < n_halfsegments; ++i) {
    for (std::uint64_t j = i; j < n_halfsegments; ++j) {
      if (cur_part_min_halfsegment_id_diff <= (j - i) + 1 &&
          (j - i) <= cur_part_max_halfsegment_id_diff + 1)
        ++different_possibly_used_pairs_of_halfsegments;
    }
  }
  static const std::uint64_t n_free_buffers = 4;
  std::uint64_t buffer_size = halfseg_buffers_ram /
    (different_possibly_used_pairs_of_halfsegments + n_free_buffers);
  typedef async_multi_stream_writer<text_offset_type> pair_multiwriter_type;
  pair_multiwriter_type *pair_multiwriter = new pair_multiwriter_type(
      different_possibly_used_pairs_of_halfsegments,
      buffer_size, n_free_buffers);

  // Create a map from used halfsegment pairs to a contiguous
  // range of integers. This is needed to use multi stream writer.
  std::uint64_t **halfseg_ids_to_file_id = new std::uint64_t*[n_halfsegments];
  {
    std::uint64_t file_counter = 0;
    for (std::uint64_t i = 0; i < n_halfsegments; ++i) {
      halfseg_ids_to_file_id[i] = new std::uint64_t[n_halfsegments];
      for (std::uint64_t j = i; j < n_halfsegments; ++j) {
        if (cur_part_min_halfsegment_id_diff <= (j - i) + 1 &&
            (j - i) <= cur_part_max_halfsegment_id_diff + 1) {
          pair_multiwriter->add_file(output_filename + ".pairs." +
              utils::intToStr(i) + "_" + utils::intToStr(j));
          halfseg_ids_to_file_id[i][j] = file_counter++;
        }
      }
    }
  }

  // Allocate buffers.
#ifdef _OPENMP
  typedef ext_buf_item<std::uint64_t, std::uint64_t, std::uint64_t>
    buf_item_type;
  std::uint64_t local_buf_size = local_buf_ram /
    (1 * sizeof(buf_item_type) +
     1 * sizeof(text_offset_type) +
     3 * sizeof(std::uint64_t));
#else
  std::uint64_t local_buf_size = local_buf_ram /
    (1 * sizeof(text_offset_type) +
     3 * sizeof(std::uint64_t));
#endif

#ifdef _OPENMP
  buf_item_type *item_buffer =
    utils::allocate_array<buf_item_type>(local_buf_size);
#endif

  text_offset_type *sa_buffer =
    utils::allocate_array<text_offset_type>(local_buf_size);
  std::uint64_t *addr_buffer =
    utils::allocate_array<std::uint64_t>(local_buf_size);
  std::uint64_t *plcp_buffer =
    utils::allocate_array<std::uint64_t>(2 * local_buf_size);

  // Start processing.
  std::uint64_t sa_items_read = 0;
  std::uint64_t step_counter = 0;
//  std::uint64_t nontrivial_pairs = 0;
  text_offset_type prev_sa = text_length;

  while (sa_items_read < text_length) {

    // Print progress message.
    ++step_counter;
    if (step_counter == (1UL << 25) / local_buf_size) {
      step_counter = 0;
      long double elapsed = utils::wclock() - start;
      std::uint64_t io_volume =
        sa_reader->bytes_read() +
        pair_multiwriter->bytes_written();
      fprintf(stderr, "\r  Compute and store LCP lower bounds: "
          "%.1Lf%%, time = %.1Lfs, I/O = %.2LfMiB/s",
          (100.L * sa_items_read) / text_length, elapsed,
          (1.L * io_volume / (1L << 20)) / elapsed);
    }

    // Fill in the buffer with SA values.
    std::uint64_t local_buf_filled =
      std::min(text_length - sa_items_read, local_buf_size);
    sa_reader->read(sa_buffer, local_buf_filled);

    // Process items in the buffer.
    // First, fetch PLCP values.
#ifdef _OPENMP
#ifdef PRODUCTION_VER
    #pragma omp parallel for
    for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
      std::uint64_t addr = sa_buffer[j] / plcp_sampling_rate;
      addr_buffer[j] = addr;
    }

    #pragma omp parallel for
    for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
      std::uint64_t addr = addr_buffer[j];
      plcp_buffer[2 * j] = sparse_plcp[addr];
      if (addr + 1 < sparse_plcp_size)
        plcp_buffer[2 * j + 1] = sparse_plcp[addr + 1];
    }
#else
    #pragma omp parallel for
    for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
      std::uint64_t addr = sa_buffer[j] / plcp_sampling_rate;
      addr_buffer[j] = addr;
      plcp_buffer[2 * j] = sparse_plcp[addr];
      if (addr + 1 < sparse_plcp_size)
        plcp_buffer[2 * j + 1] = sparse_plcp[addr + 1];
    }
#endif
#else
    for (std::uint64_t j = 0; j < local_buf_filled; ++j)
      addr_buffer[j] = sa_buffer[j] / plcp_sampling_rate;
    for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
      std::uint64_t addr = addr_buffer[j];
      plcp_buffer[2 * j] = sparse_plcp[addr];
      if (addr + 1 < sparse_plcp_size)
        plcp_buffer[2 * j + 1] = sparse_plcp[addr + 1];
    }
#endif

#ifdef _OPENMP
    {
      std::uint64_t different_halfsegments_pairs =
        (n_halfsegments * (n_halfsegments + 1)) / 2;
      std::uint64_t max_threads = omp_get_max_threads();
      std::uint64_t max_range_size =
        (local_buf_filled + max_threads - 1) / max_threads;
      std::uint64_t n_ranges =
        (local_buf_filled + max_range_size - 1) / max_range_size;

      #pragma omp parallel num_threads(n_ranges)
      {
        std::uint64_t thread_id = omp_get_thread_num();
        std::uint64_t range_beg = thread_id * max_range_size;
        std::uint64_t range_end =
          std::min(range_beg + max_range_size, local_buf_filled);
        std::uint64_t prev_src_halfseg_id =
          ((range_beg == 0) ?
           ((std::uint64_t)prev_sa) :
           ((std::uint64_t)sa_buffer[range_beg - 1])) / max_halfsegment_size;

        for (std::uint64_t i = range_beg; i < range_end; ++i) {
          std::uint64_t prevsa = (i == 0) ?
            (std::uint64_t)prev_sa :
            (std::uint64_t)sa_buffer[i - 1];
          std::uint64_t currsa = sa_buffer[i];
          std::uint64_t cur_src_halfseg_id = currsa / max_halfsegment_size;
          std::uint64_t h1 = cur_src_halfseg_id;
          std::uint64_t h2 = prev_src_halfseg_id;
          prev_src_halfseg_id = cur_src_halfseg_id;

          if (prevsa == text_length) {

            // Don't process it and don't decrement any counters.
            item_buffer[i] = buf_item_type(different_halfsegments_pairs,
                n_halfsegments, 0, 0, 0, false);
            continue;
          }

          std::uint64_t pos1 = currsa;
          std::uint64_t pos2 = prevsa;
          if (pos1 > pos2) {
            std::swap(pos1, pos2);
            std::swap(h1, h2);
          }

          // Optimization.
          if (pairs_processed[h1][h2] + (std::int64_t)local_buf_filled <= 0 ||
              pairs_processed[h1][h2] >=
              (std::int64_t)items_per_halfseg_pair[part_id][h1][h2]) {

            // Don't process it but decrement the counters.
            item_buffer[i] = buf_item_type(
                different_halfsegments_pairs, h1, h2, 0, 0, false);
            continue;
          }

          // Compute lower bound for PLCP[i].
          std::uint64_t sampled_plcp_ptr = addr_buffer[i];
          std::uint64_t sampled_plcp_idx =
            sampled_plcp_ptr * plcp_sampling_rate;
          std::uint64_t prev_sample_dist = currsa - sampled_plcp_idx;
          std::uint64_t next_sample_dist =
            plcp_sampling_rate - prev_sample_dist;
          std::uint64_t plcp_val = plcp_buffer[2 * i];
          std::uint64_t plcp_lower_bound = std::max(0L,
              (std::int64_t)plcp_val -
              (std::int64_t)prev_sample_dist);
          pos1 += plcp_lower_bound;
          pos2 += plcp_lower_bound;

          if (pos2 >= text_length) {

            // Don't process it, but decrement the counter!
            item_buffer[i] = buf_item_type(
                different_halfsegments_pairs, h1, h2, 0, 0, false);
            continue;
          }

          // Compute the maximal LCP delta in
          // case PLCP[currsa] was not sampled.
          std::uint64_t max_lcp_delta = 0;
          if (sampled_plcp_idx != currsa) {
            std::uint64_t max_lcp_delta2 = text_length -
              (std::max(currsa, prevsa) + plcp_lower_bound);
            if (sampled_plcp_idx + plcp_sampling_rate < text_length) {
              std::uint64_t plcp_upper_bound =
                plcp_buffer[2 * i + 1] + next_sample_dist;
              max_lcp_delta = std::min(max_lcp_delta2,
                 plcp_upper_bound -
                 plcp_lower_bound);
            } else max_lcp_delta = max_lcp_delta2;
          }

          if (max_lcp_delta > 0) {

            // Compute the segments containing the positions to compare and
            // determine if the pair needs to be split into multiple pairs.
            std::uint64_t seg_1_idx = pos1 / max_halfsegment_size;
            std::uint64_t seg_2_idx = pos2 / max_halfsegment_size;
            std::uint64_t seg_1_beg = seg_1_idx * max_halfsegment_size;
            std::uint64_t seg_2_beg = seg_2_idx * max_halfsegment_size;
            std::uint64_t seg_1_ext_end = std::min(text_length,
                seg_1_beg + max_halfsegment_size + max_overflow_size);
            std::uint64_t seg_2_ext_end = std::min(text_length,
                seg_2_beg + max_halfsegment_size + max_overflow_size);
            std::uint64_t seg_1_maxlcp = seg_1_ext_end - pos1;
            std::uint64_t seg_2_maxlcp = seg_2_ext_end - pos2;
            std::uint64_t seg_maxlcp = std::min(seg_1_maxlcp, seg_2_maxlcp);

            // Add the appropriate item to the buffer.
            if (max_lcp_delta > seg_maxlcp) {

              // Special pair. Will be spli into chain of pairs when
              // writing. Process and decrement the counter.
              item_buffer[i] = buf_item_type(0, h1, h2, 0, 0, true);
            } else {

              // Normal pair. Will be written as a single
              // pair. Process and cecrement the dounter.
              std::uint64_t file_id =
                halfseg_ids_to_file_id[seg_1_idx][seg_2_idx];
              item_buffer[i] = buf_item_type(file_id, h1, h2,
                  pos1 - seg_1_beg, pos2 - seg_2_beg, false);
            }
          } else {

            // Don't process it, but decrement the counter!
            item_buffer[i] = buf_item_type(
                different_halfsegments_pairs, h1, h2, 0, 0, false);
          }
        }
      }

      // Distribute the items from the pair buffer to files.
      for (std::uint64_t i = 0; i < local_buf_filled; ++i) {
        std::uint64_t addr = item_buffer[i].addr;
        if (addr != different_halfsegments_pairs) {

          // Halfsegment ids containing the positions.
          std::uint64_t h1 = item_buffer[i].id1;
          std::uint64_t h2 = item_buffer[i].id2;

          // Check if current item should be processed (in this text part).
          ++pairs_processed[h1][h2];
          if (pairs_processed[h1][h2] > 0L &&
              pairs_processed[h1][h2] <=
              (std::int64_t)items_per_halfseg_pair[part_id][h1][h2]) {
            if (item_buffer[i].special == false) {

              // Write a single lcp pair.
              pair_multiwriter->write_to_ith_file(addr, item_buffer[i].val1);
              pair_multiwriter->write_to_ith_file(addr, item_buffer[i].val2);
            } else {

              // Recompute and write the chain of pairs.
              std::uint64_t prevsa = (i == 0) ?
                (std::uint64_t)prev_sa :
                (std::uint64_t)sa_buffer[i - 1];
              std::uint64_t currsa = sa_buffer[i];

              // Compute lower bound for PLCP[i].
              std::uint64_t sampled_plcp_ptr = addr_buffer[i];
              std::uint64_t sampled_plcp_idx =
                sampled_plcp_ptr * plcp_sampling_rate;
              std::uint64_t prev_sample_dist = currsa - sampled_plcp_idx;
              std::uint64_t next_sample_dist =
                plcp_sampling_rate - prev_sample_dist;
              std::uint64_t plcp_val = plcp_buffer[2 * i];
              std::uint64_t plcp_lower_bound = std::max(0L,
                  (std::int64_t)plcp_val -
                  (std::int64_t)prev_sample_dist);
              std::uint64_t pos1 = currsa + plcp_lower_bound;
              std::uint64_t pos2 = prevsa + plcp_lower_bound;
              if (pos1 > pos2)
                std::swap(pos1, pos2);

              // Compute the maximal LCP delta in
              // case PLCP[currsa] was not sampled.
              std::uint64_t max_lcp_delta = 0;
              if (sampled_plcp_idx != currsa) {
                std::uint64_t max_lcp_delta2 = text_length -
                  (std::max(currsa, prevsa) + plcp_lower_bound);
                if (sampled_plcp_idx + plcp_sampling_rate < text_length) {
                  std::uint64_t plcp_upper_bound =
                    plcp_buffer[2 * i + 1] + next_sample_dist;
                  max_lcp_delta = std::min(max_lcp_delta2,
                     plcp_upper_bound -
                     plcp_lower_bound);
                } else max_lcp_delta = max_lcp_delta2;
              }

              // Compute the segments containing the positions to compare and
              // determine if the pair needs to be split into multiple pairs.
              std::uint64_t seg_1_idx = pos1 / max_halfsegment_size;
              std::uint64_t seg_2_idx = pos2 / max_halfsegment_size;
              std::uint64_t seg_1_beg = seg_1_idx * max_halfsegment_size;
              std::uint64_t seg_2_beg = seg_2_idx * max_halfsegment_size;
              std::uint64_t seg_1_ext_end = std::min(text_length,
                  seg_1_beg + max_halfsegment_size + max_overflow_size);
              std::uint64_t seg_2_ext_end = std::min(text_length,
                  seg_2_beg + max_halfsegment_size + max_overflow_size);
              std::uint64_t seg_1_maxlcp = seg_1_ext_end - pos1;
              std::uint64_t seg_2_maxlcp = seg_2_ext_end - pos2;
              std::uint64_t seg_maxlcp = std::min(seg_1_maxlcp, seg_2_maxlcp);

              // Split the pair into multiple pairs.
              while (max_lcp_delta > seg_maxlcp) {
                std::uint64_t file_id =
                  halfseg_ids_to_file_id[seg_1_idx][seg_2_idx];
                pair_multiwriter->write_to_ith_file(
                    file_id, pos1 - seg_1_beg);
                pair_multiwriter->write_to_ith_file(
                    file_id, pos2 - seg_2_beg);
                pos1 += seg_maxlcp;
                pos2 += seg_maxlcp;
                max_lcp_delta -= seg_maxlcp;
                seg_1_idx = pos1 / max_halfsegment_size;
                seg_2_idx = pos2 / max_halfsegment_size;
                seg_1_beg = seg_1_idx * max_halfsegment_size;
                seg_2_beg = seg_2_idx * max_halfsegment_size;
                seg_1_ext_end = std::min(text_length,
                    seg_1_beg + max_halfsegment_size + max_overflow_size);
                seg_2_ext_end = std::min(text_length,
                    seg_2_beg + max_halfsegment_size + max_overflow_size);
                seg_1_maxlcp = seg_1_ext_end - pos1;
                seg_2_maxlcp = seg_2_ext_end - pos2;
                seg_maxlcp = std::min(seg_1_maxlcp, seg_2_maxlcp);
              }

              // Write the last pair in the chain. Use one of the
              // integers in the pair to encode the upper bound for LCP.
              std::uint64_t file_id =
                halfseg_ids_to_file_id[seg_1_idx][seg_2_idx];
              std::uint64_t msb_bit =
                (1UL << (8 * sizeof(text_offset_type) - 1));
              if (pos1 == seg_1_beg + max_overflow_size) {

                // We don't have to encode pos1 - seg_1_beg. Instead, we
                // encode max_lcp_delta. Of course we also have to set
                // the MSB bit in both integers to recognize this case.
                pair_multiwriter->write_to_ith_file(
                    file_id, max_lcp_delta + msb_bit);
                pair_multiwriter->write_to_ith_file(
                    file_id, (pos2 - seg_2_beg) + msb_bit);
              } else {

                // Analogously, but we don't encode pos2 - seg_2_beg.
                pair_multiwriter->write_to_ith_file(
                    file_id, (pos1 - seg_1_beg) + msb_bit);
                pair_multiwriter->write_to_ith_file(
                    file_id, max_lcp_delta);
              }
            }
          }
        } else {

          // The pair was not processed. Check if we need to update the count.
          if ((std::uint64_t)item_buffer[i].id1 != n_halfsegments) {

            // Halfsegment ids containing the positions.
            std::uint64_t h1 = item_buffer[i].id1;
            std::uint64_t h2 = item_buffer[i].id2;

            // Check if current item should be processed (in this text part).
            ++pairs_processed[h1][h2];
          }
        }
      }
    }

#else
    {
      std::uint64_t prev_src_halfseg_id =
        (std::uint64_t)prev_sa / max_halfsegment_size;
      for (std::uint64_t i = 0; i < local_buf_filled; ++i) {
        std::uint64_t prevsa = (i == 0) ?
          (std::uint64_t)prev_sa :
          (std::uint64_t)sa_buffer[i - 1];
        std::uint64_t currsa = sa_buffer[i];
        std::uint64_t cur_src_halfseg_id = currsa / max_halfsegment_size;
        std::uint64_t h1 = cur_src_halfseg_id;
        std::uint64_t h2 = prev_src_halfseg_id;
        prev_src_halfseg_id = cur_src_halfseg_id;

        if (prevsa == text_length)
          continue;

        std::uint64_t pos1 = currsa;
        std::uint64_t pos2 = prevsa;
        if (pos1 > pos2) {
          std::swap(pos1, pos2);
          std::swap(h1, h2);
        }

        // Check if current item should be processed (in this text part).
        ++pairs_processed[h1][h2];
        if (pairs_processed[h1][h2] <= 0L ||
            pairs_processed[h1][h2] >
            (std::int64_t)items_per_halfseg_pair[part_id][h1][h2])
          continue;

        // Compute lower bound for PLCP[i].
        std::uint64_t sampled_plcp_ptr = addr_buffer[i];
        std::uint64_t sampled_plcp_idx =
          sampled_plcp_ptr * plcp_sampling_rate;
        std::uint64_t prev_sample_dist = currsa - sampled_plcp_idx;
        std::uint64_t next_sample_dist =
          plcp_sampling_rate - prev_sample_dist;
        std::uint64_t plcp_val = plcp_buffer[2 * i];
        std::uint64_t plcp_lower_bound = std::max(0L,
            (std::int64_t)plcp_val -
            (std::int64_t)prev_sample_dist);
        pos1 += plcp_lower_bound;
        pos2 += plcp_lower_bound;
        if (pos2 >= text_length)
          continue;

        // Compute the maximal LCP delta in
        // case PLCP[currsa] was not sampled.
        std::uint64_t max_lcp_delta = 0;
        if (sampled_plcp_idx != currsa) {
          std::uint64_t max_lcp_delta2 = text_length -
            (std::max(currsa, prevsa) + plcp_lower_bound);
          if (sampled_plcp_idx + plcp_sampling_rate < text_length) {
            std::uint64_t plcp_upper_bound =
              plcp_buffer[2 * i + 1] + next_sample_dist;
            max_lcp_delta = std::min(max_lcp_delta2,
                plcp_upper_bound -
                plcp_lower_bound);
          } else max_lcp_delta = max_lcp_delta2;
        }

        if (max_lcp_delta > 0) {
//          ++nontrivial_pairs;

          // Compute the segments containing the positions to compare and
          // determine if the pair needs to be split into multiple pairs.
          std::uint64_t seg_1_idx = pos1 / max_halfsegment_size;
          std::uint64_t seg_2_idx = pos2 / max_halfsegment_size;
          std::uint64_t seg_1_beg = seg_1_idx * max_halfsegment_size;
          std::uint64_t seg_2_beg = seg_2_idx * max_halfsegment_size;
          std::uint64_t seg_1_ext_end = std::min(text_length,
              seg_1_beg + max_halfsegment_size + max_overflow_size);
          std::uint64_t seg_2_ext_end = std::min(text_length,
              seg_2_beg + max_halfsegment_size + max_overflow_size);
          std::uint64_t seg_1_maxlcp = seg_1_ext_end - pos1;
          std::uint64_t seg_2_maxlcp = seg_2_ext_end - pos2;
          std::uint64_t seg_maxlcp = std::min(seg_1_maxlcp, seg_2_maxlcp);

          if (max_lcp_delta > seg_maxlcp) {

            // Split the pair into multiple pairs.
            while (max_lcp_delta > seg_maxlcp) {
              std::uint64_t file_id =
                halfseg_ids_to_file_id[seg_1_idx][seg_2_idx];
              pair_multiwriter->write_to_ith_file(file_id, pos1 - seg_1_beg);
              pair_multiwriter->write_to_ith_file(file_id, pos2 - seg_2_beg);
              pos1 += seg_maxlcp;
              pos2 += seg_maxlcp;
              max_lcp_delta -= seg_maxlcp;
              seg_1_idx = pos1 / max_halfsegment_size;
              seg_2_idx = pos2 / max_halfsegment_size;
              seg_1_beg = seg_1_idx * max_halfsegment_size;
              seg_2_beg = seg_2_idx * max_halfsegment_size;
              seg_1_ext_end = std::min(text_length,
                  seg_1_beg + max_halfsegment_size + max_overflow_size);
              seg_2_ext_end = std::min(text_length,
                  seg_2_beg + max_halfsegment_size + max_overflow_size);
              seg_1_maxlcp = seg_1_ext_end - pos1;
              seg_2_maxlcp = seg_2_ext_end - pos2;
              seg_maxlcp = std::min(seg_1_maxlcp, seg_2_maxlcp);
            }

            // Write the last pair in the chain. Use one of the
            // integers in the pair to endode the upper bound for LCP.
            std::uint64_t file_id =
              halfseg_ids_to_file_id[seg_1_idx][seg_2_idx];
            std::uint64_t msb_bit =
              (1UL << (8 * sizeof(text_offset_type) - 1));
            if (pos1 == seg_1_beg + max_overflow_size) {

              // We don't have to encode pos1 - seg_1_beg. Instead, we
              // encode max_lcp_delta. Of course we also have to set
              // the MSB bit in both integers to recognize this case.
              pair_multiwriter->write_to_ith_file(
                  file_id, max_lcp_delta + msb_bit);
              pair_multiwriter->write_to_ith_file(
                  file_id, (pos2 - seg_2_beg) + msb_bit);
            } else {

              // Analogously, but we don't encode pos2 - seg_2_beg.
              pair_multiwriter->write_to_ith_file(
                  file_id, (pos1 - seg_1_beg) + msb_bit);
              pair_multiwriter->write_to_ith_file(
                  file_id, max_lcp_delta);
            }
          } else {
            std::uint64_t file_id =
              halfseg_ids_to_file_id[seg_1_idx][seg_2_idx];
            pair_multiwriter->write_to_ith_file(file_id, pos1 - seg_1_beg);
            pair_multiwriter->write_to_ith_file(file_id, pos2 - seg_2_beg);
          }
        }
      }
    }
#endif

    // Update the number of read SA items.
    prev_sa = (std::uint64_t)sa_buffer[local_buf_filled - 1];
    sa_items_read += local_buf_filled;
  }

  // Stop I/O threads.
  sa_reader->stop_reading();

  // Update I/O volume.
  std::uint64_t io_vol =
    pair_multiwriter->bytes_written() +
    sa_reader->bytes_read();
  total_io_volume += io_vol;

  // Print summary.
  long double total_time = utils::wclock() - start;
  fprintf(stderr, "\r  Compute and store LCP lower bounds: time = %.1Lfs, "
      "I/O = %.2LfMiB/s, total I/O vol = %.2Lfbytes/input symbol\n",
      total_time, ((1.L * io_vol) / (1L << 20)) / total_time,
      (1.L * total_io_volume) / text_length);
//  fprintf(stderr, "nontrivial_pairs = %lu/%lu (%.2Lf%%)\n",
//      nontrivial_pairs, text_length,
//      (100.L * nontrivial_pairs) / text_length);

  // Clean up.
  utils::deallocate(plcp_buffer);
  utils::deallocate(addr_buffer);
  utils::deallocate(sa_buffer);

#ifdef _OPENMP
  utils::deallocate(item_buffer);
#endif

  for (std::uint64_t i = 0; i < n_halfsegments; ++i)
    delete[] halfseg_ids_to_file_id[i];
  delete[] halfseg_ids_to_file_id;
  delete pair_multiwriter;
  delete sa_reader;
  for (std::uint64_t i = 0; i < n_halfsegments; ++i)
    delete[] pairs_processed[i];
  delete[] pairs_processed;
}

}  // namespace em_sparse_phi_private

#endif  // __SRC_EM_SPARSE_PHI_SRC_COMPUTE_LCP_LOWER_BOUNDS_HPP_INCLUDED
