/**
 * @file    src/em_sparse_phi_src/compute_final_lcp.hpp
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

#ifndef __SRC_EM_SPARSE_PHI_SRC_COMPUTE_FINAL_LCP_HPP_INCLUDED
#define __SRC_EM_SPARSE_PHI_SRC_COMPUTE_FINAL_LCP_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <algorithm>
#include <omp.h>

#include "io/async_stream_reader.hpp"
#include "io/async_stream_vbyte_reader_multipart.hpp"
#include "io/async_stream_writer.hpp"
#include "io/async_multi_stream_vbyte_reader.hpp"
#include "utils.hpp"


namespace em_sparse_phi_private {

template<typename T, typename S>
struct buf_item_2 {
  buf_item_2() {}
  buf_item_2(
      T idx1,
      T idx2,
      S lb,
      std::uint8_t tp)
    : seg_1_idx(idx1),
    seg_2_idx(idx2),
    plcp_lower_bound(lb),
    type(tp) {}

  T seg_1_idx;
  T seg_2_idx;
  S plcp_lower_bound;
  std::uint8_t type;
} __attribute__ ((packed));

struct ext_buf_item_2 {
  ext_buf_item_2() {}
  ext_buf_item_2(std::uint64_t lb, std::uint8_t tp)
    : m_plcp_lower_bound(lb), m_type(tp) {}

  std::uint64_t m_plcp_lower_bound;
  std::uint8_t m_type;
} __attribute__ ((packed));


// Compute the final LCP array. Return
// max and sum of LCP values.
template<typename text_offset_type>
void compute_final_lcp(
    std::string sa_filename,
    std::string output_filename,
    std::uint64_t text_length,
    text_offset_type *sparse_plcp,
    std::uint64_t plcp_sampling_rate,
    std::uint64_t max_halfsegment_size,
    std::uint64_t halfseg_buffers_ram,
    std::uint64_t in_buf_ram,
    std::uint64_t out_buf_ram,
    std::uint64_t local_buf_ram,
    std::uint64_t max_overflow_size,
    std::string **delta_filenames,
    std::uint64_t &lcp_sum,
    std::uint64_t &max_lcp,
    std::uint64_t &total_io_volume) {

  std::uint64_t n_halfsegments =
    (text_length + max_halfsegment_size - 1) / max_halfsegment_size;
  std::uint64_t sparse_plcp_size =
    (text_length + plcp_sampling_rate - 1) / plcp_sampling_rate;

  std::uint64_t local_lcp_sum = 0;
  std::uint64_t local_max_lcp = 0;

  fprintf(stderr, "\nCompute final LCP array: ");
  long double start = utils::wclock();

#ifdef _OPENMP
  typedef buf_item_2<std::uint64_t, std::uint64_t> buf_item_type;
  std::uint64_t local_buf_size = local_buf_ram /
    (1 * sizeof(buf_item_type) +
     1 * sizeof(text_offset_type) +
     3 * sizeof(std::uint64_t));
#else
  std::uint64_t local_buf_size = local_buf_ram /
    (1 * sizeof(text_offset_type) +
     3 * sizeof(std::uint64_t));
#endif

  // Allocate buffers.
#ifdef _OPENMP
  buf_item_type *item_buffer =
    utils::allocate_array<buf_item_type>(local_buf_size);
#endif
  text_offset_type *sa_buffer =
    utils::allocate_array<text_offset_type>(local_buf_size);
  std::uint64_t *addr_buffer =
    utils::allocate_array<std::uint64_t>(local_buf_size);
  std::uint64_t *plcp_buffer =
    utils::allocate_array<std::uint64_t>(local_buf_size * 2);

  std::uint64_t **file_id = new std::uint64_t*[n_halfsegments];

  // Create reader for every pair of halfsegments.
  typedef async_multi_stream_vbyte_reader pair_multireader_type;
  std::uint64_t nonempty_halfsegment_pairs = 0;
  for (std::uint64_t i = 0; i < n_halfsegments; ++i) {
    for (std::uint64_t j = i; j < n_halfsegments; ++j) {
      std::string filename = delta_filenames[i][j];
      if (utils::file_exists(filename) && utils::file_size(filename) > 0)
        ++nonempty_halfsegment_pairs;
    }
  }

  // Initialize multireader of pairs.
  std::uint64_t buffer_size = halfseg_buffers_ram /
    std::max(1UL, nonempty_halfsegment_pairs);
  pair_multireader_type *pair_multireader = new pair_multireader_type(
      nonempty_halfsegment_pairs, buffer_size);

  // Add files to multireader.
  std::uint64_t file_id_counter = 0;
  for (std::uint64_t i = 0; i < n_halfsegments; ++i) {
    file_id[i] = new std::uint64_t[n_halfsegments];
    for (std::uint64_t j = i; j < n_halfsegments; ++j) {
      std::string filename = delta_filenames[i][j];
      if (utils::file_exists(filename) && utils::file_size(filename) > 0) {
        file_id[i][j] = file_id_counter++;
        pair_multireader->add_file(filename);
      }
    }
  }

  // Initialize SA reader.
  typedef async_stream_reader<text_offset_type> sa_reader_type;
  sa_reader_type *sa_reader = new sa_reader_type(
      sa_filename, in_buf_ram, std::max(4UL, in_buf_ram / (2UL << 20)));

  // Initialize writer of the final LCP array.
  typedef async_stream_writer<text_offset_type> lcp_writer_type;
  lcp_writer_type *lcp_writer = new lcp_writer_type(
      output_filename, out_buf_ram, std::max(4UL, out_buf_ram / (2UL << 20)));

  // Scan suffix array and with the help
  // of PLCP array compute all LCP values.
  std::uint64_t prev_sa = text_length;
  std::uint64_t sa_items_read = 0;
  std::uint64_t step_counter = 0;
  while (sa_items_read < text_length) {
    ++step_counter;
    if (step_counter == (1UL << 25) / local_buf_size) {
      step_counter = 0;
      long double elapsed = utils::wclock() - start;
      std::uint64_t io_volume =
        sa_reader->bytes_read() +
        lcp_writer->bytes_written() +
        pair_multireader->bytes_read();
      fprintf(stderr, "\rCompute final LCP array: %.1Lf%%, time = %.1Lfs, "
          "I/O = %.1LfMiB/s", (100.L * sa_items_read) / text_length,
          elapsed, (1.L * io_volume / (1L << 20)) / elapsed);
    }

    // Fill in the buffer.
    std::uint64_t local_buf_filled =
      std::min(text_length - sa_items_read, local_buf_size);
    sa_reader->read(sa_buffer, local_buf_filled);

#ifdef _OPENMP
#ifdef PRODUCTION_VER
    #pragma omp parallel for
    for (std::uint64_t j = 0; j < local_buf_filled; ++j)
      addr_buffer[j] = (std::uint64_t)sa_buffer[j] / plcp_sampling_rate;

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
      std::uint64_t addr = (std::uint64_t)sa_buffer[j] / plcp_sampling_rate;
      addr_buffer[j] = addr;
      plcp_buffer[2 * j] = sparse_plcp[addr];
      if (addr + 1 < sparse_plcp_size)
        plcp_buffer[2 * j + 1] = sparse_plcp[addr + 1];
    }
#endif
#else
    for (std::uint64_t j = 0; j < local_buf_filled; ++j)
      addr_buffer[j] = (std::uint64_t)sa_buffer[j] / plcp_sampling_rate;
    for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
      std::uint64_t addr = addr_buffer[j];
      plcp_buffer[2 * j] = sparse_plcp[addr];
      if (addr + 1 < sparse_plcp_size)
        plcp_buffer[2 * j + 1] = sparse_plcp[addr + 1];
    }
#endif

#ifdef _OPENMP

    // Process buffer.
    #pragma omp parallel for
    for (std::uint64_t i = 0; i < local_buf_filled; ++i) {
      std::uint64_t currsa = sa_buffer[i];
      std::uint64_t prevsa = (i == 0) ?
        prev_sa : (std::uint64_t)sa_buffer[i - 1];

      if (prevsa == text_length) {
        item_buffer[i] = buf_item_type(0, 0, 0, 2);
        continue;
      }

      // Find lower bound on PLCP[i] using sparse_plcp.
      std::uint64_t sampled_plcp_ptr = addr_buffer[i];
      std::uint64_t sampled_plcp_idx = sampled_plcp_ptr * plcp_sampling_rate;
      std::uint64_t plcp_val = plcp_buffer[2 * i];
      std::uint64_t sample_dist = currsa - sampled_plcp_idx;
      std::uint64_t plcp_lower_bound =
        std::max(0L, (std::int64_t)plcp_val - (std::int64_t)sample_dist);
      std::uint64_t pos1 = currsa + plcp_lower_bound;
      std::uint64_t pos2 = prevsa + plcp_lower_bound;
      if (pos1 > pos2) std::swap(pos1, pos2);
      if (pos2 >= text_length) {
        item_buffer[i] = buf_item_type(0, 0, 0, 2);
        continue;
      }

      // Compute the maximal LCP delta in
      // case PLCP[currsa] was not sampled.
      std::uint64_t max_lcp_delta = 0;
      if (sampled_plcp_idx != currsa) {
        std::uint64_t max_lcp_delta2 = text_length -
          (std::max(currsa, prevsa) + plcp_lower_bound);

        if (sampled_plcp_idx + plcp_sampling_rate < text_length) {
          std::uint64_t dist = plcp_sampling_rate - sample_dist;
          max_lcp_delta = std::min(max_lcp_delta2,
              (plcp_buffer[2 * i + 1] + dist) - plcp_lower_bound);
        } else max_lcp_delta = max_lcp_delta2;
      }

      if (max_lcp_delta == 0)
        item_buffer[i] = buf_item_type( 0, 0, plcp_lower_bound, 0);
      else {

        // Compute lcp_delta.
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

        if (max_lcp_delta > seg_maxlcp)
          item_buffer[i] = buf_item_type(0, 0, 0, 2);
        else item_buffer[i] = buf_item_type(
            seg_1_idx, seg_2_idx, plcp_lower_bound, 1);
      }
    }

    for (std::uint64_t i = 0; i < local_buf_filled; ++i) {
      if (item_buffer[i].type == 0) {
        std::uint64_t plcp_lower_bound = item_buffer[i].plcp_lower_bound;
        lcp_writer->write(plcp_lower_bound);
        local_lcp_sum += plcp_lower_bound;
        local_max_lcp = std::max(local_max_lcp, plcp_lower_bound);
      }  else if (item_buffer[i].type == 1) {

        // Common case.
        std::uint64_t seg_1_idx = item_buffer[i].seg_1_idx;
        std::uint64_t seg_2_idx = item_buffer[i].seg_2_idx;
        std::uint64_t plcp_lower_bound = item_buffer[i].plcp_lower_bound;
        std::uint64_t lcp_delta =
          pair_multireader->read_from_ith_file(file_id[seg_1_idx][seg_2_idx]);
        lcp_writer->write(plcp_lower_bound + lcp_delta);
        local_lcp_sum += plcp_lower_bound + lcp_delta;
        local_max_lcp = std::max(local_max_lcp, plcp_lower_bound + lcp_delta);
      } else {

        // Recompute everything to decode the chain of pairs.
        std::uint64_t currsa = sa_buffer[i];
        std::uint64_t prevsa = (i == 0) ?
          prev_sa : (std::uint64_t)sa_buffer[i - 1];

        if (prevsa == text_length) {
          lcp_writer->write(0);
          continue;
        }

        // Find lower bound on PLCP[i] using sparse_plcp.
        std::uint64_t sampled_plcp_ptr = addr_buffer[i];
        std::uint64_t sampled_plcp_idx =
          sampled_plcp_ptr * plcp_sampling_rate;
        std::uint64_t plcp_val = plcp_buffer[2 * i];
        std::uint64_t sample_dist = currsa - sampled_plcp_idx;
        std::uint64_t plcp_lower_bound =
          std::max(0L, (std::int64_t)plcp_val - (std::int64_t)sample_dist);
        std::uint64_t pos1 = currsa + plcp_lower_bound;
        std::uint64_t pos2 = prevsa + plcp_lower_bound;
        if (pos1 > pos2) std::swap(pos1, pos2);
        if (pos2 >= text_length) {
          lcp_writer->write(plcp_lower_bound);
          local_lcp_sum += plcp_lower_bound;
          local_max_lcp = std::max(local_max_lcp, plcp_lower_bound);
          continue;
        }

        // Compute the maximal LCP delta in
        // case PLCP[currsa] was not sampled.
        std::uint64_t max_lcp_delta = 0;
        if (sampled_plcp_idx != currsa) {
          std::uint64_t max_lcp_delta2 = text_length -
            (std::max(currsa, prevsa) + plcp_lower_bound);

          if (sampled_plcp_idx + plcp_sampling_rate < text_length) {
            std::uint64_t dist = plcp_sampling_rate - sample_dist;
            max_lcp_delta = std::min(max_lcp_delta2,
                (plcp_buffer[2 * i + 1] + dist) - plcp_lower_bound);
          } else max_lcp_delta = max_lcp_delta2;
        }

        // Compute lcp_delta.
        std::uint64_t lcp_delta = 0;
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

        // Invariant: max_lcp_delta > seg_maxlcp
        // We have to read all pairs that were produced for this
        // chain and simultaneously infer the actual value of lcp_delta
        // even if the LCP is determined in one of the earlier segments.
        bool lcp_determined = false;
        while (max_lcp_delta > seg_maxlcp) {
          std::uint64_t local_lcp_delta =
            pair_multireader->read_from_ith_file(
                file_id[seg_1_idx][seg_2_idx]);

          if (lcp_determined == false) {
            lcp_delta += local_lcp_delta;
            if (local_lcp_delta < seg_maxlcp)
              lcp_determined = true;
          }

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

        std::uint64_t local_lcp_delta =
          pair_multireader->read_from_ith_file(file_id[seg_1_idx][seg_2_idx]);

        if (lcp_determined == false)
          lcp_delta += local_lcp_delta;

        lcp_writer->write(plcp_lower_bound + lcp_delta);
        local_lcp_sum += plcp_lower_bound + lcp_delta;
        local_max_lcp = std::max(local_max_lcp, plcp_lower_bound + lcp_delta);
      }
    }

    prev_sa = (std::uint64_t)sa_buffer[local_buf_filled - 1];
    sa_items_read += local_buf_filled;
#else

    // Process buffer.
    for (std::uint64_t i = 0; i < local_buf_filled; ++i) {
      std::uint64_t currsa = sa_buffer[i];
      std::uint64_t prevsa = (i == 0) ?
        prev_sa : (std::uint64_t)sa_buffer[i - 1];

      if (prevsa == text_length) {
        lcp_writer->write(0);
        continue;
      }

      // Find lower bound on PLCP[i] using sparse_plcp.
      std::uint64_t sampled_plcp_ptr = addr_buffer[i];
      std::uint64_t sampled_plcp_idx = sampled_plcp_ptr * plcp_sampling_rate;
      std::uint64_t plcp_val = plcp_buffer[2 * i];
      std::uint64_t sample_dist = currsa - sampled_plcp_idx;
      std::uint64_t plcp_lower_bound =
        std::max(0L, (std::int64_t)plcp_val - (std::int64_t)sample_dist);
      std::uint64_t pos1 = currsa + plcp_lower_bound;
      std::uint64_t pos2 = prevsa + plcp_lower_bound;
      if (pos1 > pos2) std::swap(pos1, pos2);
      if (pos2 >= text_length) {
        lcp_writer->write(plcp_lower_bound);
        local_lcp_sum += plcp_lower_bound;
        local_max_lcp = std::max(local_max_lcp, plcp_lower_bound);
        continue;
      }

      // Compute the maximal LCP delta in
      // case PLCP[currsa] was not sampled.
      std::uint64_t max_lcp_delta = 0;
      if (sampled_plcp_idx != currsa) {
        std::uint64_t max_lcp_delta2 = text_length -
          (std::max(currsa, prevsa) + plcp_lower_bound);

        if (sampled_plcp_idx + plcp_sampling_rate < text_length) {
          std::uint64_t dist = plcp_sampling_rate - sample_dist;
          max_lcp_delta = std::min(max_lcp_delta2,
              (plcp_buffer[2 * i + 1] + dist) - plcp_lower_bound);
        } else max_lcp_delta = max_lcp_delta2;
      }

      std::uint64_t lcp_delta = 0;
      if (max_lcp_delta > 0) {

        // Compute lcp_delta.
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

          // We have to read all pairs that were produced for this
          // chain and simultaneously infer the actual value of lcp_delta
          // even if the LCP is determined in one of the earlier segments.
          bool lcp_determined = false;
          while (max_lcp_delta > seg_maxlcp) {
            std::uint64_t local_lcp_delta =
              pair_multireader->read_from_ith_file(
                  file_id[seg_1_idx][seg_2_idx]);

            if (lcp_determined == false) {
              lcp_delta += local_lcp_delta;
              if (local_lcp_delta < seg_maxlcp)
                lcp_determined = true;
            }

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

          std::uint64_t local_lcp_delta =
            pair_multireader->read_from_ith_file(
                file_id[seg_1_idx][seg_2_idx]);

          if (lcp_determined == false)
            lcp_delta += local_lcp_delta;

        } else {
          lcp_delta = pair_multireader->read_from_ith_file(
                file_id[seg_1_idx][seg_2_idx]);
        }
      }

      lcp_writer->write(plcp_lower_bound + lcp_delta);
      local_lcp_sum += plcp_lower_bound + lcp_delta;
      local_max_lcp = std::max(local_max_lcp, plcp_lower_bound + lcp_delta);
    }

    prev_sa = (std::uint64_t)sa_buffer[local_buf_filled - 1];
    sa_items_read += local_buf_filled;
#endif
  }

  // Stop I/O threads.
  sa_reader->stop_reading();
  pair_multireader->stop_reading();

  // Update I/O volume.
  std::uint64_t io_vol =
    sa_reader->bytes_read() +
    lcp_writer->bytes_written() +
    pair_multireader->bytes_read();
  total_io_volume += io_vol;

  // Print summary.
  long double total_time = utils::wclock() - start;
  fprintf(stderr, "\rCompute final LCP array: time = %.1Lfs, "
      "I/O = %.1LfMiB/s, total I/O vol = %.2Lfbytes/input symbol\n",
      total_time, (io_vol / (1L << 20)) / total_time,
      (1.L * total_io_volume) / text_length);

  // Update reference variables.
  max_lcp = local_max_lcp;
  lcp_sum = local_lcp_sum;

  // Clean up.
  delete lcp_writer;
  delete sa_reader;
  delete pair_multireader;

  for (std::uint64_t i = 0; i < n_halfsegments; ++i) {
    for (std::uint64_t j = i; j < n_halfsegments; ++j) {
      std::string filename = delta_filenames[i][j];
      if (utils::file_exists(filename))
        utils::file_delete(filename);
    }
  }

  for (std::uint64_t i = 0; i < n_halfsegments; ++i)
    delete[] file_id[i];
  delete[] file_id;

  utils::deallocate(plcp_buffer);
  utils::deallocate(addr_buffer);
  utils::deallocate(sa_buffer);

#ifdef _OPENMP
  utils::deallocate(item_buffer);
#endif
}

// Compute the final LCP array. Return
// max and sum of LCP values.
template<typename text_offset_type>
void compute_final_lcp_inplace(
    std::string sa_filename,
    std::string delta_filename,
    std::string output_filename,
    std::uint64_t text_length,
    text_offset_type *sparse_plcp,
    std::uint64_t plcp_sampling_rate,
    std::uint64_t in_buf_ram,
    std::uint64_t out_buf_ram,
    std::uint64_t local_buf_ram,
    std::uint64_t &lcp_sum,
    std::uint64_t &max_lcp,
    std::uint64_t &total_io_volume) {

  // Compute size of sparse PLCP.
  std::uint64_t sparse_plcp_size =
    (text_length + plcp_sampling_rate - 1) / plcp_sampling_rate;

  // Initialize local stats.
  std::uint64_t local_lcp_sum = 0;
  std::uint64_t local_max_lcp = 0;

  // Print initial message and start the timer.
  fprintf(stderr, "Compute final LCP array: ");
  long double start = utils::wclock();

#ifdef _OPENMP
  typedef ext_buf_item_2 buf_item_type;
  std::uint64_t local_buf_size = local_buf_ram /
    (1 * sizeof(buf_item_type) +
     1 * sizeof(text_offset_type) +
     3 * sizeof(std::uint64_t));
#else
  std::uint64_t local_buf_size = local_buf_ram /
    (1 * sizeof(text_offset_type) +
     3 * sizeof(std::uint64_t));
#endif

  // Allocate buffers.
#ifdef _OPENMP
  buf_item_type *item_buffer =
    utils::allocate_array<buf_item_type>(local_buf_size);
#endif

  text_offset_type *sa_buffer =
    utils::allocate_array<text_offset_type>(local_buf_size);
  std::uint64_t *addr_buffer =
    utils::allocate_array<std::uint64_t>(local_buf_size);
  std::uint64_t *plcp_buffer =
    utils::allocate_array<std::uint64_t>(local_buf_size * 2);

  // Initialize SA reader.
  typedef async_stream_reader<text_offset_type> sa_reader_type;
  sa_reader_type *sa_reader = new sa_reader_type(
      sa_filename, in_buf_ram / 2,
      std::max(4UL, (in_buf_ram / 2) / (2UL << 20)));

  // Initialize the reader of LCP delta values.
  typedef async_stream_vbyte_reader_multipart lcp_delta_reader_type;
  lcp_delta_reader_type *lcp_delta_reader = new lcp_delta_reader_type(
        delta_filename, in_buf_ram / 2,
        std::max(4UL, (in_buf_ram / 2) / (2UL << 20)));

  // Initialize writer of the final LCP array.
  typedef async_stream_writer<text_offset_type> lcp_writer_type;
  lcp_writer_type *lcp_writer = new lcp_writer_type(
      output_filename, out_buf_ram,
      std::max(4UL, out_buf_ram / (2UL << 20)));

  // Scan suffix array and with the help
  // of PLCP array compute all LCP values.
  std::uint64_t prev_sa = text_length;
  std::uint64_t sa_items_read = 0;
  std::uint64_t step_counter = 0;
  while (sa_items_read < text_length) {
    ++step_counter;
    if (step_counter == (1UL << 25) / local_buf_size) {
      step_counter = 0;
      long double elapsed = utils::wclock() - start;
      fprintf(stderr, "\rCompute final LCP array: %.1Lf%%, time = %.1Lfs",
          (100.L * sa_items_read) / text_length, elapsed);
    }

    // Fill in the buffer.
    std::uint64_t local_buf_filled =
      std::min(text_length - sa_items_read, local_buf_size);
    sa_reader->read(sa_buffer, local_buf_filled);

#ifdef _OPENMP
#ifdef PRODUCTION_VER
    #pragma omp parallel for
    for (std::uint64_t j = 0; j < local_buf_filled; ++j)
      addr_buffer[j] = (std::uint64_t)sa_buffer[j] / plcp_sampling_rate;

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
      std::uint64_t addr = (std::uint64_t)sa_buffer[j] / plcp_sampling_rate;
      addr_buffer[j] = addr;
      plcp_buffer[2 * j] = sparse_plcp[addr];
      if (addr + 1 < sparse_plcp_size)
        plcp_buffer[2 * j + 1] = sparse_plcp[addr + 1];
    }
#endif
#else
    for (std::uint64_t j = 0; j < local_buf_filled; ++j)
      addr_buffer[j] = (std::uint64_t)sa_buffer[j] / plcp_sampling_rate;
    for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
      std::uint64_t addr = addr_buffer[j];
      plcp_buffer[2 * j] = sparse_plcp[addr];
      if (addr + 1 < sparse_plcp_size)
        plcp_buffer[2 * j + 1] = sparse_plcp[addr + 1];
    }
#endif

#ifdef _OPENMP

    // Process buffer.
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

      for (std::uint64_t i = range_beg; i < range_end; ++i) {
        std::uint64_t currsa = sa_buffer[i];
        std::uint64_t prevsa =
          (i == 0) ? prev_sa : (std::uint64_t)sa_buffer[i - 1];

        if (prevsa == text_length) {
          item_buffer[i] = buf_item_type(0, 0);
          continue;
        }

        // Initialize positions.
        std::uint64_t pos1 = currsa;
        std::uint64_t pos2 = prevsa;

        if (pos1 > pos2)
          std::swap(pos1, pos2);

        // Find lower bound on PLCP[i] using sparse_plcp.
        std::uint64_t sampled_plcp_ptr = addr_buffer[i];
        std::uint64_t sampled_plcp_idx =
          sampled_plcp_ptr * plcp_sampling_rate;
        std::uint64_t plcp_val = plcp_buffer[2 * i];
        std::uint64_t sample_dist = currsa - sampled_plcp_idx;
        std::uint64_t plcp_lower_bound =
          std::max(0L, (std::int64_t)plcp_val - (std::int64_t)sample_dist);

        pos1 += plcp_lower_bound;
        pos2 += plcp_lower_bound;

        // Compute the maximal LCP delta in
        // case PLCP[currsa] was not sampled.
        std::uint64_t max_lcp_delta = 0;
        if (sampled_plcp_idx != currsa) {
          std::uint64_t max_lcp_delta2 = text_length -
            (std::max(currsa, prevsa) + plcp_lower_bound);

          if (sampled_plcp_idx + plcp_sampling_rate < text_length) {
            std::uint64_t dist = plcp_sampling_rate - sample_dist;
            max_lcp_delta = std::min(max_lcp_delta2,
                (plcp_buffer[2 * i + 1] + dist) - plcp_lower_bound);
          } else max_lcp_delta = max_lcp_delta2;
        }

        if (max_lcp_delta == 0)
          item_buffer[i] = buf_item_type(plcp_lower_bound, 0);
        else item_buffer[i] = buf_item_type(plcp_lower_bound, 1);
      }
    }

    for (std::uint64_t i = 0; i < local_buf_filled; ++i) {
      if (item_buffer[i].m_type == 0) {
        std::uint64_t plcp_lower_bound = item_buffer[i].m_plcp_lower_bound;
        lcp_writer->write(plcp_lower_bound);
        local_lcp_sum += plcp_lower_bound;
        local_max_lcp = std::max(local_max_lcp, plcp_lower_bound);
      } else {
        std::uint64_t plcp_lower_bound = item_buffer[i].m_plcp_lower_bound;
        std::uint64_t lcp_delta = lcp_delta_reader->read();
        lcp_writer->write(plcp_lower_bound + lcp_delta);
        local_lcp_sum += plcp_lower_bound + lcp_delta;
        local_max_lcp = std::max(local_max_lcp, plcp_lower_bound + lcp_delta);
      }
    }

    prev_sa = (std::uint64_t)sa_buffer[local_buf_filled - 1];
    sa_items_read += local_buf_filled;
#else

    // Process buffer.
    for (std::uint64_t i = 0; i < local_buf_filled; ++i) {
      std::uint64_t currsa = sa_buffer[i];
      std::uint64_t prevsa = (i == 0) ?
        prev_sa : (std::uint64_t)sa_buffer[i - 1];

      if (prevsa == text_length) {
        lcp_writer->write(0);
        continue;
      }

      std::uint64_t pos1 = currsa;
      std::uint64_t pos2 = prevsa;

      if (pos1 > pos2)
        std::swap(pos1, pos2);

      // Find lower bound on PLCP[i] using sparse_plcp.
      std::uint64_t sampled_plcp_ptr = addr_buffer[i];
      std::uint64_t sampled_plcp_idx = sampled_plcp_ptr * plcp_sampling_rate;
      std::uint64_t plcp_val = plcp_buffer[2 * i];
      std::uint64_t sample_dist = currsa - sampled_plcp_idx;
      std::uint64_t plcp_lower_bound =
        std::max(0L, (std::int64_t)plcp_val - (std::int64_t)sample_dist);

      pos1 += plcp_lower_bound;
      pos2 += plcp_lower_bound;

      if (pos2 >= text_length) {
        lcp_writer->write(plcp_lower_bound);
        local_lcp_sum += plcp_lower_bound;
        local_max_lcp = std::max(local_max_lcp, plcp_lower_bound);
        continue;
      }

      // Compute the maximal LCP delta in
      // case PLCP[currsa] was not sampled.
      std::uint64_t max_lcp_delta = 0;
      if (sampled_plcp_idx != currsa) {
        std::uint64_t max_lcp_delta2 = text_length -
          (std::max(currsa, prevsa) + plcp_lower_bound);

        if (sampled_plcp_idx + plcp_sampling_rate < text_length) {
          std::uint64_t dist = plcp_sampling_rate - sample_dist;
          max_lcp_delta = std::min(max_lcp_delta2,
              (plcp_buffer[2 * i + 1] + dist) - plcp_lower_bound);
        } else max_lcp_delta = max_lcp_delta2;
      }

      // Compute LCP delta.
      std::uint64_t lcp_delta = 0;
      if (max_lcp_delta > 0)
        lcp_delta = lcp_delta_reader->read();

      // Write final LCP value to output
      // file and update statistics.
      lcp_writer->write(plcp_lower_bound + lcp_delta);
      local_lcp_sum += plcp_lower_bound + lcp_delta;
      local_max_lcp = std::max(local_max_lcp, plcp_lower_bound + lcp_delta);
    }

    prev_sa = (std::uint64_t)sa_buffer[local_buf_filled - 1];
    sa_items_read += local_buf_filled;
#endif
  }

  // Stop I/O threads.
  sa_reader->stop_reading();
  lcp_delta_reader->stop_reading();

  // Update I/O volume.
  std::uint64_t io_vol =
    sa_reader->bytes_read() +
    lcp_delta_reader->bytes_read() +
    lcp_writer->bytes_written();
  total_io_volume += io_vol;

  // Print summary.
  long double total_time = utils::wclock() - start;
  fprintf(stderr, "\rCompute final LCP array: time = %.1Lfs, "
      "I/O = %.1LfMiB/s, total I/O vol = %.2Lfbytes/input symbol\n",
      total_time, (io_vol / (1L << 20)) / total_time,
      (1.L * total_io_volume) / text_length);

  // Update reference variables.
  max_lcp = local_max_lcp;
  lcp_sum = local_lcp_sum;

  // Clean up.
  delete lcp_writer;
  delete lcp_delta_reader;
  delete sa_reader;

  utils::deallocate(plcp_buffer);
  utils::deallocate(addr_buffer);
  utils::deallocate(sa_buffer);

#ifdef _OPENMP
  utils::deallocate(item_buffer);
#endif
}

}  // namespace em_sparse_phi_private

#endif  // __SRC_EM_SPARSE_PHI_SRC_COMPUTE_FINAL_LCP_HPP_INCLUDED
