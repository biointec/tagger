/**
 * @file    src/em_sparse_phi_src/preprocess_lcp_delta_values.hpp
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

#ifndef __SRC_EM_SPARSE_PHI_SRC_PREPROCESS_LCP_DELTA_VALUES_HPP_INCLUDED
#define __SRC_EM_SPARSE_PHI_SRC_PREPROCESS_LCP_DELTA_VALUES_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <algorithm>
#include <omp.h>

#include "io/async_stream_reader.hpp"
#include "io/async_stream_writer.hpp"
#include "io/async_stream_writer_multipart.hpp"
#include "io/async_multi_stream_vbyte_reader.hpp"
#include "utils.hpp"
#include "convert_to_vbyte_slab.hpp"


namespace em_sparse_phi_private {

template<typename T>
struct ext_buf_item_4 {
  ext_buf_item_4() {}
  ext_buf_item_4(
      T idx1,
      T idx2,
      std::uint8_t tp)
    : m_seg_1_idx(idx1),
    m_seg_2_idx(idx2),
    m_type(tp) {}

  T m_seg_1_idx;
  T m_seg_2_idx;
  std::uint8_t m_type;
} __attribute__ ((packed));

template<typename T>
struct ext_buf_item_3 {
  ext_buf_item_3() {}
  ext_buf_item_3(
      T idx1,
      T idx2,
      std::uint8_t tp,
      T orig_idx1,
      T orig_idx2)
    : m_seg_1_idx(idx1),
    m_seg_2_idx(idx2),
    m_type(tp),
    m_pos1_orig_halfseg_idx(orig_idx1),
    m_pos2_orig_halfseg_idx(orig_idx2) {}

  T m_seg_1_idx;
  T m_seg_2_idx;
  std::uint8_t m_type;
  T m_pos1_orig_halfseg_idx;
  T m_pos2_orig_halfseg_idx;
} __attribute__ ((packed));

template<typename text_offset_type>
void preprocess_lcp_delta_values_lex_partitioning(
    std::string sa_filename,
    std::string lcp_deltas_filename,
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
    std::uint64_t &total_io_volume) {

  std::uint64_t n_halfsegments =
    (text_length + max_halfsegment_size - 1) / max_halfsegment_size;
  std::uint64_t sparse_plcp_size =
    (text_length + plcp_sampling_rate - 1) / plcp_sampling_rate;

  fprintf(stderr, "\nPreprocess LCP delta values: ");
  long double start = utils::wclock();
  std::uint64_t io_volume = 0;

#ifdef _OPENMP
  typedef ext_buf_item_4<std::uint64_t> buf_item_type;
  std::uint64_t local_buf_size = local_buf_ram /
    (1 * sizeof(buf_item_type) +
     1 * sizeof(text_offset_type) +
     9 * sizeof(std::uint8_t) +
     4 * sizeof(std::uint64_t));
#else
  std::uint64_t local_buf_size = local_buf_ram /
    (1 * sizeof(text_offset_type) +
     9 * sizeof(std::uint8_t) +
     4 * sizeof(std::uint64_t));
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
  std::uint64_t *lcp_deltas =
    utils::allocate_array<std::uint64_t>(local_buf_size);
  std::uint8_t *vbyte_slab =
    utils::allocate_array<std::uint8_t>(local_buf_size * 9);

  // Shrink the buffer size to accommodate
  // extra files from text partitioning.
  std::uint64_t different_halfsegment_pairs_all_parts = 0;
  for (std::uint64_t i = 0; i < n_halfsegments; ++i) {
    for (std::uint64_t j = i; j < n_halfsegments; ++j) {
      std::string filename = delta_filenames[i][j];
      if (utils::file_exists(filename) && utils::file_size(filename) > 0)
        ++different_halfsegment_pairs_all_parts;
    }
  }

  // Create multireader for halfsegment pairs in all parts.
  typedef async_multi_stream_vbyte_reader pair_multireader_type;
  std::uint64_t buffer_size = halfseg_buffers_ram /
    std::max(1UL, different_halfsegment_pairs_all_parts);
  pair_multireader_type *pair_multireader = new pair_multireader_type(
      different_halfsegment_pairs_all_parts, buffer_size);

  // Initialize SA reader.
  typedef async_stream_reader<text_offset_type> sa_reader_type;
  sa_reader_type *sa_reader = new sa_reader_type(
      sa_filename, in_buf_ram,
      std::max(4UL, in_buf_ram / (2UL << 20)));

  // Initialize writer of the final LCP array.
  typedef async_stream_writer_multipart<std::uint8_t> lcp_delta_writer_type;
  lcp_delta_writer_type *lcp_delta_writer = new lcp_delta_writer_type(
      lcp_deltas_filename, out_buf_ram,
      std::max(4UL, out_buf_ram / (2UL << 20)));

  // Compute file ID array for the multireader
  // and the total size of all delta files.
  std::uint64_t total_lcp_delta_file_sizes = 0;
  std::uint64_t **file_id = new std::uint64_t*[n_halfsegments];
  {
    std::uint64_t file_id_counter = 0;
    for (std::uint64_t i = 0; i < n_halfsegments; ++i) {
      file_id[i] = new std::uint64_t[n_halfsegments];
      for (std::uint64_t j = i; j < n_halfsegments; ++j) {
        std::string filename = delta_filenames[i][j];
        if (utils::file_exists(filename)) {
          std::uint64_t filesize = utils::file_size(filename);
          if (filesize > 0) {
            file_id[i][j] = file_id_counter++;
            pair_multireader->add_file(filename);
            total_lcp_delta_file_sizes += filesize;
          }
        }
      }
    }
  }

  static const std::uint64_t max_allowed_extra_disk_space = (1UL << 20);

  std::uint64_t prev_sa = text_length;
  std::uint64_t sa_items_read = 0;
  std::uint64_t step_counter = 0;
  std::uint64_t used_disk_space = total_lcp_delta_file_sizes;
  std::uint64_t cur_delta_file_size = 0;
  std::uint64_t lcp_deltas_filled = 0;

  // Rewrite the LCP delta values into a single
  // sequence (but stored in multiple files).
  while (sa_items_read < text_length) {
    ++step_counter;
    if (step_counter == (1UL << 25) / local_buf_size) {
      step_counter = 0;
      long double elapsed = utils::wclock() - start;
      fprintf(stderr, "\rPreprocess LCP delta "
          "values: %.1Lf%%, time = %.1Lfs",
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

      std::uint64_t prev_src_halfseg_id = ((range_beg == 0) ? prev_sa :
          (std::uint64_t)sa_buffer[range_beg - 1]) / max_halfsegment_size;

      for (std::uint64_t i = range_beg; i < range_end; ++i) {
        std::uint64_t currsa = sa_buffer[i];
        std::uint64_t prevsa = (i == 0) ?
          prev_sa : (std::uint64_t)sa_buffer[i - 1];

        // Compute the 'source' halfsegment IDs h1 and h2.
        std::uint64_t cur_src_halfseg_id = currsa / max_halfsegment_size;
        std::uint64_t h1 = cur_src_halfseg_id;
        std::uint64_t h2 = prev_src_halfseg_id;
        prev_src_halfseg_id = cur_src_halfseg_id;

        if (prevsa == text_length) {
          item_buffer[i] = buf_item_type(0, 0, 0);
          continue;
        }

        // Initialize positions.
        std::uint64_t pos1 = currsa;
        std::uint64_t pos2 = prevsa;
        if (pos1 > pos2) {
          std::swap(pos1, pos2);
          std::swap(h1, h2);
        }

        // Find lower bound on PLCP[i] using sparse_plcp.
        std::uint64_t sampled_plcp_ptr = addr_buffer[i];
        std::uint64_t sampled_plcp_idx =
          sampled_plcp_ptr * plcp_sampling_rate;
        std::uint64_t prev_sample_dist = currsa - sampled_plcp_idx;
        std::uint64_t next_sample_dist =
          plcp_sampling_rate - prev_sample_dist;
        std::uint64_t plcp_val = plcp_buffer[2 * i];
        std::uint64_t plcp_lower_bound = std::max(0L,
            (std::int64_t)plcp_val - (std::int64_t)prev_sample_dist);

        pos1 += plcp_lower_bound;
        pos2 += plcp_lower_bound;

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

        if (max_lcp_delta == 0)
          item_buffer[i] = buf_item_type(0, 0, 0);
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
            item_buffer[i] = buf_item_type(0, 0, 2);
          else
            item_buffer[i] = buf_item_type(seg_1_idx, seg_2_idx, 1);
        }
      }
    }

    for (std::uint64_t i = 0; i < local_buf_filled; ++i) {
      if (item_buffer[i].m_type == 1) {

        // Common case.
        std::uint64_t seg_1_idx = item_buffer[i].m_seg_1_idx;
        std::uint64_t seg_2_idx = item_buffer[i].m_seg_2_idx;
        std::uint64_t lcp_delta =
          pair_multireader->read_from_ith_file(
            file_id[seg_1_idx][seg_2_idx]);
        lcp_deltas[lcp_deltas_filled++] = lcp_delta;
      } else if (item_buffer[i].m_type == 2) {

        // Recompute everything to decode the chain of pairs.
        std::uint64_t currsa = sa_buffer[i];
        std::uint64_t prevsa = (i == 0) ?
          prev_sa : (std::uint64_t)sa_buffer[i - 1];

        std::uint64_t pos1 = currsa;
        std::uint64_t pos2 = prevsa;
        if (pos1 > pos2)
          std::swap(pos1, pos2);

        // Find lower bound on PLCP[i] using sparse_plcp.
        std::uint64_t sampled_plcp_ptr = addr_buffer[i];
        std::uint64_t sampled_plcp_idx =
          sampled_plcp_ptr * plcp_sampling_rate;
        std::uint64_t prev_sample_dist = currsa - sampled_plcp_idx;
        std::uint64_t next_sample_dist =
          plcp_sampling_rate - prev_sample_dist;
        std::uint64_t plcp_val = plcp_buffer[2 * i];
        std::uint64_t plcp_lower_bound = std::max(0L,
            (std::int64_t)plcp_val - (std::int64_t)prev_sample_dist);

        pos1 += plcp_lower_bound;
        pos2 += plcp_lower_bound;

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
        // chain and simultnaously infer the actual value of lcp_delta
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

        // Add LCP delta to the buffer.
        lcp_deltas[lcp_deltas_filled++] = lcp_delta;
      }

      // XXX 32-bit and 48-bit integers?
      used_disk_space += sizeof(text_offset_type);
      if (lcp_deltas_filled == local_buf_size ||
          used_disk_space > text_length * sizeof(text_offset_type) +
          max_allowed_extra_disk_space) {

        // Convert lcp deltas to vbyte.
        std::uint64_t vbyte_slab_length =
          convert_to_vbyte_slab(lcp_deltas,
              lcp_deltas_filled, vbyte_slab);
        cur_delta_file_size += vbyte_slab_length;
        lcp_deltas_filled = 0;

        // Write lcp deltas to disk.
        lcp_delta_writer->write(vbyte_slab, vbyte_slab_length);

        // Check if need to end the current file.
        if (used_disk_space > text_length * sizeof(text_offset_type) +
            max_allowed_extra_disk_space) {
          lcp_delta_writer->end_current_file();

          // Simulate the space assuming we
          // now delete the current delta file.
          used_disk_space -= cur_delta_file_size;
          cur_delta_file_size = 0;
        }
      }
    }

    // Update local variables.
    prev_sa = (std::uint64_t)sa_buffer[local_buf_filled - 1];
    sa_items_read += local_buf_filled;

#else

    // Process buffer.
    std::uint64_t prev_src_halfseg_id = prev_sa / max_halfsegment_size;
    for (std::uint64_t i = 0; i < local_buf_filled; ++i) {
      std::uint64_t currsa = sa_buffer[i];
      std::uint64_t prevsa = (i == 0) ?
        prev_sa : (std::uint64_t)sa_buffer[i - 1];

      // Compute the 'source' halfsegment IDs h1 and h2.
      std::uint64_t cur_src_halfseg_id = currsa / max_halfsegment_size;
      std::uint64_t h1 = cur_src_halfseg_id;
      std::uint64_t h2 = prev_src_halfseg_id;
      prev_src_halfseg_id = cur_src_halfseg_id;

      if (prevsa == text_length) {
        used_disk_space += sizeof(text_offset_type);  // XXX fix this.
        continue;
      }

      std::uint64_t pos1 = currsa;
      std::uint64_t pos2 = prevsa;

      if (pos1 > pos2) {
        std::swap(pos1, pos2);
        std::swap(h1, h2);
      }

      // Find lower bound on PLCP[i] using sparse_plcp.
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

      // Compute the maximal LCP delta in case PLCP[currsa] was not sampled.
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

        // Add LCP delta value to a biffer.
        lcp_deltas[lcp_deltas_filled++] = lcp_delta;
      }

      // XXX 32-bit and 48-bit integers?
      used_disk_space += sizeof(text_offset_type);
      if (lcp_deltas_filled == local_buf_size ||
          used_disk_space > text_length * sizeof(text_offset_type) +
          max_allowed_extra_disk_space) {

        // Convert lcp deltas to vbyte.
        std::uint64_t vbyte_slab_length =
          convert_to_vbyte_slab(lcp_deltas,
              lcp_deltas_filled, vbyte_slab);
        cur_delta_file_size += vbyte_slab_length;
        lcp_deltas_filled = 0;

        // Write lcp deltas to disk.
        lcp_delta_writer->write(vbyte_slab, vbyte_slab_length);

        // Check if need to end the current file.
        if (used_disk_space > text_length * sizeof(text_offset_type) +
            max_allowed_extra_disk_space) {
          lcp_delta_writer->end_current_file();

          // Simulate the space assuming we
          // now delete the current delta file.
          used_disk_space -= cur_delta_file_size;
          cur_delta_file_size = 0;
        }
      }
    }

    // Update local variables.
    prev_sa = (std::uint64_t)sa_buffer[local_buf_filled - 1];
    sa_items_read += local_buf_filled;
#endif
  }

  if (lcp_deltas_filled > 0) {
    // Convert lcp deltas to vbyte.
    std::uint64_t vbyte_slab_length =
      convert_to_vbyte_slab(lcp_deltas, lcp_deltas_filled, vbyte_slab);

    // Write lcp deltas to disk.
    lcp_delta_writer->write(vbyte_slab, vbyte_slab_length);
  }

  // Stop I/O threads.
  sa_reader->stop_reading();
  pair_multireader->stop_reading();

  // Update I/O volume.
  io_volume +=
    sa_reader->bytes_read() +
    lcp_delta_writer->bytes_written() +
    pair_multireader->bytes_read();
  total_io_volume += io_volume;

  // Print summary.
  long double total_time = utils::wclock() - start;
  fprintf(stderr, "\rPreprocess LCP delta values: time = %.1Lfs, "
      "I/O = %.1LfMiB/s, total I/O vol = %.2Lfbytes/input symbol\n",
      total_time, (io_volume / (1L << 20)) / total_time,
      (1.L * total_io_volume) / text_length);

  // Clean up.
  for (std::uint64_t i = 0; i < n_halfsegments; ++i)
    delete[] file_id[i];
  delete[] file_id;

  delete lcp_delta_writer;
  delete sa_reader;
  delete pair_multireader;

  utils::deallocate(vbyte_slab);
  utils::deallocate(lcp_deltas);
  utils::deallocate(plcp_buffer);
  utils::deallocate(addr_buffer);
  utils::deallocate(sa_buffer);

#ifdef _OPENMP
  utils::deallocate(item_buffer);
#endif

  for (std::uint64_t i = 0; i < n_halfsegments; ++i) {
    for (std::uint64_t j = i; j < n_halfsegments; ++j) {
      std::string filename = delta_filenames[i][j];
      if (utils::file_exists(filename))
        utils::file_delete(filename);
    }
  }
}

template<typename text_offset_type>
void preprocess_lcp_delta_values_text_partitioning(
    std::string sa_filename,
    std::string lcp_deltas_filename,
    std::uint64_t text_length,
    text_offset_type *sparse_plcp,
    std::uint64_t plcp_sampling_rate,
    std::uint64_t max_halfsegment_size,
    std::uint64_t halfseg_buffers_ram,
    std::uint64_t in_buf_ram,
    std::uint64_t out_buf_ram,
    std::uint64_t local_buf_ram,
    std::uint64_t max_overflow_size,
    std::uint64_t n_parts,
    std::uint64_t ***items_per_halfseg_pair,
    std::string ***delta_filenames,
    std::uint64_t &total_io_volume) {

  // Compute basic parameters.
  std::uint64_t n_halfsegments =
    (text_length + max_halfsegment_size - 1) / max_halfsegment_size;
  std::uint64_t sparse_plcp_size =
    (text_length + plcp_sampling_rate - 1) / plcp_sampling_rate;

  // Print initial message and start timer.
  fprintf(stderr, "\nPreprocess LCP delta values: ");
  long double start = utils::wclock();
  std::uint64_t io_volume = 0;

#ifdef _OPENMP
  typedef ext_buf_item_3<std::uint64_t> buf_item_type;
  std::uint64_t local_buf_size = local_buf_ram /
    (1 * sizeof(buf_item_type) +
     1 * sizeof(text_offset_type) +
     9 * sizeof(std::uint8_t) +
     4 * sizeof(std::uint64_t));
#else
  std::uint64_t local_buf_size = local_buf_ram /
    (1 * sizeof(text_offset_type) +
     9 * sizeof(std::uint8_t) +
     4 * sizeof(std::uint64_t));
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
  std::uint64_t *lcp_deltas =
    utils::allocate_array<std::uint64_t>(local_buf_size);
  std::uint8_t *vbyte_slab =
    utils::allocate_array<std::uint8_t>(local_buf_size * 9);

  // Shrink the buffer size to accommodate
  // extra files from text partitioning.
  std::uint64_t different_halfsegment_pairs_all_parts = 0;
  for (std::uint64_t i = 0; i < n_halfsegments; ++i) {
    for (std::uint64_t j = i; j < n_halfsegments; ++j) {
      for (std::uint64_t part_id = 0; part_id < n_parts; ++part_id) {
        std::string filename = delta_filenames[part_id][i][j];
        if (utils::file_exists(filename) && utils::file_size(filename) > 0)
          ++different_halfsegment_pairs_all_parts;
      }
    }
  }

  // Create multireader for halfsegment pairs in all parts.
  typedef async_multi_stream_vbyte_reader pair_multireader_type;
  std::uint64_t buffer_size = halfseg_buffers_ram /
    std::max(1UL, different_halfsegment_pairs_all_parts);
  pair_multireader_type *pair_multireader = new pair_multireader_type(
      different_halfsegment_pairs_all_parts, buffer_size);

  // Initialize SA reader.
  typedef async_stream_reader<text_offset_type> sa_reader_type;
  sa_reader_type *sa_reader = new sa_reader_type(
      sa_filename, in_buf_ram,
      std::max(4UL, in_buf_ram / (2UL << 20)));

  // Initialize writer of the final LCP array.
  typedef async_stream_writer_multipart<std::uint8_t> lcp_delta_writer_type;
  lcp_delta_writer_type *lcp_delta_writer = new lcp_delta_writer_type(
      lcp_deltas_filename, out_buf_ram,
      std::max(4UL, out_buf_ram / (2UL << 20)));

  // Compute file ID array for the multireader
  // and the total size of all delta files.
  std::uint64_t total_lcp_delta_file_sizes = 0;
  std::uint64_t ***file_id = new std::uint64_t**[n_parts];
  {
    std::uint64_t file_id_counter = 0;
    for (std::uint64_t part_id = 0; part_id < n_parts; ++part_id) {
      file_id[part_id] = new std::uint64_t*[n_halfsegments];
      for (std::uint64_t i = 0; i < n_halfsegments; ++i) {
        file_id[part_id][i] = new std::uint64_t[n_halfsegments];
        for (std::uint64_t j = i; j < n_halfsegments; ++j) {
          std::string filename = delta_filenames[part_id][i][j];
          if (utils::file_exists(filename)) {
            std::uint64_t filesize = utils::file_size(filename);
            if (filesize > 0) {
              file_id[part_id][i][j] = file_id_counter++;
              pair_multireader->add_file(filename);
              total_lcp_delta_file_sizes += filesize;
            }
          }
        }
      }
    }
  }

  // Initialize the array that tells, for every pairs of halfsegments,
  // in which part is the next element (from that pair of halfsegments)
  // processed. This array is updated as we scan the suffix array.
  std::uint64_t **current_part = new std::uint64_t*[n_halfsegments];
  for (std::uint64_t i = 0; i < n_halfsegments; ++i) {
    current_part[i] = new std::uint64_t[n_halfsegments];
    for (std::uint64_t j = i; j < n_halfsegments; ++j) {
      current_part[i][j] = 0;
      while (current_part[i][j] < n_parts &&
          items_per_halfseg_pair[current_part[i][j]][i][j] == 0)
        ++current_part[i][j];
    }
  }

  static const std::uint64_t max_allowed_extra_disk_space = (1UL << 20);

  std::uint64_t prev_sa = text_length;
  std::uint64_t sa_items_read = 0;
  std::uint64_t step_counter = 0;
  std::uint64_t used_disk_space = total_lcp_delta_file_sizes;
  std::uint64_t cur_delta_file_size = 0;
  std::uint64_t lcp_deltas_filled = 0;

  // Rewrite the LCP delta values into a single
  // sequence (but stored in multiple files).
  while (sa_items_read < text_length) {
    ++step_counter;
    if (step_counter == (1UL << 25) / local_buf_size) {
      step_counter = 0;
      long double elapsed = utils::wclock() - start;
      fprintf(stderr, "\rPreprocess LCP delta values: "
          "%.1Lf%%, time = %.1Lfs",
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

      std::uint64_t prev_src_halfseg_id = ((range_beg == 0) ? prev_sa :
          (std::uint64_t)sa_buffer[range_beg - 1]) / max_halfsegment_size;

      for (std::uint64_t i = range_beg; i < range_end; ++i) {
        std::uint64_t currsa = sa_buffer[i];
        std::uint64_t prevsa = (i == 0) ?
          prev_sa : (std::uint64_t)sa_buffer[i - 1];

        // Compute the 'source' halfsegment IDs h1 and h2.
        std::uint64_t cur_src_halfseg_id = currsa / max_halfsegment_size;
        std::uint64_t h1 = cur_src_halfseg_id;
        std::uint64_t h2 = prev_src_halfseg_id;
        prev_src_halfseg_id = cur_src_halfseg_id;

        if (prevsa == text_length) {
          item_buffer[i] = buf_item_type(
              0, 0, 0, n_halfsegments, n_halfsegments);
          continue;
        }

        // Initialize positions.
        std::uint64_t pos1 = currsa;
        std::uint64_t pos2 = prevsa;
        if (pos1 > pos2) {
          std::swap(pos1, pos2);
          std::swap(h1, h2);
        }

        // Find lower bound on PLCP[i] using sparse_plcp.
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

        if (max_lcp_delta == 0)
          item_buffer[i] = buf_item_type(0, 0, 0, h1, h2);
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
            item_buffer[i] = buf_item_type(0, 0, 2, h1, h2);
          else
            item_buffer[i] = buf_item_type(seg_1_idx, seg_2_idx, 1, h1, h2);
        }
      }
    }

    for (std::uint64_t i = 0; i < local_buf_filled; ++i) {

      // Compute the ID of the part in which this pair was processed.
      std::uint64_t part_id = 0;
      if ((std::uint64_t)item_buffer[i].m_pos1_orig_halfseg_idx !=
          n_halfsegments) {

        // The above if check if this was not the pair (SA[0], undefined).
        std::uint64_t h1 = item_buffer[i].m_pos1_orig_halfseg_idx;
        std::uint64_t h2 = item_buffer[i].m_pos2_orig_halfseg_idx;
        part_id = current_part[h1][h2];

        // Update current_part and items_per_halfseg_pair.
        --items_per_halfseg_pair[part_id][h1][h2];
        if (items_per_halfseg_pair[part_id][h1][h2] == 0) {
          while (current_part[h1][h2] < n_parts &&
              items_per_halfseg_pair[current_part[h1][h2]][h1][h2] == 0)
            ++current_part[h1][h2];
        }
      }

      if (item_buffer[i].m_type == 1) {

        // Common case.
        std::uint64_t seg_1_idx = item_buffer[i].m_seg_1_idx;
        std::uint64_t seg_2_idx = item_buffer[i].m_seg_2_idx;
        std::uint64_t lcp_delta =
          pair_multireader->read_from_ith_file(
            file_id[part_id][seg_1_idx][seg_2_idx]);
        lcp_deltas[lcp_deltas_filled++] = lcp_delta;
      } else if (item_buffer[i].m_type == 2) {

        // Recompute everything to decode the chain of pairs.
        std::uint64_t currsa = sa_buffer[i];
        std::uint64_t prevsa = (i == 0) ?
          prev_sa : (std::uint64_t)sa_buffer[i - 1];

        std::uint64_t pos1 = currsa;
        std::uint64_t pos2 = prevsa;
        if (pos1 > pos2)
          std::swap(pos1, pos2);

        // Find lower bound on PLCP[i] using sparse_plcp.
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
        // chain and simultnaously infer the actual value of lcp_delta
        // even if the LCP is determined in one of the earlier segments.
        bool lcp_determined = false;
        while (max_lcp_delta > seg_maxlcp) {
          std::uint64_t local_lcp_delta =
            pair_multireader->read_from_ith_file(
                file_id[part_id][seg_1_idx][seg_2_idx]);

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
              file_id[part_id][seg_1_idx][seg_2_idx]);

        if (lcp_determined == false)
          lcp_delta += local_lcp_delta;

        lcp_deltas[lcp_deltas_filled++] = lcp_delta;
      }

      // XXX 32-bit and 48-bit integers?
      used_disk_space += sizeof(text_offset_type);
      if (lcp_deltas_filled == local_buf_size ||
          used_disk_space > text_length * sizeof(text_offset_type) +
          max_allowed_extra_disk_space) {

        // Convert lcp deltas to vbyte.
        std::uint64_t vbyte_slab_length =
          convert_to_vbyte_slab(lcp_deltas,
              lcp_deltas_filled, vbyte_slab);
        cur_delta_file_size += vbyte_slab_length;
        lcp_deltas_filled = 0;

        // Write lcp deltas to disk.
        lcp_delta_writer->write(vbyte_slab, vbyte_slab_length);

        // Check if need to end the current file.
        if (used_disk_space > text_length * sizeof(text_offset_type) +
            max_allowed_extra_disk_space) {
          lcp_delta_writer->end_current_file();

          // Simulate the space assuming we
          // now delete the current delta file.
          used_disk_space -= cur_delta_file_size;
          cur_delta_file_size = 0;
        }
      }
    }

    // Update local variables.
    prev_sa = (std::uint64_t)sa_buffer[local_buf_filled - 1];
    sa_items_read += local_buf_filled;

#else

    // Process buffer.
    std::uint64_t prev_src_halfseg_id = prev_sa / max_halfsegment_size;
    for (std::uint64_t i = 0; i < local_buf_filled; ++i) {
      std::uint64_t currsa = sa_buffer[i];
      std::uint64_t prevsa = (i == 0) ?
        prev_sa : (std::uint64_t)sa_buffer[i - 1];

      // Compute the 'source' halfsegment IDs h1 and h2.
      std::uint64_t cur_src_halfseg_id = currsa / max_halfsegment_size;
      std::uint64_t h1 = cur_src_halfseg_id;
      std::uint64_t h2 = prev_src_halfseg_id;
      prev_src_halfseg_id = cur_src_halfseg_id;

      if (prevsa == text_length) {
        used_disk_space += sizeof(text_offset_type);  // XXX fix this
        continue;
      }

      std::uint64_t pos1 = currsa;
      std::uint64_t pos2 = prevsa;

      if (pos1 > pos2) {
        std::swap(pos1, pos2);
        std::swap(h1, h2);
      }

      // Invariant: pos1 <= pos2.
      // Compute the ID of the part in which this pair was processed.
      std::uint64_t part_id = 0;
      {
        part_id = current_part[h1][h2];

        // Update current_part and items_per_halfseg_pair.
        --items_per_halfseg_pair[part_id][h1][h2];
        if (items_per_halfseg_pair[part_id][h1][h2] == 0) {
          while (current_part[h1][h2] < n_parts &&
              items_per_halfseg_pair[current_part[h1][h2]][h1][h2] == 0)
            ++current_part[h1][h2];
        }
      }

      // Find lower bound on PLCP[i] using sparse_plcp.
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

      // Compute the maximal LCP delta in case PLCP[currsa] was not sampled.
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

        if (max_lcp_delta > seg_maxlcp) {

          // We have to read all pairs that were produced for this
          // chain and simultaneously infer the actual value of lcp_delta
          // even if the LCP is determined in one of the earlier segments.
          bool lcp_determined = false;
          while (max_lcp_delta > seg_maxlcp) {
            std::uint64_t local_lcp_delta =
              pair_multireader->read_from_ith_file(
                  file_id[part_id][seg_1_idx][seg_2_idx]);
            
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
                file_id[part_id][seg_1_idx][seg_2_idx]);

          if (lcp_determined == false)
            lcp_delta += local_lcp_delta;
        } else {
          lcp_delta = pair_multireader->read_from_ith_file(
              file_id[part_id][seg_1_idx][seg_2_idx]);
        }

        lcp_deltas[lcp_deltas_filled++] = lcp_delta;
      }

      // XXX 32-bit and 48-bit integers?
      used_disk_space += sizeof(text_offset_type);
      if (lcp_deltas_filled == local_buf_size ||
          used_disk_space > text_length * sizeof(text_offset_type) +
          max_allowed_extra_disk_space) {

        // Convert lcp deltas to vbyte.
        std::uint64_t vbyte_slab_length =
          convert_to_vbyte_slab(lcp_deltas,
              lcp_deltas_filled, vbyte_slab);
        cur_delta_file_size += vbyte_slab_length;
        lcp_deltas_filled = 0;

        // Write lcp deltas to disk.
        lcp_delta_writer->write(vbyte_slab, vbyte_slab_length);

        // Check if need to end the current file.
        if (used_disk_space > text_length *
            sizeof(text_offset_type) +
            max_allowed_extra_disk_space) {
          lcp_delta_writer->end_current_file();

          // Simulate the space assuming we
          // now delete the current delta file.
          used_disk_space -= cur_delta_file_size;
          cur_delta_file_size = 0;
        }
      }
    }

    // Update local variables.
    prev_sa = (std::uint64_t)sa_buffer[local_buf_filled - 1];
    sa_items_read += local_buf_filled;
#endif
  }

  if (lcp_deltas_filled > 0) {

    // Convert lcp deltas to vbyte.
    std::uint64_t vbyte_slab_length =
      convert_to_vbyte_slab(lcp_deltas, lcp_deltas_filled, vbyte_slab);

    // Write lcp deltas to disk.
    lcp_delta_writer->write(vbyte_slab, vbyte_slab_length);
  }

  // Stop I/O threads.
  sa_reader->stop_reading();
  pair_multireader->stop_reading();

  // Update I/O volume.
  io_volume +=
    sa_reader->bytes_read() +
    lcp_delta_writer->bytes_written() +
    pair_multireader->bytes_read();
  total_io_volume += io_volume;

  // Print summary.
  long double total_time = utils::wclock() - start;
  fprintf(stderr, "\rPreprocess LCP delta values: time = %.1Lfs, "
      "I/O = %.1LfMiB/s, total I/O vol = %.2Lfbytes/input symbol\n",
      total_time, (io_volume / (1L << 20)) / total_time,
      (1.L * total_io_volume) / text_length);

  // Clean up.
  for (std::uint64_t i = 0; i < n_halfsegments; ++i)
    delete[] current_part[i];
  delete[] current_part;

  for (std::uint64_t part_id = 0; part_id < n_parts; ++part_id) {
    for (std::uint64_t i = 0; i < n_halfsegments; ++i)
      delete[] file_id[part_id][i];
    delete[] file_id[part_id];
  }
  delete[] file_id;

  delete lcp_delta_writer;
  delete sa_reader;
  delete pair_multireader;

  utils::deallocate(vbyte_slab);
  utils::deallocate(lcp_deltas);
  utils::deallocate(plcp_buffer);
  utils::deallocate(addr_buffer);
  utils::deallocate(sa_buffer);

#ifdef _OPENMP
  utils::deallocate(item_buffer);
#endif

  for (std::uint64_t part_id = 0; part_id < n_parts; ++part_id) {
    for (std::uint64_t i = 0; i < n_halfsegments; ++i) {
      for (std::uint64_t j = i; j < n_halfsegments; ++j) {
        std::string filename = delta_filenames[part_id][i][j];
        if (utils::file_exists(filename))
          utils::file_delete(filename);
      }
    }
  }
}

}  // namespace em_sparse_phi_private

#endif  // __SRC_EM_SPARSE_PHI_SRC_PREPROCESS_LCP_DELTA_VALUES_HPP_INCLUDED
