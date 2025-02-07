/**
 * @file    src/em_sparse_phi_src/compute_sparse_plcp.hpp
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

#ifndef __SRC_EM_SPARSE_PHI_SRC_COMPUTE_SPARSE_PHI_HPP_INCLUDED
#define __SRC_EM_SPARSE_PHI_SRC_COMPUTE_SPARSE_PHI_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <cinttypes>
#include <algorithm>
#include <parallel/algorithm>
#include <omp.h>

#include "io/async_stream_reader.hpp"
#include "io/async_stream_writer.hpp"
#include "io/async_multi_stream_writer.hpp"
#include "utils.hpp"


namespace em_sparse_phi_private {

template<typename char_type>
class text_accessor {
  public:
    text_accessor(std::string filename,
        std::uint64_t buf_size = (1UL << 20)) {
      m_buf_size = std::max(1UL, buf_size / sizeof(char_type));
      m_file_size = utils::file_size(filename) / sizeof(char_type);
      m_file = utils::file_open_nobuf(filename, "r");
      m_buf = utils::allocate_array<char_type>(m_buf_size);
      m_buf_pos = 0;
      m_buf_filled = 0;
      m_bytes_read = 0;
    }

    inline char_type access(std::uint64_t i) {
      if (!(m_buf_pos <= i && i < m_buf_pos + m_buf_filled)) {
        if (m_buf_pos + m_buf_filled != i)
          std::fseek(m_file, i * sizeof(char_type), SEEK_SET);
        m_buf_pos = i;
        m_buf_filled = std::min(m_buf_size, m_file_size - m_buf_pos);
        utils::read_from_file(m_buf, m_buf_filled, m_file);
        m_bytes_read += m_buf_filled * sizeof(char_type);
      }
      return m_buf[i - m_buf_pos];
    }

    inline std::uint64_t bytes_read() const {
      return m_bytes_read;
    }

    ~text_accessor() {
      std::fclose(m_file);
      utils::deallocate(m_buf);
    }

  private:
    std::uint64_t m_bytes_read;
    std::uint64_t m_file_size;
    std::uint64_t m_buf_size;
    std::uint64_t m_buf_pos;
    std::uint64_t m_buf_filled;
    char_type *m_buf;
    std::FILE *m_file;
};

template<typename T>
struct key_val_pair {
  key_val_pair() {}
  key_val_pair(T key, T val)
    : m_key(key), m_val(val) {}

  inline bool operator < (const key_val_pair<T> &x) const {
    return m_key < x.m_key;
  }

  T m_key;
  T m_val;
};

//============================================================================
// Note: the only assumption about text_offset_type we have is that it is
// capable of storing values in the range [0..text_length). text_length
// cannot be assigned to variable of type text_offset_type!
//============================================================================

// Returns the position, for which Phi is undefined
template<typename text_offset_type>
std::uint64_t compute_sparse_phi(
    std::string sa_filename,
    std::uint64_t text_length,
    std::uint64_t ram_use,
    std::uint64_t max_halfsegment_size,
    std::uint64_t **halfsegment_pair_count,
    text_offset_type *sparse_phi,
    std::uint64_t plcp_sampling_rate,
    std::uint64_t &total_io_volume) {

  // Print initial message.
  fprintf(stderr, "  Compute sparse Phi: ");

  // Start the timer.
  long double start = utils::wclock();

  // Initialize I/O volume.
  std::uint64_t io_volume = 0;

  // Compute basic parameters.
  std::uint64_t phi_undefined_position = 0;
  std::uint64_t n_halfsegments =
    (text_length + max_halfsegment_size - 1) / max_halfsegment_size;
  std::uint64_t sparse_phi_size =
    (text_length + plcp_sampling_rate - 1) / plcp_sampling_rate;
  std::uint64_t sparse_phi_ram = sparse_phi_size * sizeof(text_offset_type);

#ifdef _OPENMP

  // Allocate space for auxiliary data structures (buffers, etc.).
  static const std::uint64_t opt_in_buf_ram = (32UL << 20);
  static const std::uint64_t opt_par_buf_ram = (4UL << 20);

  std::uint64_t in_buf_ram = opt_in_buf_ram;
  std::uint64_t par_buf_ram = opt_par_buf_ram;

  {
    std::uint64_t total_buf_ram = in_buf_ram + par_buf_ram;
    if (sparse_phi_ram + total_buf_ram > ram_use) {
      std::uint64_t ram_budget = ram_use - sparse_phi_ram;
      long double shrink_factor =
        (long double)ram_budget / (long double)total_buf_ram;
      in_buf_ram = (std::uint64_t)((long double)in_buf_ram * shrink_factor);
      par_buf_ram = (std::uint64_t)((long double)par_buf_ram * shrink_factor);
    }
  }

  typedef async_stream_reader<text_offset_type> sa_reader_type;
  sa_reader_type *sa_reader = new sa_reader_type(sa_filename,
      in_buf_ram, std::max(4UL, in_buf_ram / (2UL << 20)));

  std::uint64_t par_buf_size = par_buf_ram / sizeof(text_offset_type);
  text_offset_type *sa_buffer =
    utils::allocate_array<text_offset_type>(par_buf_size);

  std::uint64_t sa_items_read = 0;
  std::uint64_t prev_sa = 0;
  std::uint64_t step_count = 0;

  while (sa_items_read < text_length) {

    // Print progress message.
    ++step_count;
    if (step_count == (1UL << 28) / par_buf_size) {
      step_count = 0;
      std::uint64_t io_vol = sa_items_read * sizeof(text_offset_type);
      long double elapsed = utils::wclock() - start;
      fprintf(stderr, "\r  Compute sparse Phi: %.1Lf%%, "
          "time = %.2Lfs, I/O = %.2LfMiB/s",
          (100.L * sa_items_read) / text_length,
          elapsed, (1.L * io_vol / (1L << 20)) / elapsed);
    }

    // Read the buffer of SA values.
    std::uint64_t sa_buffer_filled =
      std::min(text_length - sa_items_read, par_buf_size);
    sa_reader->read(sa_buffer, sa_buffer_filled);

    // Process first element separately.
    if (sa_items_read == 0) {
      phi_undefined_position = sa_buffer[0];
    } else {
      std::uint64_t sa_val = sa_buffer[0];
      if ((sa_val % plcp_sampling_rate) == 0)
        sparse_phi[sa_val / plcp_sampling_rate] = prev_sa;
      if (halfsegment_pair_count != NULL) {
        std::uint64_t cur_sa_halfseg_idx = sa_val / max_halfsegment_size;
        std::uint64_t prev_sa_halfseg_idx = prev_sa / max_halfsegment_size;
        ++halfsegment_pair_count[cur_sa_halfseg_idx][prev_sa_halfseg_idx];
      }
    }

    // Process remaining elements in the buffer.
    std::uint64_t max_threads = omp_get_max_threads();
    std::uint64_t max_range_size =
      (sa_buffer_filled + max_threads - 1) / max_threads;
    std::uint64_t n_ranges =
      (sa_buffer_filled + max_range_size - 1) / max_range_size;

    #pragma omp parallel num_threads(n_ranges)
    {
      std::uint64_t thread_id = omp_get_thread_num();
      std::uint64_t range_beg = thread_id * max_range_size;
      std::uint64_t range_end =
        std::min(range_beg + max_range_size, sa_buffer_filled);
      std::uint64_t **local_halfsegment_pair_count = NULL;
      if (halfsegment_pair_count != NULL) {
        local_halfsegment_pair_count = new std::uint64_t*[n_halfsegments];
        for (std::uint64_t j = 0; j < n_halfsegments; ++j) {
          local_halfsegment_pair_count[j] = new std::uint64_t[n_halfsegments];
          std::fill(local_halfsegment_pair_count[j],
              local_halfsegment_pair_count[j] + n_halfsegments, 0UL);
        }
      }

      std::uint64_t prev_segment_id =
        ((std::uint64_t)sa_buffer[std::max(1UL, range_beg) - 1])
        / max_halfsegment_size;
      for (std::uint64_t i = std::max(1UL, range_beg); i < range_end; ++i) {
        std::uint64_t sa_val = sa_buffer[i];
        if ((sa_val % plcp_sampling_rate) == 0)
          sparse_phi[sa_val / plcp_sampling_rate] = sa_buffer[i - 1];

        if (halfsegment_pair_count != NULL) {
          std::uint64_t segment_id = sa_val / max_halfsegment_size;
          ++local_halfsegment_pair_count[segment_id][prev_segment_id];
          prev_segment_id = segment_id;
        }
      }

      if (halfsegment_pair_count != NULL) {
        #pragma omp critical
        {
          for (std::uint64_t i = 0; i < n_halfsegments; ++i)
            for (std::uint64_t j = 0; j < n_halfsegments; ++j)
              halfsegment_pair_count[i][j] +=
                local_halfsegment_pair_count[i][j];
        }
      }

      if (halfsegment_pair_count != NULL) {
        for (std::uint64_t j = 0; j < n_halfsegments; ++j)
          delete[] local_halfsegment_pair_count[j];
        delete[] local_halfsegment_pair_count;
      }
    }

    // Update the number of read SA items.
    prev_sa = sa_buffer[sa_buffer_filled - 1];
    sa_items_read += sa_buffer_filled;
  }

  // Clean up.
  utils::deallocate(sa_buffer);

#else

  // Allocate space for auxiliary data structures (buffers, etc.).
  static const std::uint64_t opt_in_buf_ram = (32UL << 20);
  std::uint64_t in_buf_ram = opt_in_buf_ram;
  if (sparse_phi_ram + in_buf_ram > ram_use)
    in_buf_ram = ram_use - sparse_phi_ram;

  typedef async_stream_reader<text_offset_type> sa_reader_type;
  sa_reader_type *sa_reader = new sa_reader_type(
      sa_filename, in_buf_ram, std::max(4UL, in_buf_ram / (2UL << 20)));

  text_offset_type prev_sa = 0;
  std::uint64_t prev_sa_halfsegment_id = 0;
  std::uint64_t step_count = 0;

  for (std::uint64_t i = 0; i < text_length; ++i) {

    // Print progress message.
    ++step_count;
    if (step_count == (8UL << 20)) {
      step_count = 0;
      std::uint64_t io_vol = i * sizeof(text_offset_type);
      long double elapsed = utils::wclock() - start;
      fprintf(stderr, "\r  Compute sparse Phi: %.1Lf%%, "
          "time = %.2Lfs, I/O = %.2LfMiB/s",
          (100.L * i) / text_length, elapsed,
          (1.L * io_vol / (1L << 20)) / elapsed);
    }

    text_offset_type sa_i = sa_reader->read();
    std::uint64_t sa_i_uint64 = sa_i;

    if (i == 0) {
      phi_undefined_position = sa_i_uint64;
      prev_sa = sa_i;
      prev_sa_halfsegment_id = sa_i_uint64 / max_halfsegment_size;
      continue;
    }

    if ((sa_i_uint64 % plcp_sampling_rate) == 0)
      sparse_phi[sa_i_uint64 / plcp_sampling_rate] = prev_sa;

    if (halfsegment_pair_count != NULL) {
      std::uint64_t halfsegment_id = sa_i_uint64 / max_halfsegment_size;
      ++halfsegment_pair_count[halfsegment_id][prev_sa_halfsegment_id];
      prev_sa_halfsegment_id = halfsegment_id;
    }

    prev_sa = sa_i;
  }
#endif

  // halfsegment_pair_count == NULL means than we
  // do not employ the text-partitioning.
  if (halfsegment_pair_count != NULL) {
    for (std::uint64_t i = 0; i < n_halfsegments; ++i) {
      for (std::uint64_t j = 0; j < n_halfsegments; ++j) {
        if (i > j) {
          halfsegment_pair_count[j][i] += halfsegment_pair_count[i][j];
          halfsegment_pair_count[i][j] = 0;
        }
      }
    }
  }

  // Stop I/O threads.
  sa_reader->stop_reading();

  // Update I/O volume.
  io_volume +=
    sa_reader->bytes_read();
  total_io_volume += io_volume;

  // Print summary.
  long double total_time = utils::wclock() - start;
  fprintf(stderr, "\r  Compute sparse Phi: time = %.2Lfs, "
      "I/O = %.2LfMiB/s, total I/O vol = %.2Lfbytes/input symbol\n",
      total_time, (1.L * io_volume / (1L << 20)) / total_time,
      (1.L * total_io_volume) / text_length);

  // Clean up.
  delete sa_reader;

  // Return the result.
  return phi_undefined_position;
}

// Returns the sparse PLCP array.
template<typename char_type,
  typename text_offset_type>
text_offset_type* compute_sparse_plcp(
    std::string text_filename,
    std::string sa_filename,
    std::string output_filename,
    std::uint64_t text_length,
    std::uint64_t max_halfsegment_size,
    std::uint64_t **halfsegment_pair_count,
    std::uint64_t plcp_sampling_rate,
    std::uint64_t ram_use,
    std::uint64_t &total_io_volume) {

  // Computer basic parameters.
  std::uint64_t sparse_plcp_size =
    (text_length + plcp_sampling_rate - 1) / plcp_sampling_rate;
  std::uint64_t sparse_phi_size = sparse_plcp_size;
  std::uint64_t sparse_plcp_ram = sparse_plcp_size * sizeof(text_offset_type);
  std::uint64_t sparse_phi_ram = sparse_plcp_ram;

  // Print initial message.
  fprintf(stderr, "Compute sparse PLCP:\n");

  // Start the timer.
  long double start = utils::wclock();

#ifdef _OPENMP

  // Allocate mini sparse PLCP array.
  std::uint64_t aux_sparse_plcp_sampling_rate = plcp_sampling_rate *
    std::max(10000UL, ((64UL << 20) / plcp_sampling_rate));
  std::uint64_t aux_sparse_plcp_size =
    (text_length + aux_sparse_plcp_sampling_rate - 1) /
    aux_sparse_plcp_sampling_rate;
  std::uint64_t aux_sparse_plcp_ram =
    aux_sparse_plcp_size * sizeof(std::uint64_t);
  std::uint64_t *aux_sparse_plcp =
    utils::allocate_array<std::uint64_t>(aux_sparse_plcp_size);
#endif

  // Compute sparse Phi from SA.
  text_offset_type *sparse_phi =
    utils::allocate_array<text_offset_type>(sparse_phi_size);
  std::uint64_t phi_undefined_position =
    compute_sparse_phi(sa_filename, text_length, ram_use,
      max_halfsegment_size, halfsegment_pair_count,
      sparse_phi, plcp_sampling_rate, total_io_volume);

  //===========================================================================

  std::uint64_t max_segment_ram = 0;
  std::uint64_t max_segment_size = 0;
  std::uint64_t overflow_buf_size = 0;
  std::uint64_t n_segments = 0;

#ifdef _OPENMP

  // Compute mini sparse PLCP array.
  {
    fprintf(stderr, "  Compute mini sparse PLCP: ");
    long double mini_plcp_start = utils::wclock();
    std::uint64_t io_vol = 0;

    // Allocate mini sparse Phi array.
    std::uint64_t *aux_sparse_phi =
      utils::allocate_array<std::uint64_t>(aux_sparse_plcp_size);
    for (std::uint64_t j = 0; j < aux_sparse_plcp_size; ++j)
      aux_sparse_phi[j] = sparse_phi[j *
        (aux_sparse_plcp_sampling_rate / plcp_sampling_rate)];

    // Compute buffer sizes for text accessor.
    static const std::uint64_t opt_txt_acc_buf_ram = (1UL << 20);
    std::uint64_t txt_acc_buf_ram = opt_txt_acc_buf_ram;
    {
      std::uint64_t total_txt_acc_buf_ram = 2 * txt_acc_buf_ram;
      if (sparse_phi_ram + aux_sparse_plcp_ram +
          total_txt_acc_buf_ram > ram_use) {
        total_txt_acc_buf_ram = ram_use -
          sparse_phi_ram - aux_sparse_plcp_ram;
        txt_acc_buf_ram = total_txt_acc_buf_ram / 2;
      }
    }

    // Allocate text accessors.
    typedef text_accessor<char_type> text_accessor_type;
    text_accessor_type *a1 =
      new text_accessor_type(text_filename, txt_acc_buf_ram);
    text_accessor_type *a2 =
      new text_accessor_type(text_filename, txt_acc_buf_ram);

    // Compute mini sparse PLCP array.
    {

      // Initialize basic stats.
      std::uint64_t lcp = 0;
      std::uint64_t prev_i = 0;
      std::uint64_t prev_lcp = 0;

      // Compute LCP values.
      for (std::uint64_t j = 0; j < aux_sparse_plcp_size; ++j) {
        std::uint64_t i = j * aux_sparse_plcp_sampling_rate;
        std::uint64_t phi_i = aux_sparse_phi[j];

        if (i == phi_undefined_position) lcp = 0;
        else {
          lcp = (std::uint64_t)std::max(0L,
              (std::int64_t)(prev_i + prev_lcp) - (std::int64_t)i);
          while (i + lcp < text_length && phi_i + lcp < text_length &&
              a1->access(i + lcp) == a2->access(phi_i + lcp)) ++lcp;
        }

        aux_sparse_plcp[j] = lcp;
        prev_i = i;
        prev_lcp = lcp;
      }
    }

    // Update I/O volume.
    io_vol +=
      a1->bytes_read() +
      a2->bytes_read();
    total_io_volume += io_vol;

    // Print summary.
    long double mini_plcp_time = utils::wclock() - mini_plcp_start;
    fprintf(stderr, "%.2Lfs, I/O = %.2LfMiB/s, "
        "total I/O vol = %.2Lfbytes/input symbol\n",
        mini_plcp_time, (1.L * io_vol / (1L << 20)) / mini_plcp_time,
        (1.L * total_io_volume) / text_length);

    // Clean up.
    delete a2;
    delete a1;
    utils::deallocate(aux_sparse_phi);
  }
#endif

  // Compute sparse PLCP array.
  {

    // Set the segment RAM to minimal acceptable value.
    max_segment_ram = (std::uint64_t)((long double)ram_use * 0.9L);

    // Shrink buffers or enlarge segment to fill in the available RAM.
#ifndef _OPENMP
    static const std::uint64_t opt_in_buf_ram = (32UL << 20);
    static const std::uint64_t opt_out_buf_ram = (4UL << 20);
    static const std::uint64_t opt_txt_acc_buf_ram = (1UL << 20);
    static const std::uint64_t opt_overflow_buf_ram = (1UL << 20);

    std::uint64_t in_buf_ram = opt_in_buf_ram;
    std::uint64_t out_buf_ram = opt_out_buf_ram;
    std::uint64_t txt_acc_buf_ram = opt_txt_acc_buf_ram;
    std::uint64_t overflow_buf_ram = opt_overflow_buf_ram;

    {
      std::uint64_t total_buf_ram =
        in_buf_ram + out_buf_ram + txt_acc_buf_ram + overflow_buf_ram;
      if (max_segment_ram + total_buf_ram > ram_use) {
        std::uint64_t ram_budget = ram_use - max_segment_ram;
        long double shrink_factor =
          (long double)ram_budget / (long double)total_buf_ram;
        in_buf_ram =
          (std::uint64_t)((long double)in_buf_ram * shrink_factor);
        out_buf_ram =
          (std::uint64_t)((long double)out_buf_ram * shrink_factor);
        txt_acc_buf_ram =
          (std::uint64_t)((long double)txt_acc_buf_ram * shrink_factor);
        overflow_buf_ram =
          (std::uint64_t)((long double)overflow_buf_ram * shrink_factor);
      } else max_segment_ram = ram_use - total_buf_ram;
    }
#else
    static const std::uint64_t opt_in_buf_ram = (2UL << 20);
    static const std::uint64_t opt_out_buf_ram = (4UL << 20);
    static const std::uint64_t opt_txt_acc_buf_ram = (1UL << 20);
    static const std::uint64_t opt_local_buf_ram = (1UL << 20);
    static const std::uint64_t opt_overflow_buf_ram = (1UL << 20);

    std::uint64_t in_buf_ram = opt_in_buf_ram;
    std::uint64_t out_buf_ram = opt_out_buf_ram;
    std::uint64_t txt_acc_buf_ram = opt_txt_acc_buf_ram;
    std::uint64_t local_buf_ram = opt_local_buf_ram;
    std::uint64_t overflow_buf_ram = opt_overflow_buf_ram;

    {
      std::uint64_t max_threads = omp_get_max_threads();
      std::uint64_t total_buf_ram = out_buf_ram + overflow_buf_ram +
        max_threads * (in_buf_ram + local_buf_ram + txt_acc_buf_ram);
      if (max_segment_ram + aux_sparse_plcp_ram + total_buf_ram > ram_use) {
        std::uint64_t ram_budget = ram_use -
          aux_sparse_plcp_ram - max_segment_ram;
        long double shrink_factor =
          (long double)ram_budget / (long double)total_buf_ram;
        in_buf_ram =
          (std::uint64_t)((long double)in_buf_ram * shrink_factor);
        out_buf_ram =
          (std::uint64_t)((long double)out_buf_ram * shrink_factor);
        txt_acc_buf_ram =
          (std::uint64_t)((long double)txt_acc_buf_ram * shrink_factor);
        local_buf_ram =
          (std::uint64_t)((long double)local_buf_ram * shrink_factor);
        overflow_buf_ram =
          (std::uint64_t)((long double)overflow_buf_ram * shrink_factor);
      } else max_segment_ram = ram_use - aux_sparse_plcp_ram - total_buf_ram;
    }
#endif

    max_segment_size = std::max(1UL, max_segment_ram / sizeof(char_type));
    overflow_buf_size = std::max(1UL, overflow_buf_ram / sizeof(char_type));
    n_segments = (text_length + max_segment_size - 1) / max_segment_size;

    fprintf(stderr, "  Segment size = %lu (%.2LfMiB)\n", max_segment_size,
        (1.L * max_segment_size * sizeof(char_type)) / (1UL << 20));
    fprintf(stderr, "  Number of segments = %lu\n", n_segments);

    std::uint64_t max_possible_sparse_phi_for_segment_size =
      (max_segment_size + plcp_sampling_rate - 1) / plcp_sampling_rate;
    std::uint64_t max_possible_sparse_phi_for_segment_ram =
      2UL * sizeof(text_offset_type) *
      max_possible_sparse_phi_for_segment_size;
    std::uint64_t extra_ram_needed = 0;
#ifdef _OPENMP
    extra_ram_needed += aux_sparse_plcp_ram;
#endif

    // Two cases, depending on whether we
    // can sort sparse_phi_for_segment in RAM.
    if (max_possible_sparse_phi_for_segment_ram +
        extra_ram_needed > ram_use) {
      fprintf(stderr, "  Method = distribute by Phi[i]\n");

      // Write every pair (i, Phi[i]) into file
      // corresponding to segment containing Phi[i].
      {

        // Print initial message.
        fprintf(stderr, "  Distribute Phi into buckets: ");

        // Start the timer.
        long double phi_distr_start = utils::wclock();

        // Set the number of free buffers for the writer.
        static const std::uint64_t free_distr_bufs = 4;

        // Compute the buffer size for the writer.
#ifdef _OPENMP
        static const std::uint64_t opt_distr_buf_ram = (4UL << 20);
        std::uint64_t distr_buf_ram = opt_distr_buf_ram;

        {
          std::uint64_t total_distr_bufs_ram =
            (n_segments + free_distr_bufs) * distr_buf_ram;
          if (total_distr_bufs_ram + sparse_phi_ram +
              aux_sparse_plcp_ram > ram_use) {
            total_distr_bufs_ram = ram_use -
              sparse_phi_ram - aux_sparse_plcp_ram;
            distr_buf_ram =
              total_distr_bufs_ram / (n_segments + free_distr_bufs);
          }
        }
#else
        static const std::uint64_t opt_distr_buf_ram = (4UL << 20);
        std::uint64_t distr_buf_ram = opt_distr_buf_ram;

        {
          std::uint64_t total_distr_bufs_ram =
            (n_segments + free_distr_bufs) * distr_buf_ram;
          if (total_distr_bufs_ram + sparse_phi_ram > ram_use) {
            total_distr_bufs_ram = ram_use - sparse_phi_ram;
            distr_buf_ram =
              total_distr_bufs_ram / (n_segments + free_distr_bufs);
          }
        }
#endif

        // Write data to disk.
        typedef async_multi_stream_writer<text_offset_type>
          phi_pair_multi_stream_writer_type;
        phi_pair_multi_stream_writer_type *phi_pair_multi_stream_writer =
          new phi_pair_multi_stream_writer_type(n_segments,
              distr_buf_ram, free_distr_bufs);
        for (std::uint64_t segment_id = 0;
            segment_id < n_segments; ++segment_id)
          phi_pair_multi_stream_writer->add_file(output_filename +
              ".phi_pairs." + utils::intToStr(segment_id));

        // Process all pairs.
        {
          std::uint64_t step_counter = 0;
          for (std::uint64_t i = 0, j = 0; j < sparse_phi_size;
              ++j, i += plcp_sampling_rate) {

            // Print progress message.
            ++step_counter;
            if (step_counter == (16UL << 20)) {
              step_counter = 0;
              long double elapsed = utils::wclock() - phi_distr_start;
              std::uint64_t io_volume = 2UL * j * sizeof(text_offset_type);
              fprintf(stderr, "\r  Distribute Phi into buckets: %.1Lf%%, "
                  "time = %.1Lfs, I/O = %.2LfMiB/s",
                  (100.L * j) / sparse_phi_size, elapsed,
                  (1.L * io_volume / (1L << 20)) / elapsed);
            }

            // Handle special case.
            if (i == phi_undefined_position)
              continue;

            // Compute segment ID containing Phi[i] and write
            // the pair (i, Phi[i]) to the corresponding file.
            text_offset_type phi_i = sparse_phi[j];
            std::uint64_t segment_id = phi_i / max_segment_size;
            phi_pair_multi_stream_writer->write_to_ith_file(segment_id, i);
            phi_pair_multi_stream_writer->write_to_ith_file(segment_id, phi_i);
          }
        }

        // Update I/O volume.
        std::uint64_t io_vol = phi_pair_multi_stream_writer->bytes_written();
        total_io_volume += io_vol;

        // Print summary.
        long double phi_distr_time = utils::wclock() - phi_distr_start;
        fprintf(stderr, "\r  Distribute Phi into buckets: time = %.1Lfs, "
            "I/O = %.2LfMiB/s, total I/O vol = %.2Lfbytes/input symbol\n",
            phi_distr_time, ((1.L * io_vol) / (1L << 20)) / phi_distr_time,
            (1.L * total_io_volume) / text_length);

        // Clean up.
        delete phi_pair_multi_stream_writer;
      }

      // Clean up.
      utils::deallocate(sparse_phi);

      // Load every segment and compute PLCP for pairs
      // (i, Phi[i]), where Phi is inside the segment.
      {
        for (std::uint64_t segment_id = 0;
            segment_id < n_segments; ++segment_id) {
          std::uint64_t segment_beg = segment_id * max_segment_size;
          std::uint64_t segment_end =
            std::min(segment_beg + max_segment_size, text_length);
          std::uint64_t ext_segment_end =
            std::min(segment_end + overflow_buf_size, text_length);
          std::uint64_t ext_segment_size = ext_segment_end - segment_beg;

          fprintf(stderr, "  Process segment %lu/%lu [%lu..%lu):\n",
              segment_id + 1, n_segments, segment_beg, segment_end);

          std::string phi_pairs_filename = output_filename +
            ".phi_pairs." + utils::intToStr(segment_id);
          std::string lcp_pairs_filename = output_filename +
            ".lcp_pairs." + utils::intToStr(segment_id);

          typedef async_stream_reader<text_offset_type> phi_pair_reader_type;
          typedef async_stream_writer<text_offset_type> lcp_pair_writer_type;
          typedef async_stream_reader<char_type> text_reader_type;
          typedef text_accessor<char_type> text_accessor_type;

          // Allocate the segment.
          char_type *segment =
            utils::allocate_array<char_type>(ext_segment_size);

          // Read the segment from disk.
          {
            fprintf(stderr, "    Read segment: ");
            std::uint64_t read_io_vol = 0;
            long double read_start = utils::wclock();
            std::uint64_t offset = segment_beg * sizeof(char_type);
            utils::read_at_offset(segment, offset,
                ext_segment_size, text_filename);
            read_io_vol = ext_segment_size * sizeof(char_type);
            total_io_volume += read_io_vol;
            long double read_time = utils::wclock() - read_start;
            fprintf(stderr, "%.2Lfs, I/O = %.2LfMiB/s, total "
                "I/O vol = %.2Lfbytes/input symbol\n", read_time,
                ((1.L * read_io_vol) / (1L << 20)) / read_time,
                (1.L * total_io_volume) / text_length);
          }

          // Process the segment
          fprintf(stderr, "    Process (i, Phi[i]) pairs: ");
          long double process_start = utils::wclock();
          std::uint64_t io_vol = 0;

          // Initialize readers and writers.
#ifndef _OPENMP
          lcp_pair_writer_type *lcp_pair_writer =
            new lcp_pair_writer_type(lcp_pairs_filename,
                out_buf_ram, std::max(4UL, out_buf_ram / (2UL << 20)));
          text_accessor_type *accessor =
            new text_accessor_type(text_filename, txt_acc_buf_ram);
          phi_pair_reader_type *phi_pair_reader =
            new phi_pair_reader_type(phi_pairs_filename, in_buf_ram / 2,
                std::max(4UL, in_buf_ram / (2UL << 20)));
          text_reader_type *text_reader =
            new text_reader_type(text_filename, in_buf_ram / 2,
                std::max(4UL, in_buf_ram / (2UL << 20)));

          if (phi_pair_reader->empty() == false) {
            std::uint64_t i = (std::uint64_t)phi_pair_reader->read();
            std::uint64_t phi_i = (std::uint64_t)phi_pair_reader->read();
            std::uint64_t lcp = 0;

            std::uint64_t buf_beg = 0;
            while (buf_beg < text_length) {
              text_reader->receive_new_buffer();
              std::uint64_t buf_filled = text_reader->get_buf_filled();
              const char_type *buffer = text_reader->get_buf_ptr();

              while (true) {
                while (phi_i + lcp < text_length &&
                    i + lcp < buf_beg + buf_filled) {
                  char_type next_char = (phi_i + lcp < ext_segment_end) ?
                    segment[(phi_i + lcp) - segment_beg] :
                    accessor->access(phi_i + lcp);
                  if (next_char == buffer[(i + lcp) - buf_beg]) ++lcp;
                  else break;
                }

                if (i + lcp < buf_beg + buf_filled ||
                    phi_i + lcp == text_length ||
                    buf_beg + buf_filled == text_length) {

                  lcp_pair_writer->write((text_offset_type)i);
                  lcp_pair_writer->write((text_offset_type)lcp);
                  if (phi_pair_reader->empty() == false) {
                    std::uint64_t next_i =
                      (std::uint64_t)phi_pair_reader->read();
                    std::uint64_t next_phi_i =
                      (std::uint64_t)phi_pair_reader->read();
                    lcp = std::max(0L,
                        (std::int64_t)(i + lcp) - (std::int64_t)next_i);
                    i = next_i;
                    phi_i = next_phi_i;
                  } else {
                    i = text_length;
                    break;
                  }
                } else break;
              }

              if (i == text_length)
                break;

              buf_beg += buf_filled;
            }
          }

          // Stop I/O threads.
          text_reader->stop_reading();
          phi_pair_reader->stop_reading();

          // Update I/O volume.
          io_vol +=
            text_reader->bytes_read() +
            accessor->bytes_read() +
            lcp_pair_writer->bytes_written() +
            phi_pair_reader->bytes_read();
          total_io_volume += io_vol;

          // Print summary.
          long double process_time = utils::wclock() - process_start;
          fprintf(stderr, "\r    Process (i, Phi[i]) pairs: time = %.2Lfs, "
              "I/O = %.2LfMiB/s, total I/O vol = %.2Lfbytes/input symbol\n",
              process_time, (1.L * io_vol / (1L << 20)) / process_time,
              (1.L * total_io_volume) / text_length);

          // Clean up.
          delete text_reader;
          delete phi_pair_reader;
          delete accessor;
          delete lcp_pair_writer;

#else

          // Parallel computation of sparse PLCP values assuming each
          // thread handles the same amount of (i, Phi[i]) pairs (but
          // different threads can stream different amount of text).
          lcp_pair_writer_type *lcp_pair_writer =
            new lcp_pair_writer_type(lcp_pairs_filename, out_buf_ram, 2);

          std::uint64_t n_pairs =
            utils::file_size(phi_pairs_filename) /
            (2 * sizeof(text_offset_type));

          std::uint64_t max_threads = omp_get_max_threads();
          std::uint64_t max_range_size =
            (n_pairs + max_threads - 1) / max_threads;
          std::uint64_t n_ranges =
            (n_pairs + max_range_size - 1) / max_range_size;

          #pragma omp parallel num_threads(n_ranges)
          {
            std::uint64_t thread_id = omp_get_thread_num();
            std::uint64_t range_beg = thread_id * max_range_size;
            std::uint64_t range_end =
              std::min(n_pairs, range_beg + max_range_size);
            std::uint64_t range_size = range_end - range_beg;

            text_accessor_type *accessor =
              new text_accessor_type(text_filename, txt_acc_buf_ram);
            phi_pair_reader_type *phi_pair_reader =
              new phi_pair_reader_type(phi_pairs_filename,
                  in_buf_ram / 2, 2, 2UL * range_beg);

            std::uint64_t first_i = (std::uint64_t)phi_pair_reader->peek();
            std::uint64_t sample_plcp_addr =
              first_i / aux_sparse_plcp_sampling_rate;
            std::uint64_t sample_plcp_pos =
              sample_plcp_addr * aux_sparse_plcp_sampling_rate;
            std::uint64_t sample_dist = first_i - sample_plcp_pos;
            std::uint64_t sample_plcp_val =
              aux_sparse_plcp[sample_plcp_addr];
            std::uint64_t first_i_plcp_lower_bound =
              (std::uint64_t)std::max(0L,
                (std::int64_t)sample_plcp_val -
                (std::int64_t)sample_dist);
            std::uint64_t text_pos = first_i + first_i_plcp_lower_bound;

            text_reader_type *text_reader =
              new text_reader_type(text_filename, in_buf_ram / 2, 2, text_pos);
            std::uint64_t local_buf_size =
              local_buf_ram / (2 * sizeof(text_offset_type));
            text_offset_type *local_buf =
              utils::allocate_array<text_offset_type>(2 * local_buf_size);

            std::uint64_t lcp = 0;
            std::uint64_t pairs_processed = 0;
            while (pairs_processed < range_size) {
              std::uint64_t local_buf_filled =
                std::min(local_buf_size, range_size - pairs_processed);
              for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
                std::uint64_t i = (std::uint64_t)phi_pair_reader->read();
                std::uint64_t phi_i = (std::uint64_t)phi_pair_reader->read();
                lcp = (std::uint64_t)std::max(0L,
                    (std::int64_t)text_pos - (std::int64_t)i);
                if (i + lcp != text_pos) {
                  text_reader->skip((i + lcp) - text_pos);
                  text_pos = i + lcp;
                }

                while (text_pos < text_length &&
                    phi_i + lcp < text_length) {
                  char_type next_char = (phi_i + lcp < ext_segment_end) ?
                    segment[(phi_i + lcp) - segment_beg] :
                    accessor->access(phi_i + lcp);
                  if (next_char == text_reader->peek()) {
                    ++lcp;
                    ++text_pos;
                    text_reader->read();
                  } else break;
                }

                local_buf[2 * j] = (text_offset_type)i;
                local_buf[2 * j + 1] = (text_offset_type)lcp;
              }

              #pragma omp critical
              {
                lcp_pair_writer->write(local_buf, 2 * local_buf_filled);
              }
              pairs_processed += local_buf_filled;
            }

            #pragma omp critical
            {
              // Stop I/O threads.
              text_reader->stop_reading();
              phi_pair_reader->stop_reading();

              // Update I/O volume.
              io_vol +=
                accessor->bytes_read() +
                text_reader->bytes_read() +
                phi_pair_reader->bytes_read();
            }

            // Clean up.
            utils::deallocate(local_buf);
            delete text_reader;
            delete phi_pair_reader;
            delete accessor;
          }

          // Update I/O volume.
          io_vol +=
            lcp_pair_writer->bytes_written();
          total_io_volume += io_vol;

          // Print summary.
          long double process_time = utils::wclock() - process_start;
          fprintf(stderr, "\r    Process (i, Phi[i]) pairs: time = %.2Lfs, "
              "I/O = %.2LfMiB/s, total I/O vol = %.2Lfbytes/input symbol\n",
              process_time, (1.L * io_vol / (1L << 20)) / process_time,
              (1.L * total_io_volume) / text_length);

          // Clean up.
          delete lcp_pair_writer;
#endif

          // Clean up.
          utils::deallocate(segment);
          utils::file_delete(phi_pairs_filename);
        }
      }
    } else {
      fprintf(stderr, "  Method = distribute by max(i, Phi[i])\n");

      // Write every pair (i, Phi[i]) into file
      // corresponding to segment containing max(i, Phi[i]).
#ifdef _OPENMP

      // Bucket size used to partition (i, Phi[i]) pairs, so that
      // each thread processes roughly equal number of pairs.
      // Note: this array is assumed to be small, so we don't take it into
      // consideration when computing exact RAM usages. E.g., for 256GiB text
      // and 3.5GiB of RAM it takes only about 9MiB (XXX).
      static const std::uint64_t bucket_size_log = 24UL;
      static const std::uint64_t bucket_size = (1UL << bucket_size_log);

      // XXX utils::allocate array?
      // XXX within ram use?
      std::vector<std::uint64_t> **type_1_bucket_sizes =
        new std::vector<std::uint64_t>*[n_segments];
      for (std::uint64_t segment_id = 0;
          segment_id < n_segments; ++segment_id) {
        std::uint64_t segment_beg = segment_id * max_segment_size;
        std::uint64_t segment_end =
          std::min(text_length, segment_beg + max_segment_size);
        std::uint64_t n_buckets =
          (segment_end + bucket_size - 1) / bucket_size;
        type_1_bucket_sizes[segment_id] =
          new std::vector<std::uint64_t>(n_buckets, 0UL);
      }
#endif

      typedef key_val_pair<text_offset_type> pair_type;

      {
        fprintf(stderr, "  Distribute Phi into buckets: ");
        long double phi_distr_start = utils::wclock();

        // Initialize distr_buf_ram.
        static const std::uint64_t free_distr_bufs = 4;
        static const std::uint64_t opt_distr_buf_ram = (4UL << 20);
        std::uint64_t distr_buf_ram = opt_distr_buf_ram;

        // Shrink distr_buf_ram if necessary to fit in RAM budget.
        {
          std::uint64_t total_distr_bufs_ram =
            (2 * n_segments + free_distr_bufs) * distr_buf_ram;
#ifdef _OPEMP
          if (total_distr_bufs_ram + sparse_phi_ram +
              aux_sparse_plcp_ram > ram_use) {
            total_distr_bufs_ram = ram_use -
              sparse_phi_ram - aux_sparse_plcp_ram;
            distr_buf_ram =
              total_distr_bufs_ram / (2 * n_segments + free_distr_bufs);
          }
#else
          if (total_distr_bufs_ram + sparse_phi_ram > ram_use) {
            total_distr_bufs_ram = ram_use - sparse_phi_ram;
            distr_buf_ram =
              total_distr_bufs_ram / (2 * n_segments + free_distr_bufs);
          }
#endif
        }

        // Create the multi writer.
        typedef async_multi_stream_writer<pair_type>
          phi_pair_multi_stream_writer_type;
        phi_pair_multi_stream_writer_type *phi_pair_multi_stream_writer =
          new phi_pair_multi_stream_writer_type(2UL * n_segments,
              distr_buf_ram, free_distr_bufs);

        // Add the files to multiwriter.
        for (std::uint64_t segment_id = 0;
            segment_id < n_segments; ++segment_id) {
          phi_pair_multi_stream_writer->add_file(
              output_filename + ".phi_pairs." +
              utils::intToStr(segment_id) + ".1");
          phi_pair_multi_stream_writer->add_file(
              output_filename + ".phi_pairs." +
              utils::intToStr(segment_id) + ".2");
        }

        // Process all pairs.
        {
          std::uint64_t step_counter = 0;
          for (std::uint64_t i = 0, j = 0; j < sparse_phi_size;
              ++j, i += plcp_sampling_rate) {

            // Print the progress message.
            ++step_counter;
            if (step_counter == (16UL << 20)) {
              step_counter = 0;
              long double elapsed = utils::wclock() - phi_distr_start;
              std::uint64_t io_volume =
                phi_pair_multi_stream_writer->bytes_written();
              fprintf(stderr, "\r  Distribute Phi into buckets: %.1Lf%%, "
                  "time = %.1Lfs, I/O = %.2LfMiB/s",
                  (100.L * j) / sparse_phi_size, elapsed,
                  (1.L * io_volume / (1L << 20)) / elapsed);
            }

            // Handle special case.
            if (i == phi_undefined_position)
              continue;

            // Write data to disk.
            text_offset_type phi_i = sparse_phi[j];
            if (i < (std::uint64_t)phi_i) {
              std::uint64_t segment_id =
                (std::uint64_t)phi_i / max_segment_size;
              phi_pair_multi_stream_writer->write_to_ith_file(
                  segment_id * 2, pair_type(i, phi_i));

#ifdef _OPENMP
              std::uint64_t bucket_id = (i >> bucket_size_log);
              (*type_1_bucket_sizes[segment_id])[bucket_id] += 1;
#endif

            } else {
              std::uint64_t segment_id = i / max_segment_size;
              phi_pair_multi_stream_writer->write_to_ith_file(
                  segment_id * 2 + 1, pair_type(phi_i, i));
            }
          }
        }

        // Update I/O volume.
        std::uint64_t io_vol = phi_pair_multi_stream_writer->bytes_written();
        total_io_volume += io_vol;

        // Print summary.
        long double phi_distr_time = utils::wclock() - phi_distr_start;
        fprintf(stderr, "\r  Distribute Phi into buckets: time = %.1Lfs, "
            "I/O = %.2LfMiB/s, total I/O vol = %.2Lfbytes/input symbol\n",
            phi_distr_time, ((1.L * io_vol) / (1L << 20)) / phi_distr_time,
            (1.L * total_io_volume) / text_length);

        // Clean up.
        delete phi_pair_multi_stream_writer;
      }

      // Clean up.
      utils::deallocate(sparse_phi);

      // Load every segment and compute PLCP for pairs
      // (i, Phi[i]), where Phi is inside the segment.
      {
        for (std::uint64_t segment_id = 0;
            segment_id < n_segments; ++segment_id) {
          std::uint64_t segment_beg = segment_id * max_segment_size;
          std::uint64_t segment_end =
            std::min(segment_beg + max_segment_size, text_length);
          std::uint64_t ext_segment_end =
            std::min(segment_end + overflow_buf_size, text_length);
          std::uint64_t ext_segment_size = ext_segment_end - segment_beg;

          fprintf(stderr, "  Process segment %lu/%lu [%lu..%lu):\n",
              segment_id + 1, n_segments, segment_beg, segment_end);

          std::string phi_pairs_filename_1 =
            output_filename + ".phi_pairs." +
            utils::intToStr(segment_id) + ".1";
          std::string phi_pairs_filename_2 =
            output_filename + ".phi_pairs." +
            utils::intToStr(segment_id) + ".2";
          std::string lcp_pairs_filename =
            output_filename + ".lcp_pairs." +
            utils::intToStr(segment_id);

          std::uint64_t n_pairs_type_1 =
            utils::file_size(phi_pairs_filename_1) / sizeof(pair_type);
          std::uint64_t n_pairs_type_2 =
            utils::file_size(phi_pairs_filename_2) / sizeof(pair_type);

          typedef async_stream_reader<pair_type> phi_pair_reader_type;
          typedef async_stream_reader<char_type> text_reader_type;
          typedef async_stream_writer<text_offset_type> lcp_pair_writer_type;
          typedef text_accessor<char_type> text_accessor_type;

#ifdef _OPENMP
          std::uint64_t n_buckets =
            (segment_end + bucket_size - 1) / bucket_size;
          std::vector<std::uint64_t> type_2_bucket_sizes(n_buckets, 0UL);
#endif

          // Read the file containing pairs (i, Phi[i]) such that Phi[i] < i,
          // sort by Phi[i] and write to disk (overwriting the old file).
          {
            fprintf(stderr, "    Preprocess (i, Phi[i]) pairs: ");
            long double preprocess_start = utils::wclock();
            std::uint64_t sparse_phi_for_segment_size = n_pairs_type_2;
            pair_type *sparse_phi_for_segment =
              utils::allocate_array<pair_type>(sparse_phi_for_segment_size);
            utils::read_from_file(sparse_phi_for_segment,
                sparse_phi_for_segment_size, phi_pairs_filename_2);
            total_io_volume +=
              sparse_phi_for_segment_size * sizeof(sparse_phi_for_segment[0]);
            std::sort(sparse_phi_for_segment,
                sparse_phi_for_segment + sparse_phi_for_segment_size);

#ifdef _OPENMP
            for (std::uint64_t i = 0; i < sparse_phi_for_segment_size; ++i) {
              std::uint64_t bucket_id =
                ((std::uint64_t)(sparse_phi_for_segment[i].m_key)) >>
                bucket_size_log;
              type_2_bucket_sizes[bucket_id] += 1;
            }
#endif

            // Write sorted sparse_phi_for_segment back to disk.
            // We have to delete the file first, to make sure the functions
            // counting disk space usage correctly handle this overwrite.
            if (utils::file_exists(phi_pairs_filename_2))
              utils::file_delete(phi_pairs_filename_2);
            utils::write_to_file(sparse_phi_for_segment,
                sparse_phi_for_segment_size, phi_pairs_filename_2);

            // Update I/O volume.
            total_io_volume +=
              sparse_phi_for_segment_size * sizeof(sparse_phi_for_segment[0]);

            // Print summary.
            long double preprocess_time = utils::wclock() - preprocess_start;
            fprintf(stderr, "%.2Lfs, total I/O vol = "
                "%.2Lfbytes/input symbol\n",
                preprocess_time,
                (1.L * total_io_volume) / text_length);

            // Clean up.
            utils::deallocate(sparse_phi_for_segment);
          }

#ifdef _OPENMP
          std::uint64_t max_threads = omp_get_max_threads();
          std::uint64_t max_range_size =
            (segment_end + max_threads - 1) / max_threads;
          std::uint64_t n_ranges =
            (segment_end + max_range_size - 1) / max_range_size;
          std::vector<std::uint64_t> range_begs(n_ranges);
          std::vector<std::uint64_t> type_1_ptrs(n_ranges);
          std::vector<std::uint64_t> type_2_ptrs(n_ranges);

          std::uint64_t n_pairs = n_pairs_type_1 + n_pairs_type_2;
          std::uint64_t ideal_range_size = n_pairs / n_ranges;

          // Make the size of each bucket as close as possible to
          // ideal_range_size. Alternative to this would be to find
          // the smallest upper bound on the size of all buckets.
          // It could be binary searched (to check a single value
          // during binary search, we greedily partition the buckets
          // and then check if the number of groups is <= n_ranges).
          std::uint64_t bucket_range_beg = 0;
          std::uint64_t type_1_buckets_total_size = 0;
          std::uint64_t type_2_buckets_total_size = 0;
          for (std::uint64_t range_id = 0; range_id < n_ranges; ++range_id) {
            range_begs[range_id] = bucket_range_beg * bucket_size;
            type_1_ptrs[range_id] = type_1_buckets_total_size;
            type_2_ptrs[range_id] = type_2_buckets_total_size;

            std::uint64_t bucket_range_end = bucket_range_beg;
            std::uint64_t cur_range_size = 0;

            // Add buckets to the current range as long as adding a new bucket
            // brings the range size closer to the ideal_range_size.
            while (bucket_range_end < n_buckets &&
                (range_id + 1 == n_ranges ||
                 std::abs(
                   (std::int64_t)ideal_range_size -
                   (std::int64_t)cur_range_size) >=
                 std::abs(
                   (std::int64_t)ideal_range_size -
                   (std::int64_t)(cur_range_size +
                     (*type_1_bucket_sizes[segment_id])[bucket_range_end] +
                     type_2_bucket_sizes[bucket_range_end])))) {

              cur_range_size +=
                (*type_1_bucket_sizes[segment_id])[bucket_range_end] +
                type_2_bucket_sizes[bucket_range_end];
              type_1_buckets_total_size +=
                (*type_1_bucket_sizes[segment_id])[bucket_range_end];
              type_2_buckets_total_size +=
                type_2_bucket_sizes[bucket_range_end];
              ++bucket_range_end;
            }

            bucket_range_beg = bucket_range_end;
          }

#endif

          // Allocate the segment.
          char_type *segment =
            utils::allocate_array<char_type>(ext_segment_size);

          // Read the segment from disk.
          {

            // Print initial message.
            fprintf(stderr, "    Read segment: ");

            // Initialize I/O volume.
            std::uint64_t read_io_vol = 0;

            // Initialize the timer.
            long double read_start = utils::wclock();

            // Read the data from disk.
            std::uint64_t offset = segment_beg * sizeof(char_type);
            utils::read_at_offset(segment,
                offset, ext_segment_size, text_filename);

            // Update I/O volume.
            read_io_vol += ext_segment_size * sizeof(char_type);
            total_io_volume += read_io_vol;

            // Print summary.
            long double read_time = utils::wclock() - read_start;
            fprintf(stderr, "%.2Lfs, I/O = %.2LfMiB/s, "
                "total I/O vol = %.2Lfbytes/input symbol\n",
                read_time, ((1.L * read_io_vol) / (1L << 20)) / read_time,
              (1.L * total_io_volume) / text_length);
          }

          // Process the segment.
          fprintf(stderr, "    Process (i, Phi[i]) pairs: ");
          long double process_start = utils::wclock();
          std::uint64_t io_vol = 0;

#ifndef _OPENMP

          // Initialize readers and writers.
          text_reader_type *text_reader =
            new text_reader_type(text_filename,
                in_buf_ram / 2, std::max(4UL, (in_buf_ram / (4UL << 20))));
          lcp_pair_writer_type *lcp_pair_writer =
            new lcp_pair_writer_type(lcp_pairs_filename,
                out_buf_ram, std::max(4UL, out_buf_ram / (2UL << 20)));
          text_accessor_type *accessor_1 =
            new text_accessor_type(text_filename,
                txt_acc_buf_ram / 2);
          text_accessor_type *accessor_2 =
            new text_accessor_type(text_filename,
                txt_acc_buf_ram / 2);
          phi_pair_reader_type *phi_pair_reader_1 =
            new phi_pair_reader_type(phi_pairs_filename_1,
                in_buf_ram / 4, std::max(4UL, in_buf_ram / (4UL << 20)));
          phi_pair_reader_type *phi_pair_reader_2 =
            new phi_pair_reader_type(phi_pairs_filename_2,
                in_buf_ram / 4, std::max(4UL, in_buf_ram / (4UL << 20)));

          std::uint64_t i_1 = text_length;
          std::uint64_t phi_i_1 = 0;
          std::uint64_t lcp_1 = 0;
          if (n_pairs_type_1 > 0) {
            pair_type pp = phi_pair_reader_1->read();
            i_1 = (std::uint64_t)pp.m_key;
            phi_i_1 = (std::uint64_t)pp.m_val;
          }

          std::uint64_t i_2 = text_length;
          std::uint64_t phi_i_2 = 0;
          std::uint64_t lcp_2 = 0;
          if (n_pairs_type_2 > 0) {
            pair_type pp = phi_pair_reader_2->read();
            phi_i_2 = (std::uint64_t)pp.m_key;
            i_2 = (std::uint64_t)pp.m_val;
          }

          std::uint64_t buf_beg = 0;
          while (buf_beg < text_length &&
              (i_1 != text_length || i_2 != text_length)) {

            text_reader->receive_new_buffer();
            std::uint64_t buf_filled = text_reader->get_buf_filled();
            const char_type *buffer = text_reader->get_buf_ptr();

            while (i_1 != text_length) {
              while (phi_i_1 + lcp_1 < text_length &&
                  i_1 + lcp_1 < buf_beg + buf_filled) {
                char_type next_char = (phi_i_1 + lcp_1 < ext_segment_end) ?
                  segment[(phi_i_1 + lcp_1) - segment_beg] :
                  accessor_1->access(phi_i_1 + lcp_1);
                if (next_char == buffer[(i_1 + lcp_1) - buf_beg])
                  ++lcp_1;
                else break;
              }

              if (i_1 + lcp_1 < buf_beg + buf_filled ||
                  phi_i_1 + lcp_1 == text_length ||
                  buf_beg + buf_filled == text_length) {
                lcp_pair_writer->write((text_offset_type)i_1);
                lcp_pair_writer->write((text_offset_type)lcp_1);

                if (phi_pair_reader_1->empty() == false) {
                  pair_type pp = phi_pair_reader_1->read();
                  std::uint64_t next_i = (std::uint64_t)pp.m_key;
                  std::uint64_t next_phi_i = (std::uint64_t)pp.m_val;
                  lcp_1 = std::max(0L,
                      (std::int64_t)(i_1 + lcp_1) -
                      (std::int64_t)next_i);
                  i_1 = next_i;
                  phi_i_1 = next_phi_i;
                } else i_1 = text_length;
              } else break;
            }

            while (i_2 != text_length) {
              while (i_2 + lcp_2 < text_length &&
                  phi_i_2 + lcp_2 < buf_beg + buf_filled) {
                char_type next_char =
                  (i_2 + lcp_2 < ext_segment_end) ?
                  segment[(i_2 + lcp_2) - segment_beg] :
                  accessor_2->access(i_2 + lcp_2);
                if (next_char == buffer[(phi_i_2 + lcp_2) - buf_beg])
                  ++lcp_2;
                else break;
              }

              if (phi_i_2 + lcp_2 < buf_beg + buf_filled ||
                  i_2 + lcp_2 == text_length ||
                  buf_beg + buf_filled == text_length) {
                lcp_pair_writer->write((text_offset_type)i_2);
                lcp_pair_writer->write((text_offset_type)lcp_2);

                if (phi_pair_reader_2->empty() == false) {
                  pair_type pp = phi_pair_reader_2->read();
                  std::uint64_t next_phi_i = (std::uint64_t)pp.m_key;
                  std::uint64_t next_i = (std::uint64_t)pp.m_val;
                  lcp_2 = std::max(0L,
                      (std::int64_t)(phi_i_2 + lcp_2) -
                      (std::int64_t)next_phi_i);
                  i_2 = next_i;
                  phi_i_2 = next_phi_i;
                } else i_2 = text_length;
              } else break;
            }

            buf_beg += buf_filled;
          }

          // Stop I/O threads.
          text_reader->stop_reading();
          phi_pair_reader_1->stop_reading();
          phi_pair_reader_2->stop_reading();

          // Update I/O volume.
          io_vol +=
            text_reader->bytes_read() +
            accessor_1->bytes_read() +
            accessor_2->bytes_read() +
            lcp_pair_writer->bytes_written() +
            phi_pair_reader_1->bytes_read() +
            phi_pair_reader_2->bytes_read();
          total_io_volume += io_vol;

          // Print summary.
          long double process_time = utils::wclock() - process_start;
          fprintf(stderr, "\r    Process (i, Phi[i]) pairs: time = %.2Lfs, "
              "I/O = %.2LfMiB/s, total I/O vol = %.2Lfbytes/input symbol\n",
              process_time, (1.L * io_vol / (1L << 20)) / process_time,
              (1.L * total_io_volume) / text_length);

          // Clean up.
          delete phi_pair_reader_2;
          delete phi_pair_reader_1;
          delete accessor_2;
          delete accessor_1;
          delete lcp_pair_writer;
          delete text_reader;

#else

          // Initialize the output writer.
          lcp_pair_writer_type *lcp_pair_writer =
            new lcp_pair_writer_type(lcp_pairs_filename,
                out_buf_ram, std::max(4UL, out_buf_ram / (2UL << 20)));

          #pragma omp parallel num_threads(n_ranges)
          {

            // Parallel computation of sparse PLCP values assuming each
            // thread handles the same amount of text (but different
            // threads can handle different amount of (i, Phi[i]) pairs.
            // Each thread processes a range of text, in a sense that
            // processed are all pairs with the beginning position
            // inside that range of text.
            std::uint64_t range_id = omp_get_thread_num();
            std::uint64_t range_beg = range_begs[range_id];
            std::uint64_t range_end = 0;
            if (range_id + 1 == n_ranges) range_end = segment_end;
            else range_end = range_begs[range_id + 1];

            std::uint64_t type_1_ptr = type_1_ptrs[range_id];
            std::uint64_t type_2_ptr = type_2_ptrs[range_id];
            if (type_1_ptr != n_pairs_type_1 ||
                type_2_ptr != n_pairs_type_2) {

              // At least one pair (i, Ph[i]) with i >= range_beg was found.
              std::uint64_t local_buf_size =
                local_buf_ram / (2 * sizeof(text_offset_type));
              text_offset_type *local_buf =
                utils::allocate_array<text_offset_type>(2 * local_buf_size);

              std::uint64_t local_buf_filled = 0;

              text_reader_type *text_reader =
                new text_reader_type(text_filename,
                    in_buf_ram / 2, 2, range_beg);
              text_accessor_type *accessor_1 =
                new text_accessor_type(text_filename,
                    txt_acc_buf_ram / 2);
              text_accessor_type *accessor_2 =
                new text_accessor_type(text_filename,
                    txt_acc_buf_ram / 2);
              phi_pair_reader_type *phi_pair_reader_1 =
                new phi_pair_reader_type(phi_pairs_filename_1,
                    in_buf_ram / 4, 2, type_1_ptr);
              phi_pair_reader_type *phi_pair_reader_2 =
                new phi_pair_reader_type(phi_pairs_filename_2,
                    in_buf_ram / 4, 2, type_2_ptr);

              std::uint64_t i_1 = text_length;
              std::uint64_t phi_i_1 = 0;
              std::uint64_t lcp_1 = 0;
              if (type_1_ptr != n_pairs_type_1 &&
                  (std::uint64_t)(phi_pair_reader_1->peek().m_key) <
                  range_end) {
                pair_type pp = phi_pair_reader_1->read();
                i_1 = (std::uint64_t)pp.m_key;
                phi_i_1 = (std::uint64_t)pp.m_val;

                // Compute lower bound for PLCP[i_1]
                // using mini sparse PLCP array.
                std::uint64_t sample_plcp_addr =
                  i_1 / aux_sparse_plcp_sampling_rate;
                std::uint64_t sample_plcp_pos =
                  sample_plcp_addr * aux_sparse_plcp_sampling_rate;
                std::uint64_t sample_dist = i_1 - sample_plcp_pos;
                std::uint64_t sample_plcp_val =
                  aux_sparse_plcp[sample_plcp_addr];
                std::uint64_t first_i_plcp_lower_bound =
                  (std::uint64_t)std::max(0L,
                      (std::int64_t)sample_plcp_val -
                      (std::int64_t)sample_dist);
                lcp_1 = first_i_plcp_lower_bound;
              }

              std::uint64_t i_2 = text_length;
              std::uint64_t phi_i_2 = 0;
              std::uint64_t lcp_2 = 0;
              if (type_2_ptr != n_pairs_type_2 &&
                  (std::uint64_t)(phi_pair_reader_2->peek().m_key) <
                  range_end) {
                pair_type pp = phi_pair_reader_2->read();
                phi_i_2 = (std::uint64_t)pp.m_key;
                i_2 = (std::uint64_t)pp.m_val;

                // Compute lower bound for PLCP[i_1]
                // using mini sparse PLCP array.
                std::uint64_t sample_plcp_addr =
                  i_2 / aux_sparse_plcp_sampling_rate;
                std::uint64_t sample_plcp_pos =
                  sample_plcp_addr * aux_sparse_plcp_sampling_rate;
                std::uint64_t sample_dist = i_2 - sample_plcp_pos;
                std::uint64_t sample_plcp_val =
                  aux_sparse_plcp[sample_plcp_addr];
                std::uint64_t first_i_plcp_lower_bound =
                  (std::uint64_t)std::max(0L,
                      (std::int64_t)sample_plcp_val -
                      (std::int64_t)sample_dist);
                lcp_2 = first_i_plcp_lower_bound;
              }

              std::uint64_t buf_beg = range_beg;
              while (buf_beg < text_length &&
                  (i_1 != text_length || i_2 != text_length)) {

                text_reader->receive_new_buffer();
                std::uint64_t buf_filled = text_reader->get_buf_filled();
                const char_type *buffer = text_reader->get_buf_ptr();

                while (i_1 != text_length) {
                  while (phi_i_1 + lcp_1 < text_length &&
                      i_1 + lcp_1 < buf_beg + buf_filled) {
                    char_type next_char =
                      (phi_i_1 + lcp_1 < ext_segment_end) ?
                      segment[(phi_i_1 + lcp_1) - segment_beg] :
                      accessor_1->access(phi_i_1 + lcp_1);

                    if (next_char == buffer[(i_1 + lcp_1) - buf_beg])
                      ++lcp_1;
                    else break;
                  }

                  if (i_1 + lcp_1 < buf_beg + buf_filled ||
                      phi_i_1 + lcp_1 == text_length ||
                      buf_beg + buf_filled == text_length) {
                    local_buf[2 * local_buf_filled] = (text_offset_type)i_1;
                    local_buf[2 * local_buf_filled + 1] =
                      (text_offset_type)lcp_1;
                    ++local_buf_filled;

                    if (local_buf_filled == local_buf_size) {
                      #pragma omp critical
                      {
                        lcp_pair_writer->write(
                            local_buf, 2 * local_buf_filled);
                      }
                      local_buf_filled = 0;
                    }

                    if (phi_pair_reader_1->empty() == false &&
                        (std::uint64_t)(phi_pair_reader_1->peek().m_key) <
                        range_end) {
                      pair_type pp = phi_pair_reader_1->read();
                      std::uint64_t next_i = (std::uint64_t)pp.m_key;
                      std::uint64_t next_phi_i = (std::uint64_t)pp.m_val;
                      lcp_1 = std::max(0L,
                          (std::int64_t)(i_1 + lcp_1) -
                          (std::int64_t)next_i);
                      i_1 = next_i;
                      phi_i_1 = next_phi_i;
                    } else i_1 = text_length;
                  } else break;
                }

                while (i_2 != text_length) {
                  while (i_2 + lcp_2 < text_length &&
                      phi_i_2 + lcp_2 < buf_beg + buf_filled) {
                    char_type next_char =
                      (i_2 + lcp_2 < ext_segment_end) ?
                      segment[(i_2 + lcp_2) - segment_beg] :
                      accessor_2->access(i_2 + lcp_2);
                    if (next_char == buffer[(phi_i_2 + lcp_2) - buf_beg])
                      ++lcp_2;
                    else break;
                  }

                  if (phi_i_2 + lcp_2 < buf_beg + buf_filled ||
                      i_2 + lcp_2 == text_length ||
                      buf_beg + buf_filled == text_length) {
                    local_buf[2 * local_buf_filled] = (text_offset_type)i_2;
                    local_buf[2 * local_buf_filled + 1] =
                      (text_offset_type)lcp_2;
                    ++local_buf_filled;

                    if (local_buf_filled == local_buf_size) {
                      #pragma omp critical
                      {
                        lcp_pair_writer->write(
                            local_buf, 2 * local_buf_filled);
                      }
                      local_buf_filled = 0;
                    }

                    if (phi_pair_reader_2->empty() == false &&
                        (std::uint64_t)(phi_pair_reader_2->peek().m_key) <
                        range_end) {
                      pair_type pp = phi_pair_reader_2->read();
                      std::uint64_t next_phi_i = (std::uint64_t)pp.m_key;
                      std::uint64_t next_i = (std::uint64_t)pp.m_val;
                      lcp_2 = std::max(0L,
                          (std::int64_t)(phi_i_2 + lcp_2) -
                          (std::int64_t)next_phi_i);
                      i_2 = next_i;
                      phi_i_2 = next_phi_i;
                    } else i_2 = text_length;
                  } else break;
                }

                buf_beg += buf_filled;
              }

              if (local_buf_filled > 0) {
                #pragma omp critical
                {
                  lcp_pair_writer->write(local_buf, 2 * local_buf_filled);
                }
              }

              #pragma omp critical
              {
                // Stop I/O threads.
                text_reader->stop_reading();
                phi_pair_reader_1->stop_reading();
                phi_pair_reader_2->stop_reading();

                // Update I/O volume.
                io_vol += accessor_1->bytes_read();
                io_vol += accessor_2->bytes_read();
                io_vol += text_reader->bytes_read();
                io_vol += phi_pair_reader_1->bytes_read();
                io_vol += phi_pair_reader_2->bytes_read();
              }

              // Clean up.
              // XXX not deleted in stack order
              delete phi_pair_reader_2;
              delete phi_pair_reader_1;
              delete accessor_2;
              delete accessor_1;
              delete text_reader;
              utils::deallocate(local_buf);
            }
          }

          // Update I/O volume.
          io_vol +=
            lcp_pair_writer->bytes_written();
          total_io_volume += io_vol;

          // Print summary.
          long double process_time = utils::wclock() - process_start;
          fprintf(stderr, "\r    Process (i, Phi[i]) pairs: time = %.2Lfs, "
              "I/O = %.2LfMiB/s, total I/O vol = %.2Lfbytes/input symbol\n",
              process_time, (1.L * io_vol / (1L << 20)) / process_time,
              (1.L * total_io_volume) / text_length);

          // Clean up.
          delete lcp_pair_writer;
#endif

          // Clean up.
          utils::deallocate(segment);
          utils::file_delete(phi_pairs_filename_1);
          utils::file_delete(phi_pairs_filename_2);
        }
      }

#ifdef _OPENMP

      // Clean up.
      for (std::uint64_t segment_id = 0;
          segment_id < n_segments; ++segment_id)
        delete type_1_bucket_sizes[segment_id];
      delete[] type_1_bucket_sizes;
#endif
    }

#ifdef _OPENMP
    // Clean up.
    utils::deallocate(aux_sparse_plcp);
#endif
  }

  // Allocate sparse PLCP.
  text_offset_type *sparse_plcp =
    utils::allocate_array<text_offset_type>(sparse_plcp_size);

  // Load pairs (i, PLCP[i]) from disk and permute in RAM.
  {

    // Print initial message.
    fprintf(stderr, "  Load and permute pairs (i, PLCP[i]): ");

    // Start the timer.
    long double lcp_permute_start = utils::wclock();

    // Initialize the I/O volume.
    std::uint64_t io_vol = 0;

    // Process segments left-to-right.
    for (std::uint64_t segment_id = 0;
        segment_id < n_segments; ++segment_id) {
      std::string lcp_pairs_filename =
        output_filename + ".lcp_pairs." + utils::intToStr(segment_id);
      std::uint64_t n_pairs =
        utils::file_size(lcp_pairs_filename) / (2 * sizeof(text_offset_type));

#ifdef _OPENMP
      static const std::uint64_t opt_in_buf_ram = (32UL << 20);
      static const std::uint64_t opt_par_buf_ram = (4UL << 20);

      std::uint64_t in_buf_ram = opt_in_buf_ram;
      std::uint64_t par_buf_ram = opt_par_buf_ram;

      {
        std::uint64_t total_buf_ram = in_buf_ram + par_buf_ram;
        if (total_buf_ram + sparse_plcp_ram > ram_use) {
          std::uint64_t ram_budget = ram_use - sparse_plcp_ram;
          long double shrink_factor =
            (long double)ram_budget / (long double)total_buf_ram;
          in_buf_ram =
            (std::uint64_t)((long double)in_buf_ram * shrink_factor);
          par_buf_ram =
            (std::uint64_t)((long double)par_buf_ram * shrink_factor);
        }
      }

      typedef async_stream_reader<text_offset_type> lcp_pair_reader_type;
      lcp_pair_reader_type *lcp_pair_reader = new lcp_pair_reader_type(
          lcp_pairs_filename, in_buf_ram,
          std::max(4UL, in_buf_ram / (2UL << 20)));

      std::uint64_t par_buf_size =
        par_buf_ram / (2 * sizeof(text_offset_type));
      text_offset_type *pair_buffer =
        utils::allocate_array<text_offset_type>(2 * par_buf_size);

      std::uint64_t pairs_read = 0;
      while (pairs_read < n_pairs) {
        std::uint64_t buf_filled =
          std::min(n_pairs - pairs_read, par_buf_size);
        lcp_pair_reader->read(pair_buffer, buf_filled * 2);

        #pragma omp parallel for
        for (std::uint64_t j = 0; j < buf_filled; ++j) {
          std::uint64_t i = pair_buffer[2 * j];
          std::uint64_t lcp = pair_buffer[2 * j + 1];

          // Invariant: PLCP[i] = lcp.
          sparse_plcp[i / plcp_sampling_rate] = lcp;
        }

        pairs_read += buf_filled;
      }

      // Clean up.
      utils::deallocate(pair_buffer);

#else
      std::uint64_t opt_in_buf_ram = (32UL << 20);
      std::uint64_t in_buf_ram = opt_in_buf_ram;
      if (sparse_plcp_ram + in_buf_ram > ram_use) {
        in_buf_ram = ram_use - sparse_plcp_ram;
      }

      typedef async_stream_reader<text_offset_type> lcp_pair_reader_type;
      lcp_pair_reader_type *lcp_pair_reader = new lcp_pair_reader_type(
          lcp_pairs_filename, in_buf_ram,
          std::max(4UL, in_buf_ram / (2UL << 20)));

      for (std::uint64_t j = 0; j < n_pairs; ++j) {
        std::uint64_t i = (std::uint64_t)lcp_pair_reader->read();
        std::uint64_t lcp = (std::uint64_t)lcp_pair_reader->read();

        // Invariant: PLCP[i] = lcp.
        sparse_plcp[i / plcp_sampling_rate] = lcp;
      }
#endif

      // Stop I/O threads.
      lcp_pair_reader->stop_reading();

      // Update I/O volume.
      io_vol +=
        lcp_pair_reader->bytes_read();

      // Clean up.
      delete lcp_pair_reader;
      utils::file_delete(lcp_pairs_filename);
    }

    // Handle special case.
    if ((phi_undefined_position % plcp_sampling_rate) == 0)
      sparse_plcp[phi_undefined_position / plcp_sampling_rate] = 0;

    // Update I/O volume.
    total_io_volume += io_vol;

    // Print summary.
    long double lcp_pair_permute_time = utils::wclock() - lcp_permute_start;
    fprintf(stderr, "\r  Load and permute pairs (i, PLCP[i]): "
        "time = %.1Lfs, I/O = %.2LfMiB/s, "
        "total I/O vol = %.2Lfbytes/input symbol\n",
        lcp_pair_permute_time,
        ((1.L * io_vol) / (1L << 20)) / lcp_pair_permute_time,
        (1.L * total_io_volume) / text_length);
  }

  // Print summary.
  long double total_time = utils::wclock() - start;
  fprintf(stderr, "  Summary: time = %.2Lfs, "
      "total I/O vol = %.2Lfbytes/input symbol\n\n",
      total_time, (1.L * total_io_volume) / text_length);

  // Return the result.
  return sparse_plcp;
}

}  // namespace em_sparse_phi_private

#endif  // __SRC_EM_SPARSE_PHI_SRC_COMPUTE_SPARSE_PHI_HPP_INCLUDED
