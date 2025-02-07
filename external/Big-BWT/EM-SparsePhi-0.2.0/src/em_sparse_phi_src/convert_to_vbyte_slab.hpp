/**
 * @file    src/em_sparse_phi_src/convert_to_vbyte_slab.hpp
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

#ifndef __SRC_EM_SPARSE_PHI_SRC_CONVERT_TO_VBYTE_SLAB_HPP_INCLUDED
#define __SRC_EM_SPARSE_PHI_SRC_CONVERT_TO_VBYTE_SLAB_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <omp.h>


namespace em_sparse_phi_private {

// Encode tab[0..length) using v-byte encoding and write to dest in
// parallel. We assume that dest is sufficiently large to hold the
// output. The function returns the length of the slab.
template<typename T>
std::uint64_t convert_to_vbyte_slab(
    const T *tab,
    std::uint64_t length,
    std::uint8_t *dest) {

#ifdef _OPENMP
  std::uint64_t max_threads = omp_get_max_threads();
  std::uint64_t max_block_size = (length + max_threads - 1) / max_threads;
  std::uint64_t *slab_ptr = new std::uint64_t[max_threads + 1];

  #pragma omp parallel num_threads(max_threads)
  {
    std::uint64_t thread_id = omp_get_thread_num();
    std::uint64_t block_beg = max_block_size * thread_id;
    std::uint64_t block_end = std::min(block_beg + max_block_size, length);
    std::uint64_t slab_length = 0;
    for (std::uint64_t j = block_beg; j < block_end; ++j) {
      std::uint64_t x = tab[j];
      while (x > 127) {
        ++slab_length;
        x >>= 7;
      }
      ++slab_length;
    }
    slab_ptr[thread_id] = slab_length;

    #pragma omp barrier
    #pragma omp single
    {
      std::uint64_t total_slab_length = 0;
      for (std::uint64_t j = 0; j <= max_threads; ++j) {
        long temp = slab_ptr[j];
        slab_ptr[j] = total_slab_length;
        total_slab_length += temp;
      }
    }

    std::uint64_t dest_ptr = slab_ptr[thread_id];
    for (std::uint64_t j = block_beg; j < block_end; ++j) {
      std::uint64_t x = tab[j];
      while (x > 127) {
        dest[dest_ptr++] = ((x & 0x7f) | 0x80);
        x >>= 7;
      }
      dest[dest_ptr++] = x;
    }
  }

  std::uint64_t total_slab_length = slab_ptr[max_threads];
  delete[] slab_ptr;
  return total_slab_length;

#else
  std::uint64_t slab_length = 0;
  for (std::uint64_t j = 0; j < length; ++j) {
    std::uint64_t x = tab[j];
    while (x > 127) {
      dest[slab_length++] = ((x & 0x7f) | 0x80);
      x >>= 7;
    }
    dest[slab_length++] = x;
  }
  return slab_length;
#endif
}

}  // namespace em_sparse_phi_private

#endif  // __SRC_EM_SPARSE_PHI_SRC_CONVERT_TO_VBYTE_SLAB_HPP_INCLUDED
