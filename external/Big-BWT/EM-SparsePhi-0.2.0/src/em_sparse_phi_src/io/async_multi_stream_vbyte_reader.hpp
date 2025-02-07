/**
 * @file    src/em_sparse_phi_src/io/async_multi_stream_vbyte_reader.hpp
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

#ifndef __SRC_EM_SPARSE_PHI_SRC_IO_ASYNC_MULTI_STREAM_VBYTE_READER_HPP_INCLUDED
#define __SRC_EM_SPARSE_PHI_SRC_IO_ASYNC_MULTI_STREAM_VBYTE_READER_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <string>
#include <algorithm>
#include <thread>
#include <mutex>
#include <condition_variable>

#include "../utils.hpp"


namespace em_sparse_phi_private {

class async_multi_stream_vbyte_reader {
  private:
    template<typename T>
    struct buffer {
      buffer(std::uint64_t size, T* const mem)
        : m_content(mem), m_size(size) {
        m_filled = 0;
        m_is_filled = false;
      }

      void read_from_file(std::FILE *f) {
        utils::read_from_file(m_content, m_size, m_filled, f);
      }

      inline std::uint64_t size_in_bytes() const {
        return sizeof(T) * m_filled;
      }

      T* const m_content;
      const std::uint64_t m_size;

      std::uint64_t m_filled;
      bool m_is_filled;
    };

    template<typename T>
    struct circular_queue {
      private:
        std::uint64_t m_size;
        std::uint64_t m_filled;
        std::uint64_t m_head;
        std::uint64_t m_tail;
        T *m_data;

      public:
        circular_queue()
          : m_size(1),
            m_filled(0),
            m_head(0),
            m_tail(0),
            m_data(new T[m_size]) {}

        inline void push(T x) {
          m_data[m_head++] = x;
          if (m_head == m_size)
            m_head = 0;
          ++m_filled;
          if (m_filled == m_size)
            enlarge();
        }

        inline T &front() const {
          return m_data[m_tail];
        }

        inline void pop() {
          ++m_tail;
          if (m_tail == m_size)
            m_tail = 0;
          --m_filled;
        }

        inline bool empty() const {
          return (m_filled == 0);
        }

        inline std::uint64_t size() const {
          return m_filled;
        }

        ~circular_queue() {
          delete[] m_data;
        }

      private:
        void enlarge() {
          T *new_data = new T[2 * m_size];
          std::uint64_t left = m_filled;
          m_filled = 0;

          while (left > 0) {
            std::uint64_t tocopy = std::min(left, m_size - m_tail);
            std::copy(m_data + m_tail,
                m_data + m_tail + tocopy, new_data + m_filled);

            m_tail += tocopy;
            if (m_tail == m_size)
              m_tail = 0;
            left -= tocopy;
            m_filled += tocopy;
          }

          m_head = m_filled;
          m_tail = 0;
          m_size <<= 1;
          std::swap(m_data, new_data);
          delete[] new_data;
        }
    };

    template<typename buffer_type>
    struct request {
      request() {}
      request(buffer_type *buffer, std::uint64_t file_id) {
        m_buffer = buffer;
        m_file_id = file_id;
      }

      buffer_type *m_buffer;
      std::uint64_t m_file_id;
    };

    template<typename request_type>
    struct request_queue {
      request_queue()
        : m_no_more_requests(false) {}

      request_type get() {
        request_type ret = m_requests.front();
        m_requests.pop();
        return ret;
      }

      inline void add(request_type request) {
        std::lock_guard<std::mutex> lk(m_mutex);
        m_requests.push(request);
      }

      inline bool empty() const {
        return m_requests.empty();
      }

      circular_queue<request_type> m_requests;  // Must have FIFO property
      std::condition_variable m_cv;
      std::mutex m_mutex;
      bool m_no_more_requests;
    };

  private:
    template<typename T>
    static void async_io_thread_code(async_multi_stream_vbyte_reader *caller) {
      typedef buffer<T> buffer_type;
      typedef request<buffer_type> request_type;
      while (true) {

        // Wait for request or until 'no more requests' flag is set.
        std::unique_lock<std::mutex> lk(caller->m_read_requests.m_mutex);
        while (caller->m_read_requests.empty() &&
            !(caller->m_read_requests.m_no_more_requests))
          caller->m_read_requests.m_cv.wait(lk);

        if (caller->m_read_requests.empty() &&
            caller->m_read_requests.m_no_more_requests) {

          // No more requests -- exit.
          lk.unlock();
          break;
        }

        // Extract the buffer from the collection.
        request_type request = caller->m_read_requests.get();
        lk.unlock();

        // Process the request.
        request.m_buffer->read_from_file(caller->m_files[request.m_file_id]);
        caller->m_bytes_read += request.m_buffer->size_in_bytes();

        // Update the status of the buffer
        // and notify the waiting thread.
        std::unique_lock<std::mutex> lk2(caller->m_mutexes[request.m_file_id]);
        request.m_buffer->m_is_filled = true;
        lk2.unlock();
        caller->m_cvs[request.m_file_id].notify_one();
      }
    }

  private:
    typedef buffer<std::uint8_t> buffer_type;
    typedef request<buffer_type> request_type;

    std::uint64_t m_bytes_read;
    std::uint64_t m_items_per_buf;
    std::uint64_t n_files;
    std::uint64_t m_files_added;

    std::FILE **m_files;
    std::uint64_t *m_active_buffer_pos;
    std::uint8_t *m_mem;
    buffer_type **m_active_buffers;
    buffer_type **m_passive_buffers;
    std::mutex *m_mutexes;
    std::condition_variable *m_cvs;

    request_queue<request_type> m_read_requests;
    std::thread *m_io_thread;
    std::uint64_t m_io_thread_count;

  private:
    void issue_read_request(std::uint64_t file_id) {
      request_type req(m_passive_buffers[file_id], file_id);
      m_read_requests.add(req);
      m_read_requests.m_cv.notify_one();
    }

    void receive_new_buffer(std::uint64_t file_id) {

      // Wait for the I/O thread to finish reading passive buffer.
      std::unique_lock<std::mutex> lk(m_mutexes[file_id]);
      while (m_passive_buffers[file_id]->m_is_filled == false)
        m_cvs[file_id].wait(lk);

      // Swap active and passive buffers.
      std::swap(m_active_buffers[file_id], m_passive_buffers[file_id]);
      m_active_buffer_pos[file_id] = 0;
      m_passive_buffers[file_id]->m_is_filled = false;
      lk.unlock();

      // Issue the read request for the passive buffer.
      issue_read_request(file_id);
    }

  public:
    // Constructor, takes the number of files and a
    // size of per-file buffer (in bytes) as arguments.
    async_multi_stream_vbyte_reader(
        std::uint64_t number_of_files,
        std::uint64_t buf_size_bytes = (std::uint64_t)(1 << 20)) {

      // Sanity check.
      if (number_of_files == 0) {
        fprintf(stderr, "\nError in async_multi_stream_vbyte_reader: "
            "number_of_files == 0\n");
        std::exit(EXIT_FAILURE);
      }

      // Initialize basic parameters.
      n_files = number_of_files;
      m_files_added = 0;
      m_bytes_read = 0;

      // Computer optimal buffer size.
      buf_size_bytes = std::max((std::uint64_t)1, buf_size_bytes / 2);
      m_items_per_buf = utils::disk_block_size<std::uint8_t>(buf_size_bytes);

      // Allocate arrays storing info about each file.
      m_mutexes = new std::mutex[n_files];
      m_cvs = new std::condition_variable[n_files];
      m_active_buffer_pos = new std::uint64_t[n_files];
      m_files = new std::FILE*[n_files];
      m_active_buffers = new buffer_type*[n_files];
      m_passive_buffers = new buffer_type*[n_files];

      // Allocate buffers.
      std::uint64_t toallocate = 2 * n_files * m_items_per_buf;
      m_mem = utils::allocate_array<std::uint8_t>(toallocate);
      {
        std::uint8_t *mem = m_mem;
        for (std::uint64_t i = 0; i < n_files; ++i) {
          m_active_buffer_pos[i] = 0;
          m_active_buffers[i] = new buffer_type(m_items_per_buf, mem);
          mem += m_items_per_buf;
          m_passive_buffers[i] = new buffer_type(m_items_per_buf, mem);
          mem += m_items_per_buf;
        }
      }

      // Start the I/O thread.
      m_io_thread = new std::thread(async_io_thread_code<std::uint8_t>, this);
    }

    // The added file gets the next available ID (starting from 0).
    void add_file(std::string filename) {
      m_files[m_files_added] = utils::file_open_nobuf(filename, "r");
      issue_read_request(m_files_added);
      ++m_files_added;
    }

  private:
    // Read one byte from i-th file.
    std::uint8_t read_byte_from_ith_file(std::uint64_t i) {
      if (m_active_buffer_pos[i] == m_active_buffers[i]->m_filled)
        receive_new_buffer(i);
      return m_active_buffers[i]->m_content[m_active_buffer_pos[i]++];
    }

  public:
    // Read the next item from i-th file.
    std::uint64_t read_from_ith_file(std::uint64_t i) {
      std::uint64_t result = 0;
      std::uint64_t offset = 0;
      std::uint64_t next_char = read_byte_from_ith_file(i);
      while (next_char & 0x80) {
        result |= ((next_char & 0x7f) << offset);
        offset += 7;
        next_char = read_byte_from_ith_file(i);
      }
      result |= (next_char << offset);
      return result;
    }

    // Return performed I/O in bytes.
    inline std::uint64_t bytes_read() const {
      return m_bytes_read;
    }

    // Stop the I/O thread, now the user can
    // cafely call the bytes_read() method.
    void stop_reading() {
      if (m_io_thread != NULL) {
        std::unique_lock<std::mutex> lk(m_read_requests.m_mutex);
        m_read_requests.m_no_more_requests = true;
        lk.unlock();
        m_read_requests.m_cv.notify_one();
        m_io_thread->join();
        delete m_io_thread;
        m_io_thread = NULL;
      }
    }

    // Destructor.
    ~async_multi_stream_vbyte_reader() {
      stop_reading();

      // Delete buffers and close files.
      for (std::uint64_t i = n_files; i > 0; --i) {
        std::fclose(m_files[i - 1]);
        delete m_passive_buffers[i - 1];
        delete m_active_buffers[i - 1];
      }

      // Rest of the cleanup.
      utils::deallocate(m_mem);
      delete[] m_passive_buffers;
      delete[] m_active_buffers;
      delete[] m_files;
      delete[] m_active_buffer_pos;
      delete[] m_cvs;
      delete[] m_mutexes;
    }
};

}  // namespace em_sparse_phi_private

#endif  // __SRC_EM_SPARSE_PHI_SRC_IO_ASYNC_MULTI_STREAM_VBYTE_READER_HPP_INCLUDED
