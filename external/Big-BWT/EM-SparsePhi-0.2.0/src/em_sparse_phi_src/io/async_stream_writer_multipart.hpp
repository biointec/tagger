/**
 * @file    src/em_sparse_phi_src/io/async_stream_writer_multipart.hpp
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

#ifndef __SRC_EM_SPARSE_PHI_SRC_IO_ASYNC_STREAM_WRITER_MULTIPART_HPP_INCLUDED
#define __SRC_EM_SPARSE_PHI_SRC_IO_ASYNC_STREAM_WRITER_MULTIPART_HPP_INCLUDED

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

template<typename value_type>
class async_stream_writer_multipart {
  private:
    template<typename T>
    struct buffer {
      buffer(std::uint64_t size, T* const mem)
        : m_content(mem), m_size(size) {
        m_filled = 0;
      }

      void write_to_file(std::FILE *f) {
        utils::write_to_file(m_content, m_filled, f);
        m_filled = 0;
      }

      inline bool empty() const {
        return (m_filled == 0);
      }

      inline bool full() const {
        return (m_filled == m_size);
      }

      inline std::uint64_t size_in_bytes() const {
        return sizeof(T) * m_filled;
      }

      inline std::uint64_t free_space() const {
        return m_size - m_filled;
      }

      T* const m_content;
      const std::uint64_t m_size;

      std::uint64_t m_filled;
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

    enum signal_type {
      k_no_signal,
      k_flush_signal,
      k_stop_signal,
      k_end_file_signal
    };

    template<typename T>
    struct buffer_queue {
      typedef buffer<T> buffer_type;

      buffer_queue(
          std::uint64_t n_buffers,
          std::uint64_t items_per_buf,
          T *mem) {
        for (std::uint64_t i = 0; i < n_buffers; ++i) {
          m_queue.push(new buffer_type(items_per_buf, mem));
          mem += items_per_buf;
        }
      }

      ~buffer_queue() {
        while (!m_queue.empty()) {
          buffer_type *buf = m_queue.front();
          m_queue.pop();
          delete buf;
        }
      }

      buffer_type *pop() {
        buffer_type *ret = m_queue.front();
        m_queue.pop();
        return ret;
      }

      void push(buffer_type *buf) {
        std::lock_guard<std::mutex> lk(m_mutex);
        m_signals.push(k_flush_signal);
        m_queue.push(buf);
      }

      void send_stop_signal() {
        std::lock_guard<std::mutex> lk(m_mutex);
        m_signals.push(k_stop_signal);
      }

      void send_end_file_signal() {
        std::lock_guard<std::mutex> lk(m_mutex);
        m_signals.push(k_end_file_signal);
      }

      inline bool empty() const {
        return m_queue.empty();
      }

      circular_queue<buffer_type*> m_queue;
      circular_queue<std::uint64_t> m_signals;
      std::condition_variable m_cv;
      std::mutex m_mutex;
    };

  private:
    typedef buffer<value_type> buffer_type;
    typedef buffer_queue<value_type> buffer_queue_type;

    buffer_queue_type *m_empty_buffers;
    buffer_queue_type *m_full_buffers;

  private:
    template<typename T>
    static void io_thread_code(async_stream_writer_multipart<T> *caller) {
      typedef buffer<T> buffer_type;
      while (true) {

        // Wait for the full buffer or a signal.
        std::unique_lock<std::mutex> lk(caller->m_full_buffers->m_mutex);
        while (caller->m_full_buffers->m_signals.empty())
          caller->m_full_buffers->m_cv.wait(lk);

        // Extract the singla.
        std::uint64_t signal = caller->m_full_buffers->m_signals.front();
        caller->m_full_buffers->m_signals.pop();

        // Take action depending on the signal type.
        if (signal == k_stop_signal) {

          // Release the lock and terminate the threads.
          lk.unlock();
          break;
        } else if (signal == k_end_file_signal) {

          // Release the lock and
          // close the current file.
          lk.unlock();
          if (caller->m_file != NULL) {
            std::fclose(caller->m_file);
            caller->m_file = NULL;
            ++caller->m_cur_part;
          }
        } else {

          // Extract the buffer from the collection.
          buffer_type *buffer = caller->m_full_buffers->pop();
          lk.unlock();

          // Open new file if needed.
          if (caller->m_file == NULL) {
            std::string cur_part_filename = caller->m_filename +
              ".part" + utils::intToStr(caller->m_cur_part);
            caller->m_file = utils::file_open_nobuf(cur_part_filename, "w");
          }

          // Write data to disk.
          buffer->write_to_file(caller->m_file);

          // Add the (now empty) buffer to the collection
          // of empty buffers and notify the waiting thread.
          caller->m_empty_buffers->push(buffer);
          caller->m_empty_buffers->m_cv.notify_one();
        }
      }
    }

    // Get a free buffer from the collection of free buffers.
    buffer_type* get_empty_buffer() {
      std::unique_lock<std::mutex> lk(m_empty_buffers->m_mutex);
      while (m_empty_buffers->empty())
        m_empty_buffers->m_cv.wait(lk);
      buffer_type *ret = m_empty_buffers->pop();
      lk.unlock();
      return ret;
    }

  private:
    std::FILE *m_file;
    std::string m_filename;

    std::uint64_t m_cur_part;
    std::uint64_t m_bytes_written;
    std::uint64_t m_items_per_buf;

    value_type *m_mem;
    buffer_type *m_cur_buffer;
    std::thread *m_io_thread;

  public:
    async_stream_writer_multipart(std::string filename) {
      init(filename, (8UL << 20), 4UL);
    }

    async_stream_writer_multipart(std::string filename,
        std::uint64_t total_buf_size_bytes,
        std::uint64_t n_buffers) {
      init(filename, total_buf_size_bytes, n_buffers);
    }

    // Main initializing function.
    void init(std::string filename,
        std::uint64_t total_buf_size_bytes,
        std::uint64_t n_buffers) {

      // Sanity check.
      if (n_buffers == 0) {
        fprintf(stderr, "\nError in async_stream_writer_multipart: "
            "n_buffers == 0\n");
        std::exit(EXIT_FAILURE);
      }

      m_filename = filename;
      m_cur_part = 0;
      m_file = NULL;

      // Computer optimal buffer size.
      std::uint64_t buf_size_bytes =
        std::max((std::uint64_t)1, total_buf_size_bytes / n_buffers);
      m_items_per_buf = utils::disk_block_size<value_type>(buf_size_bytes);

      // Allocate buffers.
      m_mem = utils::allocate_array<value_type>(n_buffers * m_items_per_buf);
      m_empty_buffers = new buffer_queue_type(n_buffers,
          m_items_per_buf, m_mem);
      m_full_buffers = new buffer_queue_type(0, 0, NULL);

      // Initialize empty buffer.
      m_cur_buffer = get_empty_buffer();
      m_bytes_written = 0;

      // Start the I/O thread.
      m_io_thread = new std::thread(io_thread_code<value_type>, this);
    }

    // Write given item to the stream.
    inline void write(value_type value) {

      // We count I/O volume here (and not in the thread doing I/O) to
      // avoid the situation, where user call bytes_written(), but the
      // I/O thread is still writing the last buffer.
      m_bytes_written += sizeof(value_type);
      if (m_cur_buffer->full())
        flush();

      m_cur_buffer->m_content[m_cur_buffer->m_filled++] = value;
    }

    // Write values[0..length) to the stream.
    inline void write(const value_type *values, std::uint64_t length) {
      m_bytes_written += length * sizeof(value_type);
      while (length > 0) {
        if (m_cur_buffer->full())
          flush();

        std::uint64_t tocopy = std::min(length, m_cur_buffer->free_space());
        std::copy(values, values + tocopy,
            m_cur_buffer->m_content + m_cur_buffer->m_filled);
        m_cur_buffer->m_filled += tocopy;
        values += tocopy;
        length -= tocopy;
      }
    }

    // Return performed I/O in bytes.
    inline std::uint64_t bytes_written() const {
      return m_bytes_written;
    }

    // It's safe to call if the buffer is not full, though
    // in principle should only be called internally. Calling
    // it too often will lead to poor I/O performance.
    void flush() {
      if (!m_cur_buffer->empty()) {
        m_full_buffers->push(m_cur_buffer);
        m_full_buffers->m_cv.notify_one();
        m_cur_buffer = get_empty_buffer();
      }
    }

    void end_current_file() {

      // Flush current buffer.
      flush();

      // Let the I/O thread know that
      // we want to end current file.
      m_full_buffers->send_end_file_signal();
      m_full_buffers->m_cv.notify_one();
    }

    ~async_stream_writer_multipart() {
      end_current_file();

      // Let the I/O thread know that we're done.
      m_full_buffers->send_stop_signal();
      m_full_buffers->m_cv.notify_one();

      // Wait for the I/O thread to finish.
      m_io_thread->join();

      // Clean up.
      delete m_empty_buffers;
      delete m_full_buffers;
      delete m_io_thread;
      delete m_cur_buffer;
      utils::deallocate(m_mem);
    }
};

}  // namespace em_sparse_phi_private

#endif  // __SRC_EM_SPARSE_PHI_SRC_IO_ASYNC_STREAM_WRITER_MULTIPART_HPP_INCLUDED
