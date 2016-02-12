/* Copyright Peter Andrews 2013 CSHL */

#ifndef LONGMEM_LOCKED_H_
#define LONGMEM_LOCKED_H_

#include <pthread.h>
#include <stdarg.h>
#include <stdint.h>

#include <cstdio>
#include <vector>

#include "./util.h"
#include "./error.h"

const bool checked = false;

inline void lock(pthread_mutex_t * const mutex) {
  if (checked) {
    if (pthread_mutex_lock(mutex)) throw paa::Error("lock error");
  } else {
    pthread_mutex_lock(mutex);
  }
}
inline void unlock(pthread_mutex_t * const mutex) {
  if (checked) {
    if (pthread_mutex_unlock(mutex)) throw paa::Error("unlock error");
  } else {
    pthread_mutex_unlock(mutex);
  }
}

template<class T>
class RingBuffer {
 public:
  RingBuffer(const uint64_t N_ = 1000, const T & t = T(),
             bool use_lock_ = true, bool use_space_ = true,
             bool use_ready_ = true) :
      N(N_), begin(0), end(0), buffer(N, t),
      use_lock(use_lock_), use_space(use_space_), use_ready(use_ready_) {
    // printf("created ring buffer of size %ld\n", N);
    init();
  }
  RingBuffer(const RingBuffer & other) :
      N(other.N), begin(other.begin), end(other.end), buffer(other.buffer),
      use_lock(other.use_lock), use_space(other.use_space),
      use_ready(other.use_ready) {
    // printf("copied ring buffer of size %ld\n", N);
    init();
  }
  ~RingBuffer() {
    if (checked) {
      lock();
      if (end != begin)
        throw paa::Error("buffer was not empty when finished")
            << begin << " " << end;
      unlock();
    }
    buffer.clear();
    if (checked) {
      if (use_lock && pthread_mutex_destroy(&mutex))
        throw paa::Error("mutex destroy error");
      if (use_space && pthread_cond_destroy(&space))
        throw paa::Error("space cond destroy error");
      if (use_ready && pthread_cond_destroy(&ready))
        throw paa::Error("ready cond destroy error");
    } else {
      if (use_lock) pthread_mutex_destroy(&mutex);
      if (use_space) pthread_cond_destroy(&space);
      if (use_ready) pthread_cond_destroy(&ready);
    }
  }
  void push(const T & val) {
    lock();
    if (use_space) while (end - begin == N) pthread_cond_wait(&space, &mutex);
    buffer[end++ % N] = val;
    if (use_ready) pthread_cond_signal(&ready);
    unlock();
  }
  const T pop() {
    lock();
    if (use_ready) while (end - begin == 0) pthread_cond_wait(&ready, &mutex);
    const T val = buffer[begin++ % N];
    if (use_space) pthread_cond_signal(&space);
    unlock();
    return val;
  }
  T * next(const unsigned int extra = 1) {
    T * n = nullptr;
    lock();
    if (checked && end < begin)
      throw paa::Error("unexpected begin") << begin << "or end"
                                      << end << "or extra" << extra;
    if (use_space) {
      while (end - begin + extra >= N) pthread_cond_wait(&space, &mutex);
      n = &buffer[end % N];
    } else if (end - begin + extra < N) {
      n = &buffer[end % N];
    }
    unlock();
    return n;
  }
  void commit() {
    lock();
    ++end;
    if (use_ready) pthread_cond_signal(&ready);
    unlock();
  }
  T * first() {
    T * f = nullptr;
    lock();
    if (use_ready) {
      while (end == begin) pthread_cond_wait(&ready, &mutex);
      f = &buffer[begin % N];
    } else if (end != begin) {
      f = &buffer[begin % N];
    }
    unlock();
    return f;
  }
  void yield() {
    lock();
    ++begin;
    if (use_space) pthread_cond_signal(&space);
    unlock();
  }

 private:
  void init() {
    if (checked) {
      if (pthread_mutex_init(&mutex, nullptr))
        throw paa::Error("mutex init error");
      if (use_space && pthread_cond_init(&space, nullptr))
        throw paa::Error("space cond init error");
      if (use_ready && pthread_cond_init(&ready, nullptr))
        throw paa::Error("ready cond init error");
    } else {
      pthread_mutex_init(&mutex, nullptr);
      if (use_space) pthread_cond_init(&space, nullptr);
      if (use_ready) pthread_cond_init(&ready, nullptr);
    }
  }
  void lock() { if (use_lock) ::lock(&mutex); }
  void unlock() { if (use_lock) ::unlock(&mutex); }
  uint64_t N;
  uint64_t begin;
  uint64_t end;
  std::vector<T> buffer;
  bool use_lock;
  bool use_space;
  bool use_ready;
  pthread_mutex_t mutex;
  pthread_cond_t space;
  pthread_cond_t ready;
  RingBuffer & operator=(const RingBuffer & disabled_assignment_operator);
};

class OutputBuffer {
 public:
  OutputBuffer(const uint64_t max_size_ = 100000,
               const uint64_t max_line_size_ = 10000,
               FILE * out_ = stdout) :
      max_size(max_size_), max_line_size(max_line_size_), end(0),
      buffer(max_line_size), out(out_) {}
  void set_max_line_size(const uint64_t max_line_size_) {
    max_line_size = max_line_size_;
  }
  void prepare_for_line(const uint64_t expected_line_size) {
    const uint64_t needed_size = end + expected_line_size;
    if (needed_size > buffer.size()) buffer.resize(needed_size);
  }
  void flush() {
    if (end) {
      lock(&mutex);
      fwrite(&buffer[0], sizeof(buffer[0]), end, stdout);
      unlock(&mutex);
      end = 0;
    }
  }
  void printf(const char * format, ...) {
    const uint64_t needed_size = end + max_line_size;
    if (needed_size > buffer.size()) buffer.resize(needed_size);
    va_list args;
    va_start(args, format);
    if (checked) {
      const int64_t written = vsprintf(&buffer[end], format, args);
      if (written < 0) throw paa::Error("vsprintf output error");
      if (written >= (int64_t)max_line_size)
        throw paa::Error("too much written to output buffer");
      end += written;
      if (end >= buffer.size())
        throw paa::Error("too much written to output buffer");
    } else {
      end += vsprintf(&buffer[end], format, args);
    }
    va_end(args);
    if (needed_size > max_size && end && buffer[end - 1] == '\n')
      flush();
  }

 private:
  uint64_t max_size;
  uint64_t max_line_size;
  uint64_t end;
  std::vector<char> buffer;
  FILE * out;
  static pthread_mutex_t mutex;
};

#endif  // LONGMEM_LOCKED_H_
