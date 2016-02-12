
#ifndef LONGMEM_SIZE_H_
#define LONGMEM_SIZE_H_

#include <stdint.h>
// Copyright Peter Andrews 2013 CSHL

// Signed integer type
#ifdef SINTS
typedef int64_t SINT;
#define SINT_MAX LONG_MAX
#else
typedef int32_t SINT;
#define SINT_MAX INT_MAX
#endif

// Unsigned integer type
#ifdef UINTS
typedef uint64_t ANINT;
#else
typedef uint32_t ANINT;
#endif

#endif  // LONGMEM_SIZE_H_
