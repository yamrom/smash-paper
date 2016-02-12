/* Copyright Peter Andrews 2013 CSHL */

#include "./locked.h"

#pragma GCC diagnostic ignored "-Wzero-as-null-pointer-constant"
pthread_mutex_t OutputBuffer::mutex = PTHREAD_MUTEX_INITIALIZER;
#pragma GCC diagnostic pop

