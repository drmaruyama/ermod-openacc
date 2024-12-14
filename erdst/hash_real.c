#include <stdint.h>
#include <stddef.h>
#include <string.h>

void hash_double_(const double *v_, const int *elms_, uint64_t *o_)
{
  uint64_t x = 0, buf;
  size_t i, elms;
  elms = (size_t)*elms_;
  for(i = 0; i < elms; ++i){
    memcpy(&buf, &v_[i], sizeof(double));
    x = ((x << 7) | (x >> 57)) ^ buf;
  }
  *o_ = x;
}

void hash_float_(const float *v_, const int *elms_, uint64_t *o_)
{
  uint64_t x = 0;
  uint32_t buf;
  size_t i, elms;
  elms = (size_t)*elms_;
  for(i = 0; i < elms; ++i){
    memcpy(&buf, &v_[i], sizeof(float));
    x = ((x << 7) | (x >> 57)) ^ (uint64_t) buf;
  }
  *o_ = x;
}
