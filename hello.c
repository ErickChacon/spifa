#include "R.h"

void hello(int *n) {
  int i;
  for (i = 0; i < *n; ++i) {
    Rprintf("hello %d \n", i);
  }
  *n += 1;
}

