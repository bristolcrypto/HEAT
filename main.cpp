#include <typeinfo>
#include <gmpxx.h>

#include "API.h"
#include "tests.hpp"

int main() {
  int success = 0;

#ifndef FXPT
  #ifdef MPZ
    success |= testInteger_mpz<10>();
  #else
    success |= testInteger<10>();
    success |= testInteger_polynomial<10>();
    success |= testInteger_vector<10>();
  #endif
#else
    success |= testFixedpt<10>();
    success |= testInteger<10>();
    success |= testInteger_polynomial<10>();
#endif

  return 0;
}
