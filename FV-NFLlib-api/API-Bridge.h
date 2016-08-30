#pragma once

#include "../API.h"
#include "params.h"
#include <nfl.hpp>

namespace FV {
namespace params {
// Definition of the polynomial type (with 64-bit integers, 16 coefficients, a (5*62)-bit modulus).
using poly_t = nfl::poly_from_modulus<uint64_t, 1 << 4, 620>;
template <typename T>
struct plaintextModulus;
template <>
struct plaintextModulus<mpz_class> {
  static mpz_class value() { return mpz_class(PLAINTEXT_MODULUS_MPZ); }
};
using gauss_struct = nfl::gaussian<uint16_t, uint64_t, 2>;
using gauss_t = nfl::FastGaussianNoise<uint16_t, uint64_t, 2>;
gauss_t fg_prng_sk(8.0, 128, 1 << 14);
gauss_t fg_prng_evk(8.0, 128, 1 << 14);
gauss_t fg_prng_pk(8.0, 128, 1 << 14);
gauss_t fg_prng_enc(8.0, 128, 1 << 14);
}
}  // namespace FV::params


#include "FV.hpp"
