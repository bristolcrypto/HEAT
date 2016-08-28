#pragma once

#ifdef BINARY
#define PLAINTEXT_MODULUS 2
#define PLAINTEXT_MODULUS_MPZ "2"
#else 
#ifdef MPZ
#define PLAINTEXT_MODULUS 0
#define PLAINTEXT_MODULUS_MPZ "987543987543987543987543987543987543"
#else // ULONG
#define PLAINTEXT_MODULUS 28276
#define PLAINTEXT_MODULUS_MPZ "28276"
#endif
#endif

#define SECURITY_LEVEL 80
#define MULTIPLICATIVE_DEPTH 1

// Definition of the polynomial type (with 64-bit integers, 16 coefficients, a
// (5*62)-bit modulus)
using poly_t = nfl::poly_from_modulus<uint64_t, 1 << 4, 620>;

// Auto-generated poly type to enable multiplication over ZZ
// (twice the size (+1) of the polynomials we manipulate)
using polyZ_t =
    nfl::poly<poly_t::value_type, poly_t::degree, poly_t::nmoduli * 2 + 1>;

// FHE types
using sk_t = SecretKey<poly_t>;
using evk_t = EvaluationKey<poly_t, polyZ_t>;
using pk_t = PublicKey<poly_t, polyZ_t>;
using ciphertext_t = Ciphertext<poly_t, polyZ_t>;
