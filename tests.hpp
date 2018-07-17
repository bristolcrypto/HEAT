#ifndef HE_TESTS
#define HE_TESTS

#include <cassert>
#include <iostream>
#include "timing.h"
#include "API.h"
#include <cmath>

template <unsigned long nb_tests>
int testInteger() {
  // Determine if the test should run
  #ifdef BINARY
  unsigned long run = HE::supports_bit_encryption;
  #else
  unsigned long run = HE::supports_unsigned_encryption;
  #endif

  if (!run) {
    return 0;
  }

  timing t;

  gmp_randclass prng(gmp_randinit_default);
//  prng.seed(0);
  prng.seed(time(NULL)); // To set different seeds

  void* parameters = nullptr;
  void* sk = nullptr;
  void* pk = nullptr;
  void* evk = nullptr;

  // Init
  t.start();
  HE::init(&parameters);
  t.stop("Init");

  // Keygen
  t.start();
  HE::keygen(parameters, &sk, &pk, &evk);
  t.stop("Keygen");
  // HE::serialize_sk("sk.bin", sk);
  // std::cout << "serialized?" << std::endl;

  // Random messages
  unsigned long* messages1 = new unsigned long[nb_tests];
  unsigned long* messages2 = new unsigned long[nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    messages1[i] =
        mpz_class(prng.get_z_range(HE::plaintext_modulus)).get_ui();
    messages2[i] =
        mpz_class(prng.get_z_range(HE::plaintext_modulus)).get_ui();
  }

  // Encrypt
  void** ciphertexts1 = new void* [nb_tests];
  void** ciphertexts2 = new void* [nb_tests];
  t.start();
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::encryptInteger(pk, &(ciphertexts1[i]), messages1[i]);
    HE::encryptInteger(pk, &(ciphertexts2[i]), messages2[i]);
  }
  t.stop("Encrypt Integer", nb_tests * 2);

  // Decrypt
  unsigned long* messages1_decrypted = new unsigned long[nb_tests];
  unsigned long* messages2_decrypted = new unsigned long[nb_tests];
  t.start();
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::decryptInteger(sk, pk, ciphertexts1[i], &(messages1_decrypted[i]));
    HE::decryptInteger(sk, pk, ciphertexts2[i], &(messages2_decrypted[i]));
  }
  t.stop("Decrypt Integer", nb_tests * 2);

  // Correctness of decryption
  for (unsigned long i = 0; i < nb_tests; i++) {
    assert(messages1[i] == messages1_decrypted[i]);
    assert(messages2[i] == messages2_decrypted[i]);
  }

  // Homomorphic additions
  unsigned long* messages_added = new unsigned long[nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    messages_added[i] = (messages1[i] + messages2[i]) % HE::plaintext_modulus;
  }
  void** ciphertexts_added = new void* [nb_tests];
  t.start();
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::add(pk, evk, &(ciphertexts_added[i]), ciphertexts1[i], ciphertexts2[i]);
  }
  t.stop("Homomorphic Addition", nb_tests);

  // Correctness of addition
  unsigned long* messages_added_decrypted = new unsigned long[nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::decryptInteger(sk, pk, ciphertexts_added[i],
                     &(messages_added_decrypted[i]));
    assert(messages_added_decrypted[i] == messages_added[i]);
  }

  // Homomorphic multiplications
  unsigned long* messages_multiplied = new unsigned long[nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    messages_multiplied[i] = (messages1[i] * messages2[i]) % HE::plaintext_modulus;
  }
  void** ciphertexts_multiplied = new void* [nb_tests];
  t.start();
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::mul(pk, evk, &(ciphertexts_multiplied[i]), ciphertexts1[i],
          ciphertexts2[i]);
  }
  t.stop("Homomorphic Multiplication", nb_tests);

  // Correctness of multiplication
  unsigned long* messages_multiplied_decrypted = new unsigned long[nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::decryptInteger(sk, pk, ciphertexts_multiplied[i],
                     &(messages_multiplied_decrypted[i]));
    assert(messages_multiplied_decrypted[i] == messages_multiplied[i]);
  }

  // Homomorphic scalar multiplications
  void** ciphertexts_sc_multiplied = new void* [nb_tests];
  t.start();
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::mulByConstant(pk, evk, &(ciphertexts_sc_multiplied[i]), ciphertexts1[i],
    		messages2[i]);
  }
  t.stop("Homomorphic Scalar Multiplication", nb_tests);

  // Correctness of scalar multiplication
  unsigned long* messages_sc_multiplied_decrypted = new unsigned long[nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::decryptInteger(sk, pk, ciphertexts_sc_multiplied[i],
                     &(messages_sc_multiplied_decrypted[i]));
    assert(messages_sc_multiplied_decrypted[i] == messages_multiplied[i]);
  }

  delete[] messages1;
  delete[] messages2;
  delete[] messages_added;
  delete[] messages_added_decrypted;
  delete[] messages_multiplied;
  delete[] messages_multiplied_decrypted;
  delete[] messages_sc_multiplied_decrypted;

  for(long i=0; i< nb_tests; i++)
  {
	  HE::freeup_ciphertext(pk,ciphertexts1[i]);
	  HE::freeup_ciphertext(pk,ciphertexts2[i]);
	  HE::freeup_ciphertext(pk,ciphertexts_added[i]);
	  HE::freeup_ciphertext(pk,ciphertexts_multiplied[i]);
	  HE::freeup_ciphertext(pk,ciphertexts_sc_multiplied[i]);
  }

  delete[] ciphertexts1;
  delete[] ciphertexts2;
  delete[] ciphertexts_added;
  delete[] ciphertexts_multiplied;
  delete[] ciphertexts_sc_multiplied;

  HE::freeup_keys(parameters,sk,pk,evk);

  return 0;
}

template <unsigned long nb_tests>
int testInteger_polynomial() {
  // Determine if the test should run
    #ifdef BINARY
  unsigned long run = HE::supports_bit_encryption;
  #else
  unsigned long run = HE::supports_unsigned_encryption;
  #endif
  run *= HE::supports_polynomial_encryption;

  if (!run) {
    return 0;
  }

  timing t;

  gmp_randclass prng(gmp_randinit_default);
//  prng.seed(0);
  prng.seed(time(NULL)); // To set different seeds

  void* parameters = nullptr;
  void* sk = nullptr;
  void* pk = nullptr;
  void* evk = nullptr;

  // Init
  t.start();
  HE::init(&parameters);
  t.stop("Init");

  // Keygen
  t.start();
  HE::keygen(parameters, &sk, &pk, &evk);
  t.stop("Keygen");

  // Size of polynomials
  size_t number_of_polynomial_slots =
      HE::number_of_polynomial_slots(pk);

  // Random polynomials
  unsigned long** polynomials1 = new unsigned long* [nb_tests];
  unsigned long** polynomials2 = new unsigned long* [nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    polynomials1[i] = new unsigned long[number_of_polynomial_slots];
    polynomials2[i] = new unsigned long[number_of_polynomial_slots];
    for (unsigned long j = 0; j < number_of_polynomial_slots; j++) {
      polynomials1[i][j] =
          mpz_class(prng.get_z_range(HE::plaintext_modulus)).get_ui();
      polynomials2[i][j] =
          mpz_class(prng.get_z_range(HE::plaintext_modulus)).get_ui();
    }
  }

  // Encrypt
  void** ciphertexts1 = new void* [nb_tests];
  void** ciphertexts2 = new void* [nb_tests];
  t.start();
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::encryptPolynomial(pk, &(ciphertexts1[i]), polynomials1[i],
                        number_of_polynomial_slots);
    HE::encryptPolynomial(pk, &(ciphertexts2[i]), polynomials2[i],
                        number_of_polynomial_slots);
  }
  t.stop("Encrypt Polynomial", nb_tests * 2);

  // Decrypt
  unsigned long** polynomials1_decrypted = new unsigned long* [nb_tests];
  unsigned long** polynomials2_decrypted = new unsigned long* [nb_tests];
  t.start();
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::decryptPolynomial(sk, pk, ciphertexts1[i], &(polynomials1_decrypted[i]),
                        number_of_polynomial_slots);
    HE::decryptPolynomial(sk, pk, ciphertexts2[i], &(polynomials2_decrypted[i]),
                        number_of_polynomial_slots);
  }
  t.stop("Decrypt Polynomial", nb_tests * 2);

  // Correctness of decryption
  for (unsigned long i = 0; i < nb_tests; i++) {
    for (unsigned long j = 0; j < number_of_polynomial_slots; j++) {
      assert(polynomials1[i][j] == polynomials1_decrypted[i][j]);
      assert(polynomials2[i][j] == polynomials2_decrypted[i][j]);
    }
  }

  // Homomorphic additions
  unsigned long** polynomials_added = new unsigned long* [nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    polynomials_added[i] = new unsigned long[number_of_polynomial_slots];
    for (unsigned long j = 0; j < number_of_polynomial_slots; j++) {
      polynomials_added[i][j] =
          (polynomials1[i][j] + polynomials2[i][j]) % HE::plaintext_modulus;
    }
  }
  void** ciphertexts_added = new void* [nb_tests];
  t.start();
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::add(pk, evk, &(ciphertexts_added[i]), ciphertexts1[i], ciphertexts2[i]);
  }
  t.stop("Homomorphic Addition", nb_tests);

  // Correctness of addition
  unsigned long** polynomials_added_decrypted = new unsigned long* [nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    polynomials_added_decrypted[i] =
        new unsigned long[number_of_polynomial_slots];
    HE::decryptPolynomial(sk, pk, ciphertexts_added[i],
                        &(polynomials_added_decrypted[i]),
                        number_of_polynomial_slots);
    for (unsigned long j = 0; j < number_of_polynomial_slots; j++) {
      assert(polynomials_added_decrypted[i][j] == polynomials_added[i][j]);
    }
  }

  // Homomorphic multiplications
  void** ciphertexts_multiplied = new void* [nb_tests];
  t.start();
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::mul(pk, evk, &(ciphertexts_multiplied[i]), ciphertexts1[i],
          ciphertexts2[i]);
  }
  t.stop("Homomorphic Multiplication", nb_tests);

  delete[] polynomials1;
  delete[] polynomials2;
  delete[] polynomials_added;
  delete[] polynomials_added_decrypted;

  for(long i=0; i< nb_tests; i++)
  {
	  HE::freeup_ciphertext(pk,ciphertexts1[i]);
	  HE::freeup_ciphertext(pk,ciphertexts2[i]);
	  HE::freeup_ciphertext(pk,ciphertexts_added[i]);
	  HE::freeup_ciphertext(pk,ciphertexts_multiplied[i]);
  }

  delete[] ciphertexts1;
  delete[] ciphertexts2;
  delete[] ciphertexts_added;
  delete[] ciphertexts_multiplied;

  HE::freeup_keys(parameters,sk,pk,evk);

  return 0;
}

template <unsigned long nb_tests>
int testInteger_vector() {
  // Determine if the test should run
  #ifdef BINARY
  unsigned long run = HE::supports_bit_encryption;
  #else
  unsigned long run = HE::supports_unsigned_encryption;
  #endif
  run *= HE::supports_vector_encryption;

  if (!run) {
    return 0;
  }

  timing t;

  gmp_randclass prng(gmp_randinit_default);
//  prng.seed(0);
  prng.seed(time(NULL)); // To set different seeds

  void* parameters = nullptr;
  void* sk = nullptr;
  void* pk = nullptr;
  void* evk = nullptr;

  // Init
  t.start();
  HE::init(&parameters);
  t.stop("Init");

  // Keygen
  t.start();
  HE::keygen(parameters, &sk, &pk, &evk);
  t.stop("Keygen");

  // Number of slots
  size_t number_of_vector_slots = HE::number_of_vector_slots(pk);

  // Random polynomials
  unsigned long** vectors1 = new unsigned long* [nb_tests];
  unsigned long** vectors2 = new unsigned long* [nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    vectors1[i] = new unsigned long[number_of_vector_slots];
    vectors2[i] = new unsigned long[number_of_vector_slots];
    for (unsigned long j = 0; j < number_of_vector_slots; j++) {
      vectors1[i][j] =
          mpz_class(prng.get_z_range(HE::plaintext_modulus)).get_ui();
      vectors2[i][j] =
          mpz_class(prng.get_z_range(HE::plaintext_modulus)).get_ui();
    }
  }

  // Encrypt
  void** ciphertexts1 = new void* [nb_tests];
  void** ciphertexts2 = new void* [nb_tests];
  t.start();
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::encryptVector(pk, &(ciphertexts1[i]), vectors1[i],
                    number_of_vector_slots);
    HE::encryptVector(pk, &(ciphertexts2[i]), vectors2[i],
                    number_of_vector_slots);
  }
  t.stop("Encrypt Vector", nb_tests * 2);

  // Decrypt
  unsigned long** vectors1_decrypted = new unsigned long* [nb_tests];
  unsigned long** vectors2_decrypted = new unsigned long* [nb_tests];
  t.start();
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::decryptVector(sk, pk, ciphertexts1[i], &(vectors1_decrypted[i]),
                    number_of_vector_slots);
    HE::decryptVector(sk, pk, ciphertexts2[i], &(vectors2_decrypted[i]),
                    number_of_vector_slots);
  }
  t.stop("Decrypt Vector", nb_tests * 2);

  // Correctness of decryption
  for (unsigned long i = 0; i < nb_tests; i++) {
    for (unsigned long j = 0; j < number_of_vector_slots; j++) {
      assert(vectors1[i][j] == vectors1_decrypted[i][j]);
      assert(vectors2[i][j] == vectors2_decrypted[i][j]);
    }
  }

  // Homomorphic additions
  unsigned long** vectors_added = new unsigned long* [nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    vectors_added[i] = new unsigned long[number_of_vector_slots];
    for (unsigned long j = 0; j < number_of_vector_slots; j++) {
      vectors_added[i][j] =
          (vectors1[i][j] + vectors2[i][j]) % HE::plaintext_modulus;
    }
  }
  void** ciphertexts_added = new void* [nb_tests];
  t.start();
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::add(pk, evk, &(ciphertexts_added[i]), ciphertexts1[i], ciphertexts2[i]);
  }
  t.stop("Homomorphic Addition", nb_tests);

  // Correctness of addition
  unsigned long** vectors_added_decrypted = new unsigned long* [nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::decryptVector(sk, pk, ciphertexts_added[i], &(vectors_added_decrypted[i]),
                    number_of_vector_slots);
    for (unsigned long j = 0; j < number_of_vector_slots; j++) {
      assert(vectors_added_decrypted[i][j] == vectors_added[i][j]);
    }
  }

  // Homomorphic multiplications
  unsigned long** vectors_multiplied = new unsigned long* [nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    vectors_multiplied[i] = new unsigned long[number_of_vector_slots];
    for (unsigned long j = 0; j < number_of_vector_slots; j++) {
      vectors_multiplied[i][j] =
          (vectors1[i][j] * vectors2[i][j]) % HE::plaintext_modulus;
    }
  }
  void** ciphertexts_multiplied = new void* [nb_tests];
  t.start();
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::mul(pk, evk, &(ciphertexts_multiplied[i]), ciphertexts1[i],
          ciphertexts2[i]);
  }
  t.stop("Homomorphic Multiplication", nb_tests);

  // Correctness of multiplication
  unsigned long** vectors_multiplied_decrypted = new unsigned long* [nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::decryptVector(sk, pk, ciphertexts_multiplied[i],
                    &(vectors_multiplied_decrypted[i]), number_of_vector_slots);
    for (unsigned long j = 0; j < number_of_vector_slots; j++) {
      assert(vectors_multiplied_decrypted[i][j] == vectors_multiplied[i][j]);
    }
  }

  delete[] vectors1;
  delete[] vectors2;
  delete[] vectors_added;
  delete[] vectors_added_decrypted;
  delete[] vectors_multiplied;
  delete[] vectors_multiplied_decrypted;

  for(long i=0; i< nb_tests; i++)
  {
	  HE::freeup_ciphertext(pk,ciphertexts1[i]);
	  HE::freeup_ciphertext(pk,ciphertexts2[i]);
	  HE::freeup_ciphertext(pk,ciphertexts_added[i]);
	  HE::freeup_ciphertext(pk,ciphertexts_multiplied[i]);
  }

  delete[] ciphertexts1;
  delete[] ciphertexts2;
  delete[] ciphertexts_added;
  delete[] ciphertexts_multiplied;

  HE::freeup_keys(parameters,sk,pk,evk);

  return 0;
}

template <unsigned long nb_tests>
int testInteger_mpz() {
  // Determine if the test should run
  #ifdef MPZ
  unsigned long run = HE::supports_mpz_encryption;
  #else
  unsigned long run = 0;
  #endif

  if (!run) {
    return 0;
  }

  timing t;

  gmp_randclass prng(gmp_randinit_default);
//  prng.seed(0);
  prng.seed(time(NULL)); // To set different seeds

  void* parameters = nullptr;
  void* sk = nullptr;
  void* pk = nullptr;
  void* evk = nullptr;

  // Init
  t.start();
  HE::init(&parameters);
  t.stop("Init");

  // Keygen
  t.start();
  HE::keygen(parameters, &sk, &pk, &evk);
  t.stop("Keygen");

  // Random messages
  mpz_class* messages1 = new mpz_class[nb_tests];
  mpz_class* messages2 = new mpz_class[nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    messages1[i] = prng.get_z_range(HE::plaintext_modulus_mpz);
    messages2[i] = prng.get_z_range(HE::plaintext_modulus_mpz);
  }

  // Encrypt
  void** ciphertexts1 = new void* [nb_tests];
  void** ciphertexts2 = new void* [nb_tests];
  t.start();
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::encryptInteger(pk, &(ciphertexts1[i]), messages1[i]);
    HE::encryptInteger(pk, &(ciphertexts2[i]), messages2[i]);
  }
  t.stop("Encrypt Integer", nb_tests * 2);

  // Decrypt
  mpz_class* messages1_decrypted = new mpz_class[nb_tests];
  mpz_class* messages2_decrypted = new mpz_class[nb_tests];
  t.start();
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::decryptInteger(sk, pk, ciphertexts1[i], &(messages1_decrypted[i]));
    HE::decryptInteger(sk, pk, ciphertexts2[i], &(messages2_decrypted[i]));
  }
  t.stop("Decrypt Integer", nb_tests * 2);

  // Correctness of decryption
  for (unsigned long i = 0; i < nb_tests; i++) {
    assert(messages1[i] == messages1_decrypted[i]);
    assert(messages2[i] == messages2_decrypted[i]);
  }

  // Homomorphic additions
  mpz_class* messages_added = new mpz_class[nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    messages_added[i] = (messages1[i] + messages2[i]) % HE::plaintext_modulus_mpz;
  }
  void** ciphertexts_added = new void* [nb_tests];
  t.start();
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::add(pk, evk, &(ciphertexts_added[i]), ciphertexts1[i], ciphertexts2[i]);
  }
  t.stop("Homomorphic Addition", nb_tests);

  // Correctness of addition
  mpz_class* messages_added_decrypted = new mpz_class[nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::decryptInteger(sk, pk, ciphertexts_added[i],
                     &(messages_added_decrypted[i]));
    assert(messages_added_decrypted[i] == messages_added[i]);
  }

  // Homomorphic multiplications
  mpz_class* messages_multiplied = new mpz_class[nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    messages_multiplied[i] = (messages1[i] * messages2[i]) % HE::plaintext_modulus_mpz;
  }
  void** ciphertexts_multiplied = new void* [nb_tests];
  t.start();
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::mul(pk, evk, &(ciphertexts_multiplied[i]), ciphertexts1[i],
          ciphertexts2[i]);
  }
  t.stop("Homomorphic Multiplication", nb_tests);

  // Correctness of multiplication
  mpz_class* messages_multiplied_decrypted = new mpz_class[nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::decryptInteger(sk, pk, ciphertexts_multiplied[i],
                     &(messages_multiplied_decrypted[i]));
    assert(messages_multiplied_decrypted[i] == messages_multiplied[i]);
  }

  // Homomorphic scalar multiplications
  void** ciphertexts_sc_multiplied = new void* [nb_tests];
  t.start();
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::mulByConstant(pk, evk, &(ciphertexts_sc_multiplied[i]), ciphertexts1[i],
    		messages2[i]);
  }
  t.stop("Homomorphic Scalar Multiplication", nb_tests);

  // Correctness of scalar multiplication
  mpz_class* messages_sc_multiplied_decrypted = new mpz_class[nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::decryptInteger(sk, pk, ciphertexts_sc_multiplied[i],
                     &(messages_sc_multiplied_decrypted[i]));
    assert(messages_sc_multiplied_decrypted[i] == messages_multiplied[i]);
  }

  delete[] messages1;
  delete[] messages2;
  delete[] messages_added;
  delete[] messages_added_decrypted;
  delete[] messages_multiplied;
  delete[] messages_multiplied_decrypted;
  delete[] messages_sc_multiplied_decrypted;

  for(unsigned long i=0; i< nb_tests; i++)
  {
	  HE::freeup_ciphertext(pk,ciphertexts1[i]);
	  HE::freeup_ciphertext(pk,ciphertexts2[i]);
	  HE::freeup_ciphertext(pk,ciphertexts_added[i]);
	  HE::freeup_ciphertext(pk,ciphertexts_multiplied[i]);
	  HE::freeup_ciphertext(pk,ciphertexts_sc_multiplied[i]);
  }

  delete[] ciphertexts1;
  delete[] ciphertexts2;
  delete[] ciphertexts_added;
  delete[] ciphertexts_multiplied;
  delete[] ciphertexts_sc_multiplied;

  HE::freeup_keys(parameters,sk,pk,evk);

  return 0;
}

#ifdef FXPT

#define MAX_PRECISION_DIFF 3 // Make sure that this value is suitably chosen depending upon the values of MAX_PRECISION in params.h files and the input interval.

#define checkPreciMacro(diff,preci, PreciInt, out, in)\
{\
	if( (PreciInt) < 0)\
	{\
		if( (preci) >= (long) floor(log2((long) fabs(in))+1) )\
			diff =  (long)(fabs((out)-(in)) * (1L << ((preci)- (long) floor(log2((long) fabs(in))+1) )) );\
		else\
			diff =  ((long)(fabs((out)-(in)))) >> (-1* ((preci)- (long) floor(log2((long) fabs(in))+1) )) ;\
	}\
	else\
	{\
		if( (preci) >= (PreciInt))\
				diff =  (long)(fabs((out)-(in)) * (1L << ((preci)-(PreciInt))));\
		else\
			diff =  ((long)(fabs((out)-(in)))) >> (-1* ((preci)-(PreciInt) )) ;\
	}\
}


template <unsigned long nb_tests>
int testFixedpt() {

	unsigned long run = HE::supports_fixedpt_encryption;

	if (!run) {
    return 0;
  }

  timing t;

  gmp_randclass prng(gmp_randinit_default);
//  prng.seed(0);
  prng.seed(time(NULL)); // To set different seeds

  void* parameters = nullptr;
  void* sk = nullptr;
  void* pk = nullptr;
  void* evk = nullptr;

  // Init
  t.start();
  HE::init(&parameters);
  t.stop("Init");

  // Keygen
  t.start();
  HE::keygen(parameters, &sk, &pk, &evk);
  t.stop("Keygen");

  // Random messages
  double* messages1 = new double[nb_tests];
  double* messages2 = new double[nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    messages1[i] =
        mpf_class(prng.get_f(sizeof(double))).get_d();
    if(messages1[i] > 0.5)
    	messages1[i] -= 1;
    messages1[i] *= 16;		//Scaling by an arbitrary constant.
    messages2[i] =
        mpf_class(prng.get_f(sizeof(double))).get_d();
    if(messages2[i] > 0.5)
    	messages2[i] -= 1;
    messages2[i] *= 16;		//Scaling by an arbitrary constant.
  }

  // Encrypt
  void** ciphertexts1 = new void* [nb_tests];
  void** ciphertexts2 = new void* [nb_tests];
  t.start();
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::encryptFixedpt(pk, &(ciphertexts1[i]), messages1[i]);
    HE::encryptFixedpt(pk, &(ciphertexts2[i]), messages2[i]);
  }
  t.stop("Encrypt fixed-point number", nb_tests * 2);

  // Decrypt
  double* messages1_decrypted = new double[nb_tests];
  double* messages2_decrypted = new double[nb_tests];
  t.start();
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::decryptFixedpt(sk, pk, ciphertexts1[i], &(messages1_decrypted[i]));
    HE::decryptFixedpt(sk, pk, ciphertexts2[i], &(messages2_decrypted[i]));
  }
  t.stop("Decrypt fixed-point number", nb_tests * 2);

  // Correctness of decryption
  for (unsigned long i = 0; i < nb_tests; i++) {
	long diff;
	checkPreciMacro(diff, MAX_PRECISION_DIFF, -1, messages1_decrypted[i], messages1[i]);
    assert(diff == 0);
    checkPreciMacro(diff, MAX_PRECISION_DIFF, -1, messages2_decrypted[i], messages2[i]);
    assert(diff == 0);
  }

  // Homomorphic additions
  double* messages_added = new double[nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    messages_added[i] = (messages1[i] + messages2[i]);
  }
  void** ciphertexts_added = new void* [nb_tests];
  t.start();
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::add(pk, evk, &(ciphertexts_added[i]), ciphertexts1[i], ciphertexts2[i]);
  }
  t.stop("Homomorphic Addition", nb_tests);

  // Correctness of addition
  double* messages_added_decrypted = new double[nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::decryptFixedpt(sk, pk, ciphertexts_added[i],
                     &(messages_added_decrypted[i]));
	long diff;
	checkPreciMacro(diff, MAX_PRECISION_DIFF, -1, messages_added_decrypted[i], messages_added[i]);
    assert(diff == 0);
  }

  // Homomorphic multiplications
  double* messages_multiplied = new double[nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    messages_multiplied[i] = (messages1[i] * messages2[i]);
  }
  void** ciphertexts_multiplied = new void* [nb_tests];
  t.start();
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::mul(pk, evk, &(ciphertexts_multiplied[i]), ciphertexts1[i],
          ciphertexts2[i]);
  }
  t.stop("Homomorphic Multiplication", nb_tests);

  // Correctness of multiplication
  double* messages_multiplied_decrypted = new double[nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::decryptFixedpt(sk, pk, ciphertexts_multiplied[i],
                     &(messages_multiplied_decrypted[i]));
	long diff;
	checkPreciMacro(diff, MAX_PRECISION_DIFF, -1, messages_multiplied_decrypted[i], messages_multiplied[i]);
    assert(diff == 0);
  }


  // Homomorphic scalar multiplications
  void** ciphertexts_sc_multiplied = new void* [nb_tests];
  t.start();
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::mulByConstant(pk, evk, &(ciphertexts_sc_multiplied[i]), ciphertexts1[i],
    		messages2[i]);
  }
  t.stop("Homomorphic Scalar Multiplication", nb_tests);

  // Correctness of scalar multiplication
  double* messages_sc_multiplied_decrypted = new double[nb_tests];
  for (unsigned long i = 0; i < nb_tests; i++) {
    HE::decryptFixedpt(sk, pk, ciphertexts_sc_multiplied[i],
                     &(messages_sc_multiplied_decrypted[i]));
	long diff;
	checkPreciMacro(diff, MAX_PRECISION_DIFF, -1, messages_sc_multiplied_decrypted[i], messages_multiplied[i]);
    assert(diff == 0);
  }


  delete[] messages1;
  delete[] messages2;
  delete[] messages_added;
  delete[] messages_added_decrypted;
  delete[] messages_multiplied;
  delete[] messages_multiplied_decrypted;
  delete[] messages_sc_multiplied_decrypted;

  for(unsigned long i=0; i< nb_tests; i++)
  {
	  HE::freeup_ciphertext(pk,ciphertexts1[i]);
	  HE::freeup_ciphertext(pk,ciphertexts2[i]);
	  HE::freeup_ciphertext(pk,ciphertexts_added[i]);
	  HE::freeup_ciphertext(pk,ciphertexts_multiplied[i]);
	  HE::freeup_ciphertext(pk,ciphertexts_sc_multiplied[i]);
  }

  delete[] ciphertexts1;
  delete[] ciphertexts2;
  delete[] ciphertexts_added;
  delete[] ciphertexts_multiplied;
  delete[] ciphertexts_sc_multiplied;

  HE::freeup_keys(parameters,sk,pk,evk);

  return 0;
}

#endif // FXPT

#endif // HE_TESTS

