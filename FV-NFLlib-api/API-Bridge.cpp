#include <iostream>
#include <cassert>
#include <vector>

#include <fstream>
#include <cereal/archives/binary.hpp>

#include "API-Bridge.h"

using namespace std;
using namespace FV;

namespace HE {

mpz_class plaintext_modulus_mpz(PLAINTEXT_MODULUS_MPZ);
unsigned long plaintext_modulus(PLAINTEXT_MODULUS);
unsigned long security_level(SECURITY_LEVEL);
unsigned long multiplicative_depth(MULTIPLICATIVE_DEPTH);

bool supports_bit_encryption = true;
bool supports_unsigned_encryption = true;
bool supports_mpz_encryption = true;
bool supports_polynomial_encryption = true;
bool supports_vector_encryption = false;

size_t number_of_polynomial_slots(void* pk) { return params::poly_p::degree; }
size_t number_of_vector_slots(void* pk) { return 0; }

/**
 * Init function
 */
int init(void** parameters) {
  return 0;
}

/**
 * Key Generation
 */
int keygen(void* parameters, void** sk, void** pk, void** evk) {
  *sk = (void*)new sk_t();
  *evk = (void*)new evk_t(*((sk_t*)*sk), 32);
  *pk = (void*)new pk_t(*((sk_t*)*sk), *((evk_t*)*evk));
  return 0;
}

/**
 * Encrypt in integer (usigned or mpz_class) using the public key
 */
int encryptInteger(void* pk, void** ciphertext, unsigned long message,
                      unsigned long level) {
  encryptInteger(pk, ciphertext, mpz_class(message), level);
  return 0;
}
int encryptInteger(void* pk, void** ciphertext, mpz_class const& message,
                      unsigned long level) {
  pk_t* PK = (pk_t*)pk;
  mess_t m(message);
  *ciphertext = (void*)new ciphertext_t();
//  encrypt<pk_t,ciphertext_t,mess_t>(*((ciphertext_t*)*ciphertext), *PK, m);
  encrypt(*((ciphertext_t*)*ciphertext), *PK, m);
  return 0;
}

/**
 * Encrypt in integer (usigned or mpz_class) using the secret key
 */
int encryptIntegerWithSK(void* sk, void* pk, void** ciphertext,
                            unsigned long message, unsigned long level) {
  encryptInteger(pk, ciphertext, message, level);
  return 0;
}
int encryptIntegerWithSK(void* sk, void* pk, void** ciphertext,
                            mpz_class const& message, unsigned long level) {
  encryptInteger(pk, ciphertext, message, level);
  return 0;
}

/**
 * Encrypt a vector (of usigned or mpz_class) using the public key [BATCHING]
 */
int encryptVector(void* pk, void** ciphertext, unsigned long* message_v,
                     size_t size, unsigned long level) {
  assert(false);
  return 1;
}
int encryptVector(void* pk, void** ciphertext, mpz_class* message_v,
                     size_t size, unsigned long level) {
  assert(false);
  return 1;
}

/**
 * Encrypt a vector (of usigned or mpz_class) using the secret key [BATCHING]
 */
int encryptVectorWithSK(void* sk, void* pk, void** ciphertext,
                           unsigned long* message_v, size_t size,
                           unsigned long level) {
  assert(false);
  return 1;
}
int encryptVectorWithSK(void* sk, void* pk, void** ciphertext,
                           mpz_class* message_v, size_t size,
                           unsigned long level) {
  assert(false);
  return 1;
}

/**
 * Encrypt a polynomial (of usigned or mpz_class) using the public key
 */
int encryptPolynomial(void* pk, void** ciphertext, unsigned long* message_p,
                         size_t size, unsigned long level) {
  assert(size == number_of_polynomial_slots(pk));
  pk_t* PK = (pk_t*)pk;

  std::vector<unsigned long> m_v;
  for (size_t i = 0; i < size; i++) {
    m_v.push_back(message_p[i]);
  }

  params::poly_p* m = (params::poly_p*) malloc(sizeof(params::poly_p));
  m->set(m_v.begin(), m_v.end());

  *ciphertext = (void*)new ciphertext_t();
  encrypt_poly(*((ciphertext_t*)*ciphertext), *PK, *m);

  free(m);
  return 0;
}
int encryptPolynomial(void* pk, void** ciphertext, mpz_class* message_p,
                         size_t size, unsigned long level) {
  assert(size == number_of_polynomial_slots(pk));
  pk_t* PK = (pk_t*)pk;

  std::vector<mpz_class> m_v;
  for (size_t i = 0; i < size; i++) {
    m_v.push_back(message_p[i]);
  }

  params::poly_p* m = (params::poly_p*) malloc(sizeof(params::poly_p));
  m->set(m_v.begin(), m_v.end());

  *ciphertext = (void*)new ciphertext_t();
  encrypt_poly(*((ciphertext_t*)*ciphertext), *PK, *m);

  free(m);
  return 0;
}

/**
 * Encrypt a polynomial (of usigned or mpz_class) using the secret key
 */
int encryptPolynomialWithSK(void* sk, void* pk, void** ciphertext,
                               unsigned long* message_p, size_t size,
                               unsigned long level) {
  encryptPolynomial(pk, ciphertext, message_p, size, level);
  return 0;
}
int encryptPolynomialWithSK(void* sk, void* pk, void** ciphertext,
                               mpz_class* message_p, size_t size,
                               unsigned long level) {
  encryptPolynomial(pk, ciphertext, message_p, size, level);
  return 0;
}

/**
 * Decrypt a ciphertext to an integer (unsigned long or mpz_class)
 */
int decryptInteger(void* sk, void* pk, void* ciphertext,
                      unsigned long* message, unsigned long level) {
  pk_t* PK = (pk_t*)pk;
  sk_t* SK = (sk_t*)sk;
  ciphertext_t* c = (ciphertext_t*)ciphertext;

  mess_t m;
  decrypt<sk_t,pk_t,ciphertext_t,mess_t>(m, *SK, *PK, *c);
  *message = m.getValue().get_ui();
  return 0;
}
int decryptInteger(void* sk, void* pk, void* ciphertext, mpz_class* message,
                      unsigned long level) {
  pk_t* PK = (pk_t*)pk;
  sk_t* SK = (sk_t*)sk;
  ciphertext_t* c = (ciphertext_t*)ciphertext;

  mess_t m;
  decrypt<sk_t,pk_t,ciphertext_t,mess_t>(m, *SK, *PK, *c);
  *message = m.getValue();
  return 0;
}

/**
 * Decrypt a ciphertext to a polynomial (of unsigned long of mpz_class)
 */
int decryptPolynomial(void* sk, void* pk, void* ciphertext,
                         unsigned long** message_p, size_t size,
                         unsigned long level) {
  assert(size <= number_of_polynomial_slots(pk));
  
  pk_t* PK = (pk_t*)pk;
  sk_t* SK = (sk_t*)sk;
  ciphertext_t* c = (ciphertext_t*)ciphertext;
  
  std::array<mpz_t, params::poly_p::degree> polym;
  for (size_t i = 0; i < params::poly_p::degree; i++) {
    mpz_init(polym[i]);
  }

  decrypt_poly(polym, *SK, *PK, *c);

  *message_p = new unsigned long[size];
  for (size_t i=0; i<size; i++) {
    (*message_p)[i] = mpz_get_ui(polym[i]);
  }

  for (size_t i = 0; i < params::poly_p::degree; i++) {
    mpz_clear(polym[i]);
  }
  return 0;
}
int decryptPolynomial(void* sk, void* pk, void* ciphertext,
                         mpz_class** message_p, size_t size,
                         unsigned long level) {
  assert(size <= number_of_polynomial_slots(pk));
  
  pk_t* PK = (pk_t*)pk;
  sk_t* SK = (sk_t*)sk;
  ciphertext_t* c = (ciphertext_t*)ciphertext;
  
  std::array<mpz_t, params::poly_p::degree> polym;
  for (size_t i = 0; i < params::poly_p::degree; i++) {
    mpz_init(polym[i]);
  }

  decrypt_poly(polym, *SK, *PK, *c);

  *message_p = new mpz_class[size];
  for (size_t i=0; i<size; i++) {
    (*message_p)[i] = mpz_class(polym[i]);
  }

  for (size_t i = 0; i < params::poly_p::degree; i++) {
    mpz_clear(polym[i]);
  }
  return 0;
}

/**
 * Decrypt a ciphertext to a vector (of unsigned long of mpz_class) [BATCHING]
 */
int decryptVector(void* sk, void* pk, void* ciphertext,
                     unsigned long** message_v, size_t size,
                     unsigned long level) {
  assert(false);
  return 1;
}
int decryptVector(void* sk, void* pk, void* ciphertext,
                     mpz_class** message_v, size_t size, unsigned long level) {
  assert(false);
  return 1;
}

/**
 * Homomorphic addition of two ciphertexts
 */
int add(void* pk, void* evk, void** output, void* input1, void* input2,
           unsigned long level1, unsigned long level2) {
  ciphertext_t* c1 = (ciphertext_t*)input1;
  ciphertext_t* c2 = (ciphertext_t*)input2;
  ciphertext_t* c_add = new ciphertext_t((*c1) + (*c2));
  *output = (void*)c_add;
  return 0;
}

/**
 * Homomorphic multiplication of two ciphertexts
 */
int mul(void* pk, void* evk, void** output, void* input1, void* input2,
           unsigned long level1, unsigned long level2) {
  ciphertext_t* c1 = (ciphertext_t*)input1;
  ciphertext_t* c2 = (ciphertext_t*)input2;
  ciphertext_t* c_add = new ciphertext_t((*c1) * (*c2));
  *output = (void*)c_add;
  return 0;
}

/**
 * Homomorphic multiplication by a constant of a ciphertext
 */
int mulByConstant(void* pk, void* evk, void** output, void* input,
                     unsigned long constant, unsigned long level) {

  ciphertext_t* c_in = (ciphertext_t*)input;
  mess_t multiplier(constant);
  ciphertext_t* c_out = new ciphertext_t((*c_in) * multiplier);
  *output = (void*)c_out;
  return 0;
}
int mulByConstant(void* pk, void* evk, void** output, void* input,
                     mpz_class const& constant, unsigned long level) {

  ciphertext_t* c_in = (ciphertext_t*)input;
  mess_t multiplier(constant);
  ciphertext_t* c_out = new ciphertext_t((*c_in) * multiplier);
  *output = (void*)c_out;
  return 0;
}

/**
 * Deallocate the memory occupied by a ciphertext
 */

int freeup_ciphertext(void* pk, void* ciphertext)
{

	delete (static_cast<ciphertext_t *>(ciphertext));
	return 0;
}

/**
 * Deallocate the memory occupied by parameters generated during HE::init() and HE::keygen()
 */

int freeup_keys(void* parameters, void* sk, void* pk, void* evk)
{
	delete (static_cast<pk_t *>(pk));
	delete (static_cast<evk_t *>(evk));
	delete (static_cast<sk_t *>(sk));

	return 0;
}

/**
 * Serialize functions (back and forth)
 */
// int serialize_parameters  (char* filename, void* parameters);
//int serialize_sk          (const char* filename, void* sk) {
//  sk_t *SK = (sk_t *)sk;
//  std::ofstream file(filename, std::ios::binary);
//  cereal::BinaryOutputArchive oarchive(file);
//  oarchive (*SK);
//  return 0;
//}
// int serialize_evk         (char* filename, void* evk);
// int serialize_pk          (char* filename, void* pk);
// int serialize_ciphertext  (char* filename, void* ciphertext);

// int deserialize_parameters(char* filename, void** parameters);
// int deserialize_sk        (char* filename, void** sk);
// int deserialize_evk       (char* filename, void** evk);
// int deserialize_pk        (char* filename, void** pk);
// int deserialize_ciphertext(char* filename, void** ciphertext);

}
