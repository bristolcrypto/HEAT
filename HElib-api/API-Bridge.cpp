/* Note: the definitions of functions for fixed-point arithmetic are conditionally included.*/


#include <iostream>
#include <cassert>
#include <vector>
#include <fstream>

#include "API-Bridge.h"

namespace HE {

unsigned long plaintext_modulus(PLAINTEXT_MODULUS);
mpz_class plaintext_modulus_mpz(PLAINTEXT_MODULUS_MPZ);

#ifdef FXPT
bool supports_bit_encryption = false;
bool supports_fixedpt_encryption = true;
#else
bool supports_bit_encryption = true;
#endif //FXPT

bool supports_unsigned_encryption = true;
bool supports_mpz_encryption = false;
bool supports_polynomial_encryption = true;
bool supports_vector_encryption = false;

unsigned long security_level(SECURITY_LEVEL);
unsigned long multiplicative_depth(MULTIPLICATIVE_DEPTH);

size_t number_of_polynomial_slots(void* pk) { return (size_t) ((static_cast<const FHEPubKey *>(pk))->getContext().zMStar.getPhiM()); }
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

	  long p = plaintext_modulus;
	  long r = 1;
	  long d = 1;
	  long w = 64;
	  long c = 2;
	  long k = security_level;
	  long L= multiplicative_depth;
	  long s = 0;
	  long m = 0;
	  Vec<long> mvec;
	  vector<long> *gens = new vector<long>;
	  vector<long> *ords = new vector<long>;
	  if (mvec.length()>0)
		  m = computeProd(mvec);

	  #ifdef FXPT
	  	  bool fracrep = FRAC_REP_ON;
	  	  long p_fin, m_fin, d_fin;
	  	  d = MIN_RING_DEG;
	  	  FindPM(p_fin, m_fin, d_fin, k, L, c, p, d, s, m, fracrep, true, false);
		  p=p_fin; m=m_fin;
		  plaintext_modulus = p;
	  #else
		  long m_fin = FindM(k, L, c, p, d, s, m, false);
		  m=m_fin;
	  #endif	//FXPT

	  vector<long> *gens1 = new vector<long>;
	  vector<long> *ords1 = new vector<long>;
	  convert(*gens1, *gens);
	  convert(*ords1, *ords);

	  FHEcontext *context = new FHEcontext(m, p, r, *gens1, *ords1);
	  buildModChain(*context, L, c);

	  context->zMStar.printout();

	  FHESecKey* secretKey = new FHESecKey(*context);
	  secretKey->GenSecKey(w); // A Hamming-weight-w secret key

	  const FHEPubKey* publicKey = new const FHEPubKey(*secretKey);

	  *sk = (void *) secretKey;
	  *pk = (void *) publicKey;

	  addSome1DMatrices(*secretKey); // compute key-switching matrices that we need

	  return 0;
}

#ifndef FXPT

/**
 * Encrypt in integer (usigned or mpz_class) using the public key
 */
int encryptInteger(void* pk, void** ciphertext, unsigned long message,
                      unsigned long level) {
	ZZX ply((long)message);
	Ctxt *ct = new Ctxt(*(static_cast<const FHEPubKey *>(pk)));
	(static_cast<const FHEPubKey *>(pk))->Encrypt(*ct,ply);
	*ciphertext = (void *) ct;

	return 1;
}

int encryptInteger(void* pk, void** ciphertext, mpz_class const& message,
                      unsigned long level) {
	  assert(false);
	  return 1;
}

/**
 * Encrypt in integer (usigned or mpz_class) using the secret key
 */
int encryptIntegerWithSK(void* sk, void* pk, void** ciphertext,
                            unsigned long message, unsigned long level) {
	ZZX ply((long)message);
	Ctxt *ct = new Ctxt(*(static_cast<FHESecKey *>(sk)));
	(static_cast<FHESecKey *>(sk))->Encrypt(*ct,ply);
	*ciphertext = (void *) ct;

	return 1;
}

int encryptIntegerWithSK(void* sk, void* pk, void** ciphertext,
                            mpz_class const& message, unsigned long level) {
	  assert(false);
	  return 1;
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

	assert(size <= number_of_polynomial_slots(pk));

	ZZX ply;
	for(size_t i=0; i<size; i++)
		SetCoeff(ply,i,(long) message_p[i]);
	Ctxt *ct = new Ctxt(*(static_cast<const FHEPubKey *>(pk)));
	(static_cast<const FHEPubKey *>(pk))->Encrypt(*ct,ply);
	*ciphertext = (void *) ct;

	return 1;
}

int encryptPolynomial(void* pk, void** ciphertext, mpz_class* message_p,
                         size_t size, unsigned long level) {
	  assert(false);
	  return 1;
}

/**
 * Encrypt a polynomial (of usigned or mpz_class) using the secret key
 */
int encryptPolynomialWithSK(void* sk, void* pk, void** ciphertext,
                               unsigned long* message_p, size_t size,
                               unsigned long level) {

	assert(size <= number_of_polynomial_slots(pk));

	ZZX ply;
	for(size_t i=0; i<size; i++)
		SetCoeff(ply,i,(long) message_p[i]);
	Ctxt *ct = new Ctxt(*(static_cast<FHESecKey *>(sk)));
	(static_cast<FHESecKey *>(sk))->Encrypt(*ct,ply);
	*ciphertext = (void *) ct;

    return 1;
}

int encryptPolynomialWithSK(void* sk, void* pk, void** ciphertext,
                               mpz_class* message_p, size_t size,
                               unsigned long level) {
	  assert(false);
	  return 1;
}


/**
 * Decrypt a ciphertext to an integer (unsigned long of mpz_class)
 */
int decryptInteger(void* sk, void* pk, void* ciphertext,
                      unsigned long* message, unsigned long level) {
	ZZX ply;
    (static_cast<FHESecKey *>(sk))->Decrypt(ply, *(static_cast<Ctxt *>(ciphertext)));
    *message = conv<unsigned long>(ConstTerm(ply));

	return 0;
}

int decryptInteger(void* sk, void* pk, void* ciphertext, mpz_class* message,
                      unsigned long level) {
	  assert(false);
	  return 1;
}


/**
 * Decrypt a ciphertext to a polynomial (of unsigned long of mpz_class)
 */
int decryptPolynomial(void* sk, void* pk, void* ciphertext,
                         unsigned long** message_p, size_t size,
                         unsigned long level) {
  assert(size <= number_of_polynomial_slots(pk));
  ZZX ply;
  (static_cast<FHESecKey *>(sk))->Decrypt(ply, *(static_cast<Ctxt *>(ciphertext)));
  assert((size_t)deg(ply) <= size);

  *message_p = new unsigned long[size];
  for (size_t i=0; i<= (size_t)size; i++) {
    (*message_p)[i] = conv<unsigned long>(coeff(ply,i));
  }

  return 0;
}

int decryptPolynomial(void* sk, void* pk, void* ciphertext,
                         mpz_class** message_p, size_t size,
                         unsigned long level) {
	  assert(false);
	  return 1;
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
  Ctxt* ct1 = static_cast<Ctxt *>(input1);
  Ctxt* ct2 = static_cast<Ctxt *>(input2);
  Ctxt *ct = new Ctxt(*ct1);
  *ct += *ct2;
  *output = (void*) ct;

  return 0;
}

/**
 * Homomorphic multiplication of two ciphertexts
 */
int mul(void* pk, void* evk, void** output, void* input1, void* input2,
           unsigned long level1, unsigned long level2) {
	  Ctxt* ct1 = static_cast<Ctxt *>(input1);
	  Ctxt* ct2 = static_cast<Ctxt *>(input2);
	  Ctxt *ct = new Ctxt(*ct1);
	  *ct *= (*ct2);
	  *output = (void*) ct;

	  return 0;
}

/**
 * Homomorphic multiplication by a constant of a ciphertext
 */
int mulByConstant(void* pk, void* evk, void** output, void* input,
                     unsigned long constant, unsigned long level) {
	  Ctxt* ct1 = static_cast<Ctxt *>(input);
	  Ctxt *ct = new Ctxt(*ct1);
	  ct->multByConstant(ZZ(constant));
	  *output = (void*) ct;

	  return 0;
}

int mulByConstant(void* pk, void* evk, void** output, void* input,
                     mpz_class const& constant, unsigned long level) {
  assert(false);
  return 1;
}

#else

/**
 * Encrypt a fixed-point number (double) using the public key
*/
int encryptFixedpt(void* pk, void** ciphertext, double message,
                      unsigned long level) {
	  ZZX ply(0);
	  long fplvl=0;
	  double in_fin;
	  long preci = MAX_PRECISION;
	  long PreciInt = -1;
	  long dg = (long) number_of_polynomial_slots(pk);

      flt2plyEncode(ply, fplvl, in_fin, message, preci, PreciInt, base, dg, fracrep);
      CtxtExt *ct = new CtxtExt(*(static_cast<const FHEPubKey *>(pk)));
      ct->lvl = fplvl;
  	  (static_cast<const FHEPubKey *>(pk))->Encrypt(*ct,ply);
  	  *ciphertext = (void *) ct;

  	  return 1;
}

/**
 * Encrypt a fixed-point number (double) using the secret key
*/
int encryptFixedptWithSK(void* sk, void* pk, void** ciphertext,
                            double message, unsigned long level) {
	  ZZX ply;
	  long fplvl=0;
	  double in_fin;
	  long preci = MAX_PRECISION;
	  long PreciInt = -1;
	  long dg = (long) number_of_polynomial_slots(pk);

      flt2plyEncode(ply, fplvl, in_fin, message, preci, PreciInt, base, dg,fracrep);
      CtxtExt *ct = new CtxtExt(*(static_cast<FHESecKey *>(sk)));
      ct->lvl = fplvl;
  	  (static_cast<FHESecKey *>(sk))->Encrypt(*ct,ply);
  	  *ciphertext = (void *) ct;

      return 1;
}

/**
 * Decrypt a ciphertext to a fixed-point number (double)
*/

int decryptFixedpt(void* sk, void* pk, void* ciphertext,
                      double* message, unsigned long level) {
	  ZZX ply;
	  long fplvl=0;
	  long dg = (long) number_of_polynomial_slots(pk);

      (static_cast<FHESecKey *>(sk))->Decrypt(ply, *(static_cast<Ctxt *>(ciphertext)));
      fplvl = (static_cast<CtxtExt *>(ciphertext))->lvl;
      PolyRed(ply,  plaintext_modulus, false);
      ply2fltDecode(*message, ply, fplvl, base, dg, dg, fracrep);

      return 0;
}

/**
 * Homomorphic addition of two ciphertexts
 */
int add(void* pk, void* evk, void** output, void* input1, void* input2,
           unsigned long level1, unsigned long level2) {
  CtxtExt* ct1 = static_cast<CtxtExt *>(input1);
  CtxtExt* ct2 = static_cast<CtxtExt *>(input2);

  if(fracrep)
  {
	  CtxtExt *ct = new CtxtExt(*ct1);
	  *(static_cast<Ctxt *>(ct)) += *(static_cast<Ctxt *>(ct2));
	  ct->lvl = (ct1->lvl >= ct2->lvl)? ct1->lvl: ct2->lvl;			// Values are guaranteed to be non-negative
	  *output = (void*) ct;
  }
  else 				// 'lvl' values are guaranteed to be non-positive
  {
	  if(ct1->lvl <= ct2->lvl)
	  {
		  ZZX ply(0L);
		  SetCoeff(ply,ct2->lvl - ct1->lvl,1);
		  CtxtExt *ct = new CtxtExt(*ct2);
		  (static_cast<Ctxt *>(ct))->multByConstant(ply);
		  *(static_cast<Ctxt *>(ct)) += *(static_cast<Ctxt *>(ct1));
		  ct->lvl = ct1->lvl;
		  *output = (void*) ct;
	  }
	  else
	  {
		  ZZX ply(0L);
		  SetCoeff(ply,ct1->lvl - ct2->lvl,1);
		  CtxtExt *ct = new CtxtExt(*ct1);
		  (static_cast<Ctxt *>(ct))->multByConstant(ply);
		  *(static_cast<Ctxt *>(ct)) += *(static_cast<Ctxt *>(ct2));
		  ct->lvl = ct2->lvl;
		  *output = (void*) ct;
	  }
  }

  return 0;
}

/**
 * Homomorphic multiplication of two ciphertexts
 */
int mul(void* pk, void* evk, void** output, void* input1, void* input2,
           unsigned long level1, unsigned long level2) {
	  CtxtExt* ct1 = static_cast<CtxtExt *>(input1);
	  CtxtExt* ct2 = static_cast<CtxtExt *>(input2);
	  CtxtExt *ct = new CtxtExt(*ct1);
	  *(static_cast<Ctxt *>(ct)) *= *(static_cast<Ctxt *>(ct2));
	  ct->lvl = ct1->lvl + ct2->lvl;

	  *output = (void*) ct;

	  return 0;
}

/**
 * Homomorphic multiplication by a constant of a ciphertext
 */
int mulByConstant(void* pk, void* evk, void** output, void* input,
                     unsigned long constant, unsigned long level) {
	  mulByConstant(pk, evk, output, input, (double)constant, level);
	  return 0;
}

int mulByConstant(void* pk, void* evk, void** output, void* input,
                     mpz_class const& constant, unsigned long level) {
  assert(false);
  return 1;
}

/**
 * Homomorphic multiplication of a ciphertext by a double
 */
int mulByConstant(void* pk, void* evk, void** output, void* input,
                     double constant, unsigned long level) {
	  ZZX ply;
	  long fplvl=0;
	  double in_fin;
	  long preci = MAX_PRECISION;
	  long PreciInt = -1;
	  long dg = (long) number_of_polynomial_slots(pk);

      flt2plyEncode(ply, fplvl, in_fin, constant, preci, PreciInt, base, dg, fracrep);

      CtxtExt* ct1 = static_cast<CtxtExt *>(input);
      CtxtExt *ct = new CtxtExt(*ct1);
      (static_cast<Ctxt *>(ct))->multByConstant(ply);
      ct->lvl += fplvl;
	  *output = (void*) ct;

	  return 0;
}

/**
 * Encrypt in integer (usigned or mpz_class) using the public key
 */
int encryptInteger(void* pk, void** ciphertext, unsigned long message,
                      unsigned long level) {
	  encryptFixedpt(pk, ciphertext, (double) message, level);
	  return 1;
}

int encryptInteger(void* pk, void** ciphertext, mpz_class const& message,
                      unsigned long level) {
	assert(false);
	return 1;
}

/**
 * Encrypt in integer (usigned or mpz_class) using the secret key
 */
int encryptIntegerWithSK(void* sk, void* pk, void** ciphertext,
                            unsigned long message, unsigned long level) {
      encryptFixedptWithSK(sk, pk, ciphertext, (double) message, level);
	  return 1;
}

int encryptIntegerWithSK(void* sk, void* pk, void** ciphertext,
                            mpz_class const& message, unsigned long level) {
	  assert(false);
	  return 1;
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
	  assert(false);
	  return 1;
}

int encryptPolynomial(void* pk, void** ciphertext, mpz_class* message_p,
                         size_t size, unsigned long level) {
	  assert(false);
	  return 1;
}

/**
 * Encrypt a polynomial (of usigned or mpz_class) using the secret key
 */
int encryptPolynomialWithSK(void* sk, void* pk, void** ciphertext,
                               unsigned long* message_p, size_t size,
                               unsigned long level) {
	  assert(false);
	  return 1;
}

int encryptPolynomialWithSK(void* sk, void* pk, void** ciphertext,
                               mpz_class* message_p, size_t size,
                               unsigned long level) {
	  assert(false);
	  return 1;
}


/**
 * Decrypt a ciphertext to an integer (unsigned long of mpz_class)
 */
int decryptInteger(void* sk, void* pk, void* ciphertext,
                      unsigned long* message, unsigned long level) {
	  double tmp_message = (double) *message;
	  decryptFixedpt(sk, pk, ciphertext, &tmp_message, level);
	  *message = (unsigned long) tmp_message;
	  return 1;
}

int decryptInteger(void* sk, void* pk, void* ciphertext, mpz_class* message,
                      unsigned long level) {
	  assert(false);
	  return 1;
}


/**
 * Decrypt a ciphertext to a polynomial (of unsigned long of mpz_class)
 */
int decryptPolynomial(void* sk, void* pk, void* ciphertext,
                         unsigned long** message_p, size_t size,
                         unsigned long level) {
	  assert(false);
	  return 1;
}

int decryptPolynomial(void* sk, void* pk, void* ciphertext,
                         mpz_class** message_p, size_t size,
                         unsigned long level) {
	  assert(false);
	  return 1;
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

#endif

/**
 * Serialize functions (back and forth)
 */

int serialize_parameters  (const char* filename, void* parameters) {  assert(false); return 1;}
int serialize_sk          (const char* filename, void* sk) {  assert(false); return 1;}
int serialize_evk         (const char* filename, void* evk) {  assert(false); return 1;}
int serialize_pk          (const char* filename, void* pk) {  assert(false); return 1;}
int serialize_ciphertext  (const char* filename, void* ciphertext) {  assert(false); return 1;}

int deserialize_parameters(const char* filename, void** parameters) {  assert(false); return 1;}
int deserialize_sk        (const char* filename, void** sk) {  assert(false); return 1;}
int deserialize_evk       (const char* filename, void** evk) {  assert(false); return 1;}
int deserialize_pk        (const char* filename, void** pk) {  assert(false); return 1;}
int deserialize_ciphertext(const char* filename, void** ciphertext) {  assert(false); return 1;}

}
