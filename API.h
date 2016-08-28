#ifndef _HE_API
#define _HE_API

#include <gmpxx.h>

/**
 * /!\ WARNING /!\
 * 
 * Constant parameters need to be defined in a params.h file
 * 
 * PLAINTEXT_MODULUS (unsigned long)
 * PLAINTEXT_MODULUS_MPZ (string)
 * SECURITY_LEVEL (unsigned long)
 * MULTIPLICATIVE_DEPTH (unsigned long)
 *
 * API for fixed-point arithmetic is conditionally included (see below).
 *
 */

namespace HE {

extern mpz_class plaintext_modulus_mpz;
extern unsigned long plaintext_modulus;
extern unsigned long security_level;
extern unsigned long multiplicative_depth;

extern bool supports_bit_encryption;
extern bool supports_unsigned_encryption;
extern bool supports_mpz_encryption;

extern bool supports_polynomial_encryption;
extern bool supports_vector_encryption;

size_t number_of_polynomial_slots(void* pk);
size_t number_of_vector_slots(void* pk);

/**
 * Init function (uses values in params.h)
 *
 * @param parameters        reference to the pointer to the parameters
 *
 * @return error code
 */
int init(void** parameters);

/**
 * Key Generation
 *
 * @param parameters        pointer to the parameters
 * @param sk                reference to the pointer to the secret key
 * @param pk                reference to the pointer to the public key
 * @param evk               reference to the pointer to the evaluation key
 *
 * @return error code
 */
int keygen(void* parameters, void** sk, void** pk, void** evk);

/**
 * Encrypt in integer (usigned or mpz_class) using the public key
 *
 * @param pk         pointer to the public key
 * @param ciphertext reference to the pointer to the ciphertext
 * @param message    integer message (unsigned long or mpz_class)
 * @param level      (optional) level to encrypt
 *
 * @return error code
 */
int encryptInteger(void* pk, void** ciphertext, unsigned long message,
                      unsigned long level = 1);
int encryptInteger(void* pk, void** ciphertext, mpz_class const& message,
                      unsigned long level = 1);

/**
 * Encrypt in integer (usigned or mpz_class) using the secret key
 *
 * @param sk         pointer to the secret key
 * @param pk         pointer to the public key
 * @param ciphertext reference to the pointer to the ciphertext
 * @param message    integer message (unsigned long or mpz_class)
 * @param level      (optional) level to encrypt
 *
 * @return error code
 */
int encryptIntegerWithSK(void* sk, void* pk, void** ciphertext,
                            unsigned long message, unsigned long level = 1);
int encryptIntegerWithSK(void* sk, void* pk, void** ciphertext,
                            mpz_class const& message, unsigned long level = 1);

/**
 * Encrypt a polynomial (of usigned or mpz_class) using the public key
 *
 * @param pk         pointer to the public key
 * @param ciphertext reference to the pointer to the ciphertext
 * @param message_p  pointer to the polynomial coefficients (of type unsigned
 *long or mpz_class)
 * @param size     number of coefficients in the polynomial
 * @param level      (optional) level to encrypt
 *
 * @return error code
 */
int encryptPolynomial(void* pk, void** ciphertext, unsigned long* message_p,
                         size_t size, unsigned long level = 1);
int encryptPolynomial(void* pk, void** ciphertext, mpz_class* message_p,
                         size_t size, unsigned long level = 1);

/**
 * Encrypt a polynomial (of usigned or mpz_class) using the secret key
 *
 * @param sk         pointer to the secret key
 * @param pk         pointer to the public key
 * @param ciphertext reference to the pointer to the ciphertext
 * @param message_p  pointer to the polynomial coefficients (of type unsigned
 *long or mpz_class)
 * @param size     number of coefficients in the polynomial
 * @param level      (optional) level to encrypt
 *
 * @return error code
 */
int encryptPolynomialWithSK(void* sk, void* pk, void** ciphertext,
                               unsigned long* message_p, size_t size,
                               unsigned long level = 1);
int encryptPolynomialWithSK(void* sk, void* pk, void** ciphertext,
                               mpz_class* message_p, size_t size,
                               unsigned long level = 1);

/**
 * Encrypt a vector (of usigned or mpz_class) using the public key [BATCHING]
 *
 * @param pk         pointer to the public key
 * @param ciphertext reference to the pointer to the ciphertext
 * @param message_v  pointer to the vector coefficients (of type unsigned long
 *or mpz_class)
 * @param size     number of coefficients
 * @param level      (optional) level to encrypt
 *
 * @return error code
 */
int encryptVector(void* pk, void** ciphertext, unsigned long* message_v,
                     size_t size, unsigned long level = 1);
int encryptVector(void* pk, void** ciphertext, mpz_class* message_v,
                     size_t size, unsigned long level = 1);

/**
 * Encrypt a vector (of usigned or mpz_class) using the secret key [BATCHING]
 *
 * @param sk         pointer to the secret key
 * @param pk         pointer to the public key
 * @param ciphertext reference to the pointer to the ciphertext
 * @param message_v  pointer to the vector coefficients (of type unsigned long
 *or mpz_class)
 * @param size     number of coefficients
 * @param level      (optional) level to encrypt
 *
 * @return error code
 */
int encryptVectorWithSK(void* sk, void* pk, void** ciphertext,
                           unsigned long* message_v, size_t size,
                           unsigned long level = 1);
int encryptVectorWithSK(void* sk, void* pk, void** ciphertext,
                           mpz_class* message_v, size_t size,
                           unsigned long level = 1);

/**
 * Decrypt a ciphertext to an integer (unsigned long of mpz_class)
 *
 * @param sk         pointer to the secret key
 * @param pk         pointer to the public key
 * @param ciphertext pointer to the ciphertext
 * @param message    reference to the message
 * @param level      (optional) level to decrypt
 *
 * @return error code
 */
int decryptInteger(void* sk, void* pk, void* ciphertext,
                      unsigned long* message, unsigned long level = 1);
int decryptInteger(void* sk, void* pk, void* ciphertext, mpz_class* message,
                      unsigned long level = 1);

/**
 * Decrypt a ciphertext to a polynomial (of unsigned long of mpz_class)
 *
 * @param sk         pointer to the secret key
 * @param pk         pointer to the public key
 * @param ciphertext pointer to the ciphertext
 * @param message_p  reference to the pointer to the polynomial coefficients
 * @param size     number of coefficients in the polynomial
 * @param level      (optional) level to decrypt
 *
 * @return error code
 */
int decryptPolynomial(void* sk, void* pk, void* ciphertext,
                         unsigned long** message_p, size_t size,
                         unsigned long level = 1);
int decryptPolynomial(void* sk, void* pk, void* ciphertext,
                         mpz_class** message_p, size_t size,
                         unsigned long level = 1);

/**
 * Decrypt a ciphertext to a vector (of unsigned long of mpz_class) [BATCHING]
 *
 * @param sk         pointer to the secret key
 * @param pk         pointer to the public key
 * @param ciphertext pointer to the ciphertext
 * @param message_p  reference to the pointer to the vector
 * @param size     number of coefficients
 * @param level      (optional) level to decrypt
 *
 * @return error code
 */
int decryptVector(void* sk, void* pk, void* ciphertext,
                     unsigned long** message_v, size_t size,
                     unsigned long level = 1);
int decryptVector(void* sk, void* pk, void* ciphertext,
                     mpz_class** message_v, size_t size,
                     unsigned long level = 1);

/**
 * Homomorphic addition of two ciphertexts
 *
 * @param pk     pointer to the public key
 * @param evk    pointer to the evaluation key
 * @param output reference to the pointer to the output ciphertext
 * @param input1 pointer to the input1 ciphertext
 * @param input2 pointer to the input2 ciphertext
 * @param level1 (optional) level of the first input
 * @param level2 (optional) level of the second input
 *
 * @return error code
 */
int add(void* pk, void* evk, void** output, void* input1, void* input2,
           unsigned long level1 = 1, unsigned long level2 = 1);

/**
 * Homomorphic multiplication of two ciphertexts
 *
 * @param pk     pointer to the public key
 * @param evk    pointer to the evaluation key
 * @param output reference to the pointer to the output ciphertext
 * @param input1 pointer to the input1 ciphertext
 * @param input2 pointer to the input2 ciphertext
 * @param level1 (optional) level of the first input
 * @param level2 (optional) level of the second input
 *
 * @return error code
 */
int mul(void* pk, void* evk, void** output, void* input1, void* input2,
           unsigned long level1 = 1, unsigned long level2 = 1);

/**
 * Homomorphic multiplication by a constant of a ciphertext
 *
 * @param pk       pointer to the public key
 * @param evk      pointer to the evaluation key
 * @param output   reference to the pointer to the output ciphertext
 * @param input    pointer to the input ciphertext
 * @param constant constant (unsigned long or mpz_class)
 * @param level    (optional) level of the input
 *
 * @return error code
 */
int mulByConstant(void* pk, void* evk, void** output, void* input,
                     unsigned long constant, unsigned long level = 1);
int mulByConstant(void* pk, void* evk, void** output, void* input,
                     mpz_class const& constant, unsigned long level = 1);

/**
 * Serialize functions (back and forth)
 *
 * @param file name
 * @param element to serialize or reference to element to deserialize
 * 
 * @return error code
 */
int serialize_parameters  (const char* filename, void* parameters);
int serialize_sk          (const char* filename, void* sk);
int serialize_evk         (const char* filename, void* evk);
int serialize_pk          (const char* filename, void* pk);
int serialize_ciphertext  (const char* filename, void* ciphertext);

int deserialize_parameters(const char* filename, void** parameters);
int deserialize_sk        (const char* filename, void** sk);
int deserialize_evk       (const char* filename, void** evk);
int deserialize_pk        (const char* filename, void** pk);
int deserialize_ciphertext(const char* filename, void** ciphertext);

#ifdef FXPT

extern bool supports_fixedpt_encryption;

/**
 * Encrypt a fixed-point number (double) using the public key
 *
 * @param pk         pointer to the public key
 * @param ciphertext reference to the pointer to the ciphertext
 * @param message    fixed-point number message (double)
 * @param level      (optional) level to encrypt
 *
 * @return error code
 */
int encryptFixedpt(void* pk, void** ciphertext, double message,
                      unsigned long level = 1);

/**
 * Encrypt a fixed-point number (double) using the secret key
 *
 * @param sk         pointer to the secret key
 * @param pk         pointer to the public key
 * @param ciphertext reference to the pointer to the ciphertext
 * @param message    fixed-point number message (double)
 * @param level      (optional) level to encrypt
 *
 * @return error code
 */
int encryptFixedptWithSK(void* sk, void* pk, void** ciphertext,
                            double message, unsigned long level = 1);

/**
 * Decrypt a ciphertext to a fixed-point number (double)
 *
 * @param sk         pointer to the secret key
 * @param pk         pointer to the public key
 * @param ciphertext pointer to the ciphertext
 * @param message    reference to the fixed-point number message
 * @param level      (optional) level to decrypt
 *
 * @return error code
 */
int decryptFixedpt(void* sk, void* pk, void* ciphertext,
                      double* message, unsigned long level = 1);


int mulByConstant(void* pk, void* evk, void** output, void* input,
                     double constant, unsigned long level);


#endif //FXPT

}

#endif
