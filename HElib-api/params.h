#ifndef API_PARAMS_H_
#define API_PARAMS_H_

#ifdef BINARY
#define PLAINTEXT_MODULUS 2
#define PLAINTEXT_MODULUS_MPZ "2"
#else
#define PLAINTEXT_MODULUS 5003
#define PLAINTEXT_MODULUS_MPZ "5003"
#endif

#define SECURITY_LEVEL 80
#define MULTIPLICATIVE_DEPTH 20

#define	MAX_PRECISION 8		// precision (in bits) for the absolute value (excludes the sign bit)
#define FRAC_REP_ON	1	//currently only the "fractional representation" of fixed-point numbers is supported
#define BAL_BASE 3

#ifdef FXPT
#define MIN_RING_DEG 100
#else
#define MIN_RING_DEG 1
#endif

#endif /* API_PARAMS_H_ */


