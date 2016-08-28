#include "findRing.h"

void FindPM(long& p, long& m, long& d, long k, long L, long c, long p_min, long d_min, long s, long m_min, bool fracrep, bool chosen_p_only, bool verbose)
{
  // In addition to the bounds satisfied in FindM, this function chooses suitable values
  // of p, d, m satisfying given lower bounds on p, d. 's'=1.

  // get a lower-bound on the parameter N=phi(m):
  // 1. Each level in the modulus chain corresponds to pSize=p2Size/2
  //    bits (where we have one prime of this size, and all the others are of
  //    size p2Size).
  //    When using DoubleCRT, we need 2m to divide q-1 for every prime q.
  // 2. With L levels, the largest modulus for "fresh ciphertexts" has size
  //          Q0 ~ p^{L+1} ~ 2^{(L+1)*pSize}
  // 3. We break each ciphertext into upto c digits, do each digit is as large
  //    as    D=2^{(L+1)*pSize/c}
  // 4. The added noise variance term from the key-switching operation is
  //    c*N*sigma^2*D^2, and this must be mod-switched down to w*N (so it is
  //    on par with the added noise from modulus-switching). Hence the ratio
  //    P that we use for mod-switching must satisfy c*N*sigma^2*D^2/P^2<w*N,
  //    or    P > sqrt(c/w) * sigma * 2^{(L+1)*pSize/c}
  // 5. With this extra P factor, the key-switching matrices are defined
  //    relative to a modulus of size
  //          Q0 = q0*P ~ sqrt{c/w} sigma 2^{(L+1)*pSize*(1+1/c)}
  // 6. To get k-bit security we need N>log(Q0/sigma)(k+110)/7.2, i.e. roughly
  //          N > (L+1)*pSize*(1+1/c)(k+110) / 7.2

  // Compute a bound on m, and make sure that it is not too large.

  s=1;
  double cc = 1.0+(1.0/(double)c);
  double dN = ceil((L+1)*FHE_pSize*cc*(k+110)/7.2);
  long N = NTL_SP_BOUND;
  if (N > dN && N > d_min) N = (dN > d_min)? dN : d_min;	// Check on d_min
  else {
    cerr << "Cannot support a bound of " << (dN > d_min)? dN : d_min;
    Error(", aborting.\n");
  }

  m=0; p=0; d=0;

  // find an m and the first p satisfying phi(m)>=N >= d_min, p >= p_min.

  if (m_min)
  {
	  p = p_min;	// It is assumed that in this case a suitable values for p is provided as p_min. The check phi(m)>=N is not done either nor other tests are performed.
	  m = m_min;
	  d = phi_N(m);

	  return;
  }

  if (m==0 && fracrep == false && chosen_p_only == true) {
    // search only for odd values of m, to make phi(m) a little closer to m
	// prime must be used as specified by p_min or else do nothing.
	long candidate_p = (p_min);
//      long lb_cand_m = max(N,candidate_p + 1)|1;
	      long lb_cand_m = max(N,N)|1;		// Ignoring the condition p < m

		for (long candidate_m=lb_cand_m; candidate_m<10*lb_cand_m; candidate_m+=2) {

       if (GCD(candidate_p,candidate_m)!=1) continue;

       d = phi_N(candidate_m); // compute phi(m)
       if (d < N) continue;       // phi(m) too small

       p = candidate_p;  // all tests passed, return this value of p
       m = candidate_m;  // all tests passed, return this value of m
       break;
    }
  }

  if (m==0 && fracrep == false && chosen_p_only == false) {
    // search only for odd values of m, to make phi(m) a little closer to m
	long candidate_p = NextPrime(p_min);
	while (candidate_p < NTL_SP_BOUND)
	{
//      long lb_cand_m = max(N,candidate_p + 1)|1;
	      long lb_cand_m = max(N,N)|1;		// Ignoring the condition p < m

		for (long candidate_m=lb_cand_m; candidate_m<10*lb_cand_m; candidate_m+=2) {

       if (GCD(candidate_p,candidate_m)!=1) continue;

       d = phi_N(candidate_m); // compute phi(m)
       if (d < N) continue;       // phi(m) too small

       p = candidate_p;  // all tests passed, return this value of p
       m = candidate_m;  // all tests passed, return this value of m
       break;
    }
	if (p>0 && m>0) break;
    candidate_p = NextPrime(candidate_p+1);
   }
  }

  if (m==0 && fracrep==true  && chosen_p_only == false) {
    // search only for values of m that are powers of 2 such that ordP=phi(m)
	long candidate_p = NextPrime(p_min);
	while (candidate_p < NTL_SP_BOUND)
	{
//      long candidate_m = NextPowerOfTwo(max(N,candidate_p + 1));
      long candidate_m = (1L << NextPowerOfTwo(max(N,N)));		// Ignoring the condition p < m

      while (candidate_m <= (NTL_SP_BOUND)) {

       if (GCD(candidate_p,candidate_m)!=1) { candidate_m <<= 1; continue;}

       d = phi_N(candidate_m); // compute phi(m)
       if (d < N)
       	   { candidate_m <<= 1; continue;}       // phi(m) too small

       p = candidate_p;  // all tests passed, return this value of p
       m = candidate_m;  // all tests passed, return this value of m
       break;
    }
	if (p>0 && m>0) break;
    candidate_p = NextPrime(candidate_p+1);
   }
  }

  if (m==0 && fracrep==true  && chosen_p_only == true) {
    // search only for values of m that are powers of 2 such that ordP=phi(m)
      long candidate_p = (p_min);
//      long candidate_m = NextPowerOfTwo(max(N,candidate_p + 1));
      long candidate_m = (1L << NextPowerOfTwo(max(N,N)));		// Ignoring the condition p < m

      while (candidate_m <= (NTL_SP_BOUND)) {

       if (GCD(candidate_p,candidate_m)!=1) { candidate_m <<= 1; continue;}

       d = phi_N(candidate_m); // compute phi(m)
       if (d < N)
       	   { candidate_m <<= 1; continue;}       // phi(m) too small

       p = candidate_p;  // all tests passed, return this value of p
       m = candidate_m;  // all tests passed, return this value of m
       break;
    }
  }

  if (verbose) {
    cerr << "*** Bound N="<<N<<", choosing m="<<m <<", p="<< p
         << ", d="<< d << endl;
  }

  return;
}
