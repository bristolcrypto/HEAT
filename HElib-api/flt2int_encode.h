#ifndef FIXEDPOINT_FLT2INT_ENCODE_H_
#define FIXEDPOINT_FLT2INT_ENCODE_H_

#include <NTL/ZZX.h>
#include <NTL/RR.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>
#include "FHE.h"
#include "timing.h"
#include <vector>
#include <cassert>
#include <cstdio>
#include <cmath>
#include <ctime>

#define _USE_MATH_DEFINES

/*If y<=2, then x = x%y, else x will be in (-y/2,...,y/2] */
#define BINT_REP(x,y){\
		x %=(y);\
		if((y)>2) if ((x)>((y)>>1)) x-=(y);\
		}

/* Polynomial of type ZZX is converted to type ZZ_pX */
#define ZZX2PX(a,b,c) {\
		a = ZZ_pX(0L);\
        for(long def_i=0; def_i <= deg(b); def_i++)\
	       SetCoeff(a,def_i,coeff(b,def_i)%c);\
        }

/* Polynomial of type ZZ_pX is converted to type ZZX */
#define ZZPX2X(a,b) {\
		a = ZZX(0L);\
        for(long def_i=0; def_i <= deg(b); def_i++)\
	       SetCoeff(a,def_i,rep(coeff(b,def_i)));\
        }

void BbaseDigits(vector<long>& out, long in, long base);

void flt2plyEncode(ZZX& outply, long& fplvl, double& in_fin, double in, long preci, long PreciInt, long base, long dg, bool fracrep)
;

void ply2fltDecode(double& out, const ZZX& in, long fplvl, long base, long dg, long PrecBbase, bool fracrep);

long checkPreci(long preci, long PreciInt, double out, double in);

void Test_flt2plyEncode();

void LatfltEncode(long n, ZZX& outply, double inR, double inI, long preci, long dg);

void LatfltDecode(double& outR, double& outI, const ZZX& in, long dg);

void Test_LatfltEncode();

#endif /* FIXEDPOINT_FLT2INT_ENCODE_H_ */
