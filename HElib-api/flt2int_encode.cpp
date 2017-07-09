#include "flt2int_encode.h"

void BbaseDigits(vector<long>& out, long in, long base)
{
	if(base<=0 || in<0)
	{
		cerr << "Invalid base or input" << endl;
		exit(0);
	}

	long x = in;

	while(x > 0)
	{
		BINT_REP(x,base);
		out.push_back(x);
		in -= x;
		in /= base;
		x = in;
	}

	return;
}

/* On input a double value 'in', the value will be rounded off
 * (actually truncated) to the desired precision 'preci'.
 * 'preci' is the number of precision bits for the absolute value (this number must be at most 52 bits, which is the
 * number of significant bits of the mantissa in a double, and two+length(base) many bits less than the length of
 * long integers) and does not include the sign bit. 'PreceiInt' is the number of significant bits reserved for the
 * integer part (excluding the sign bit or for the absolute value of the integer part). If PreciInt < 0, then the routine
 * determines an optimal value for this parameter.
 * The rounded off values will be returned in 'in_fin'.
 * If fracrep=1, then the modulus polynomial used for encoding is X^dg+1. The parameter 'base', usually 3,
 * is the base of the balanced-base encoding to be used. The 'outply' contains the "integer encodings"
 * if fracrep=0, else contains the "fractional" encoding polynomials modulo X^dg+1.
 * The parameter 'fplvl' contains the (non-positive) scaling factors that are needed for the "integer" encoding, or contains an upper bound on the degree
 * of the integer part in case of fractional representation.
 * If fracrep=0, then an error is raised when the degree of the encoding exceeds 'dg'-1. */

void flt2plyEncode(ZZX& outply, long& fplvl, double& in_fin, double in, long preci, long PreciInt, long base, long dg, bool fracrep)
{
		double x = fabs(in);
		long y = long(x);			// Assuming that the precision provided by 'long' suffices
		long ibits, fbits;
		int pos=1;
		vector<long> digits;
		long PrecBbase=0;			// Number of significant fractional 'base' digits

		if(PreciInt > preci)
		{
			cerr << endl << "Error: incorrect precision input" << endl;
			return;
		}

		if(in < 0) pos = -1;

		if(PreciInt < 0)					// If the condition is true, then variables ibits and fbits are alloted minimal possible value
		{
			ibits = NumBits(y);
			fbits = 0L;

			if(ibits >= preci)
			{
				y >>= (ibits-preci);		// Consider only the higher order bits. No round off, only truncate.
				ibits = preci;
				in_fin = pos*y;
			}
			else
				if ((double) y == x)
					in_fin = pos*((double) y);
				else
				{
					while(fbits < preci-ibits && (double)y != x)
					{
						fbits++;
						x *= (1L<<1);
						y = (long) x;
					}
					in_fin = pos*((double) y) / ((double)(1L << fbits));
				}
		}
		else
		{
			ibits = NumBits(y);
			fbits = preci-PreciInt;

			if(ibits > PreciInt)
			{
				y >>= (ibits-PreciInt);		// Consider only the higher order bits. No round off, only truncate. The bits of the fractional part will be set to zero.
				ibits = PreciInt;
				in_fin = pos*y;
				y <<= fbits;
			}
			else
				if ((double) y == x)
				{
					in_fin = pos*((double) y);
					y <<= fbits;
				}
				else
				{
					for(long i=0; i < fbits; i++)
						x *= (1L<<1);
					y = (long) x;
					in_fin = pos*((double) y) / ((double)(1L << fbits));
				}
			ibits = PreciInt;
		}

		double z = in_fin;
		while(fbits < preci-ibits && z != floor(z))
		{
			z *= base;
			fbits++;
		}
		PrecBbase = ceil(log(1L << fbits)/log(base));
		y = pos*in_fin * pow(base,PrecBbase);
		BbaseDigits(digits, y, base);

		if(digits.size() > (size_t)dg)
		{
			cerr << "Error: the number of coefficients exceeds the degree bound." << endl;
			return;
		}

		if(fracrep)
			fplvl = ceil(log((1L << (ibits+1)))/log(base))-1;
		else
			fplvl = -1*PrecBbase;		// If PreciInt < 0, then this value depends only 'preci' and 'preciInt'.

		clear(outply);				// Set all coefficients to 0.

		if(fracrep)
		{
			long flen = ((long) digits.size() < PrecBbase)? digits.size() : PrecBbase;
			for(long j=0; j < flen; j++)
				SetCoeff(outply,dg-PrecBbase+j,-1*pos*digits[j]);
			flen = digits.size() - PrecBbase;
			for(long j=0; j < flen; j++)
				SetCoeff(outply,j,pos*digits[j+PrecBbase]);
		}
		else
			for(size_t j=0; j < digits.size(); j++)
				SetCoeff(outply,j,pos*digits[j]);

return;
}

/* This routine on input an encoded polynomial decodes it into a double floating point type
 * returned in 'out'.
 * If fracrep=true, then 'fplvl' contains a degree upper bound (possibly inclusive) for the integer part, and
 * the modulus polynomial used for encoding is X^dg+1. If fracrep=false, then 'fplvl'
 * contains the (non-positive) scaling factor for the integer encoding. 'PrecBbase' denotes the maximum
 * precision w.r.t. the given balanced base that is tolerated (i.e., the total number of significant digits).
 * The maximum value for PrecBbase = dg.
 * Unutilised significant digits of the integer part are used first for the integer part in case PrecBbase < fplvl,
 * else is utilised when decoding the fractional part.
 * An error is raised when the degree of the encoding exceeds 'dg'-1.*/
void ply2fltDecode(double& out, const ZZX& in, long fplvl, long base, long dg, long PrecBbase, bool fracrep)
{
	if(PrecBbase > dg)
		{
		cerr << endl << "Number of desired significant digits exceeds the number of available coefficients" << endl;
		return;
		}

	if(deg(in) > dg-1 || abs(fplvl) > dg)
	{
		cerr << "Error: the number of coefficients exceeds the degree bound." << endl;
		return;
	}

		long finpos=0, iprec=0, fprec=0;			//fprec undefined for "integer encoding."
		ZZ x = ZZ(0);
		long degi = fplvl;					//Contains an upper bound the degree of the integer part when fracrep=true.
		RR y;

		if(fracrep)
		{
			iprec = fplvl+1;
			if(iprec > PrecBbase)
			{
				iprec = PrecBbase;
				fprec = 0;
			}
			fprec = PrecBbase-iprec;
			if(dg-fplvl-1 < fprec)
				fprec = dg-fplvl-1;
		}
		else
			iprec = (deg(in)+1 > PrecBbase)? PrecBbase : deg(in)+1;

		for(long j=1; j<=iprec; j++)			//Scanning the integer part from the MSD
		{
			ZZ c;

			if(fracrep)
				c = coeff(in,fplvl+1-j);
			else
				c = coeff(in,deg(in)+1-j);

			if(c != ZZ(0))
			{
				while(finpos < j)
				{
					finpos++;
					x *= base;
				}
				x += c;
			}

			if(x==0)										// This will happen only when fracrep=true
			{
				degi--;
				if(iprec >= PrecBbase && iprec < fplvl+1)
					iprec++;								//Unutilised significant digits are transferred to integer part
				else
					fprec++;								//Unutilised significant digits are transferred to fractional part
			}

			if(j==iprec)								    // For trailing zeroes in the integer part
				while(finpos < j)
				{
					finpos++;
					x *= base;
				}
		}

		finpos = 0;

		if(fracrep)											//Scanning the fractional part from the digit after the period when fracrep=true
		{
			for(long j=1; j<=fprec; j++)
			{
				ZZ c = -1*coeff(in,dg-j);
				if(c != ZZ(0))
				{
					while(finpos < j)
					{
						finpos++;
						x *= base;
					}
					x += c;
				}
			}

			out = to_double( to_RR(x) / (power(to_RR(base),(finpos))) );
			if (iprec < degi+1)
				out = to_double( to_RR(x) * (power(to_RR(base), (degi+1-iprec))) );
		}
		else
		{
			out = to_double( to_RR(x) * (power(to_RR(base), (fplvl))) );
			if (iprec < deg(in)+1)
				out *= to_double( (power(to_RR(base), (deg(in)+1-iprec))) );
		}
return;
}

long checkPreci(long preci, long PreciInt, double out, double in)
{
	long x;

	if( PreciInt < 0)
	{
		if(preci >= NumBits((long) fabs(in)))
			x =  (long)(fabs(out-in) * (1L << (preci-NumBits((long) fabs(in)))));
		else
			x =  ((long)(fabs(out-in))) >> (-1 * (preci-NumBits((long) fabs(in))));
	}
	else
	{
		if(preci >= PreciInt)
			x =  (long)(fabs(out-in) * (1L << (preci-PreciInt)));
		else
			x =  ((long)(fabs(out-in))) >> (-1 * (preci-PreciInt));
	}

	return x;
}


void Test_flt2plyEncode()
{
	vector<ZZX> plyEnc;
	vector<long> fplvl;
	vector<double> in, in_fin, out;
	long preci, PreciInt, base, dg, count, precBbase;
	bool fracrep;

	cerr << "Input the precision required (in bits): ";
	cin >> preci;

	cerr << "Input the number of significant bits for the absolute value of the integer part: ";
	cin >> PreciInt;

	cerr << "Input the balanced base: ";
	cin >> base;

	cerr << "Input the degree of the modulus polynomial: ";
	cin >> dg;

	cerr << "Should fractional representation be used (1) or not (0)?: ";
	cin >> fracrep;

	cerr << "Input the vector of doubles: ";

	cin >> in;

	count = in.size();

	plyEnc.resize(count);
	fplvl.resize(count);
	in_fin.resize(count);
	out.resize(count);

	for(long i=0; i<count; i++)
		flt2plyEncode(plyEnc[i], fplvl[i], in_fin[i], in[i], preci, PreciInt, base, dg, fracrep);

	cerr << endl << "Inputs after approximation to the given precision: " << in_fin << endl;

	cerr << endl << "The polynomial encodings are: " << plyEnc << endl;
	cerr << endl << "The fractional part levels are: " << fplvl << endl;

	if(fracrep)
	{
		cerr << endl << "Input degree upper bounds for the integer part: ";
		cin >> fplvl;
	}

	precBbase = dg;

	for(long i=0; i<count; i++)
		ply2fltDecode(out[i],plyEnc[i],fplvl[i],base,dg,precBbase,fracrep);

	cerr << endl << "Inputs after decoding: " << out << endl;

	for(long i=0; i<count; i++)
	{
		long diff = checkPreci(preci, PreciInt, out[i], in_fin[i]);
		if( diff != 0)
			cerr << "Decoding precision error at index " << i << ". Diff = " << diff << endl;
	}

	return;
}


/*
 * This (HEURISTIC) routine encodes a double value by finding the integer coefficients corresponding to the 'n'-th roots of unity.
 * The encoding ring assumed is x^dg+1, where dg is a power of two >= n/2. It is observed in practise that the
 * routine below works best when n=16, 32 or 64 and preci <= 32.
 * The root of unity is chosen to be the canonical value. The other parameters needed for the lattice reduction are also
 * hardcoded. The routine is MOST EFFECTIVE when -1 <= inR, inI <= 1.
 */
void LatfltEncode(long n, ZZX& outply, double inR, double inI, long preci, long dg)
{
	ZZ C;
	if(n < 32)
		C = ZZ(1) << (preci+6);		// Lattice parameters - heuristic
	else
		C = ZZ(1) << (preci+3);		// Lattice parameters - heuristic

	ZZ T  = ZZ(100);				// Lattice parameters - heuristic

	vector<ZZ> A(n,ZZ(0)), B(n,ZZ(0));

	for (long i=0; i<n; i++)
	{
		A[i] = RoundToZZ(to_RR(C)* cos((2*M_PI*i)/(double)n));
		B[i] = RoundToZZ(to_RR(C)* sin((2*M_PI*i)/(double)n));
	}

	ZZ a = RoundToZZ(to_RR(C)*inR);
	ZZ b = RoundToZZ(to_RR(C)*inI);
	Mat<ZZ> M, MT;
	M.SetDims(n+3,n+1);
	MT.SetDims(n+1,n+3);

	for (long i=0; i<n; i++)
	{
		M[i][i] = ZZ(1);
		M[n+1][i] = A[i];
		M[n+2][i] = B[i];
	}
	M[n][n] = T;
	M[n+1][n] = 0-a;
	M[n+2][n] = 0-b;

	MT = transpose(M);
	ZZ dt;
	LLL(dt,MT);			// This NTL routine takes basis vectors as row vectors of the input matrix.
	M = transpose(MT);

	long i=0;
	while (i<=n && (M[n][i] != T && M[n][i]!= 0-T))
		i+=1;
  	if (i > n)
		printf("No desired lattice encoding possible!\n");
	else
	{
        ZZ sig = M[n][i]/T;
        outply = ZZX(0L);
    	for (long j=0; j<n/2; j++)
    		SetCoeff(outply, j*(dg<<1)/n, sig*M[j][i]-sig*M[j+n/2][i]);		// Here we need n <= 2*dg.
	}

	return;
}

/*
 * In this routine the encoding polynomial is evaluated at the canonical 2*dg-th rooth of unity.
 */
void LatfltDecode(double& outR, double& outI, const ZZX& in, long dg)
{
	RR xR(0), xI(0);
	for (long i=0; i<dg; i++)
	{
		xR += to_RR(coeff(in,i))*cos((2*M_PI*i)/(double)(dg <<1));
		xI += to_RR(coeff(in,i))*sin((2*M_PI*i)/(double)(dg <<1));
	}
	outR = to_double(xR);
	outI = to_double(xI);

	return;
}


void Test_LatfltEncode()
{
	ZZX plyEnc;
	double inR, inI, outR, outI;
	long nLat, preci, dg, PreciInt=-1;

	cerr << "Input the precision required (in bits): ";
	cin >> preci;

	cerr << "Input the lattice dimension for lattice reduction: ";
	cin >> nLat;

	cerr << "Input the degree of the modulus polynomial: ";
	cin >> dg;

	cerr << "Input the fractional real part: ";
	cin >> inR;

	cerr << "Input the fractional imaginary part: ";
	cin >> inI;

	LatfltEncode(nLat, plyEnc, inR, inI, preci, dg);

	LatfltDecode(outR, outI, plyEnc, dg);

	cerr << endl << "Real part after decoding: " << outR << ",\tImg part after decoding: " << outI << endl;

	long diff = checkPreci(preci, PreciInt, outR, inR);
	if( diff != 0)
			cerr << "Real part: decoding precision error. Scaled Diff = " << diff << endl;
	diff = checkPreci(preci, PreciInt, outI, inI);
	if( diff != 0)
			cerr << "Img part: decoding precision error. Scaled Diff = " << diff << endl;

	return;
}
