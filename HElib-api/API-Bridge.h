#ifndef API_BRIDGE_H_
#define API_BRIDGE_H_

#include "API.h"
#include "flt2int_encode.h"
#include "findRing.h"
#include "params.h"

#ifdef FXPT

	long base = BAL_BASE;
	bool fracrep = FRAC_REP_ON;

#endif	//FXPT


class CtxtExt: public Ctxt
{
  friend class FHEPubKey;
  friend class FHESecKey;

  public:
	long lvl;

  CtxtExt(const FHEPubKey& newPubKey, long newPtxtSpace=0, long lvlInit=0) // constructor
	: Ctxt(newPubKey, newPtxtSpace),
	  lvl(lvlInit)
	{
	}
};


#endif /* API_BRIDGE_H_ */
