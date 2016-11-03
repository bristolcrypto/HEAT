This Homomorphic Encryption API library integrates the FV-NFLlib and the HElib homomorphic encryption libraries.

To compile this HE-API library, execute the following commands in the current directory:
	./init.sh
	make -f Makefile-HElib-only (currently, 'make' runs into error while compiling with FV-NFLlib)

See the "API.h" file for a brief description of the common API. See the "FV-NFLlib-api/params.h" and the "HElib-api/params.h" files to set the parameters to be used by the FV-NFLlib and the HElib libraries, respectively. 

Contributors:	Tancrède Lepoint (formerly at CryptoExperts) and Srinivas Vivek (University of Bristol). This HE-API library is a contribution of the HEAT project (http://heat-project.eu/).
