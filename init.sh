#!/bin/sh

if [[ -d "FV-NFLlib" && ! -L "FV-NFLlib" ]] ; then
	echo -e "\n*****Compiling the NFLlib library and the bridge API files.*****\n" 
	cd FV-NFLlib-api
	make
	cd ..
else
	echo -e "\nFV-NFLlib library does not exist.\n If you gave used Git to download this API library, then execute "git submodule init" first and then "git submodule update".\n"
fi

if [[ -d "HElib" && ! -L "HElib" ]] ; then
	echo -e "\n*****Compiling the HElib library and the bridge API files.*****\n"
	cd HElib/src
	make
	cd ../../HElib-api
	make
	cd ../
else
	echo -e "\nHElib library does not exist.\n If you gave used Git to download this API library, then execute "git submodule init" first and then "git submodule update".\n"
fi

