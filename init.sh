#!/bin/sh

cd helib/
make
cd api-bridge/
make
cd ../..
make

