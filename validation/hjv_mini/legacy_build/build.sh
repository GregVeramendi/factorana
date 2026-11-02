#!/bin/bash
# Build the legacy factorana C++ code on macOS with Homebrew ROOT 6, Ipopt (MUMPS) and Open MPI.
set -u
cd "$(dirname "$0")"
CXX=mpicxx
ROOTCFLAGS=$(root-config --cflags)
IPOPT_INC=/opt/homebrew/opt/ipopt/include/coin-or
IPOPT_LIB=/opt/homebrew/opt/ipopt/lib
CXXFLAGS="-O2 -fPIC -pthread -Wall -Wno-unused -Wno-sign-compare $ROOTCFLAGS -Iinclude -I$IPOPT_INC -DIPOPT_BUILD"
build_status=0
# (the archive is rebuilt below from the five library objects)
for f in BitFunc.cc gaushermite.C TModel.cc minlkhd_nlp.cc TMinLkhd.cc probmain.cc probmain_mpi.cc; do
  o=obj/${f%.*}.o
  $CXX $CXXFLAGS -c src/$f -o $o > log/${f%.*}.log 2>&1
  rc=$?
  printf "%-16s rc=%d errors=%s\n" "$f" "$rc" "$(grep -c 'error:' log/${f%.*}.log)"
  [ $rc -ne 0 ] && build_status=1
done
exit $build_status
# Archive the library objects (not the probmain* entry points, which build_prob.sh links per problem)
ar rcs libfactorana_legacy.a obj/BitFunc.o obj/gaushermite.o obj/TModel.o obj/minlkhd_nlp.o obj/TMinLkhd.o && echo "archived libfactorana_legacy.a"
