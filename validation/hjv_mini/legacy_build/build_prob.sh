#!/bin/bash
# Compile and link one legacy problem file against the patched library.
# Usage: build_prob.sh prob_x.cc [mpi]   -> writes ./prob_x (or ./prob_x_mpi)
set -u
B="$(cd "$(dirname "$0")" && pwd)"
src="$1"; variant="${2:-serial}"; shift 2 2>/dev/null; EXTRA=("$@")
base="$(basename "${src%.*}")"; [ -n "${OUT:-}" ] && base="$OUT"
main="$B/obj/probmain.o"; out="$base"
if [ "$variant" = "mpi" ]; then main="$B/obj/probmain_mpi.o"; out="${base}_mpi"; fi
read -r -a ROOTCFLAGS <<< "$(root-config --cflags)"
read -r -a ROOTLIBS   <<< "$(root-config --libs)"
CXXFLAGS=(-O2 -fPIC -pthread "${ROOTCFLAGS[@]}" -I"$B/include" -I/opt/homebrew/opt/ipopt/include/coin-or -DIPOPT_BUILD)
mpicxx "${CXXFLAGS[@]}" "${EXTRA[@]}" -c "$src" -o "$base.o" || exit 1
mpicxx -O2 "$base.o" "$main" "$B/libfactorana_legacy.a" "${ROOTLIBS[@]}" -lMinuit -lMathCore \
  -L/opt/homebrew/opt/ipopt/lib -lipopt -Wl,-rpath,/opt/homebrew/opt/ipopt/lib -o "$out" 2>&1 | grep -v "^ld: warning"
[ -x "$out" ] && echo "built $out"
