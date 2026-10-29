#!/bin/sh
# libira includes lap.f with every routine prefixed ira_.
# `nm -D --defined-only` must not list the unprefixed BLAS/LAPACK names.
# With no extra argument the script also links a dgemm_ caller against
# -lira -lopenblas and runs it.
set -eu
symbols_only=0
if [ "${1:-}" = "--symbols-only" ]; then
  symbols_only=1
  shift
fi
if [ $# -ne 1 ]; then
  echo "usage: $0 [--symbols-only] /path/to/libira.so" >&2
  exit 2
fi
lib=$1
here=$(CDPATH= cd -- "$(dirname "$0")" && pwd)
src=$here/link_external_blas.c
libdir=$(CDPATH= cd -- "$(dirname "$lib")" && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

nm -D --defined-only "$lib" | awk '{print $NF}' > "$work/defs"
bad=0
while read -r name; do
  [ -n "$name" ] || continue
  if grep -qx "${name}_" "$work/defs" || grep -qx "${name}__" "$work/defs"; then
    echo "error: $lib defines unprefixed ${name}_" >&2
    bad=1
  fi
done << 'EOF'
daxpy
dbdsqr
dcombssq
dcopy
ddot
dgebak
dgebal
dgebd2
dgebrd
dgeev
dgehd2
dgehrd
dgelq2
dgelqf
dgemm
dgemv
dgeqr2
dgeqrf
dger
dgesvd
dhseqr
disnan
dlabad
dlabrd
dlacpy
dladiv
dladiv1
dladiv2
dlaexc
dlahqr
dlahr2
dlaisnan
dlaln2
dlange
dlanv2
dlapy2
dlaqr0
dlaqr1
dlaqr2
dlaqr3
dlaqr4
dlaqr5
dlarf
dlarfb
dlarfg
dlarft
dlarfx
dlartg
dlas2
dlascl
dlaset
dlasq1
dlasq2
dlasq3
dlasq4
dlasq5
dlasq6
dlasr
dlasrt
dlassq
dlasv2
dlasy2
dnrm2
dorg2r
dorgbr
dorghr
dorgl2
dorglq
dorgqr
dorm2r
dormbr
dormhr
dorml2
dormlq
dormqr
drot
dscal
dswap
dtrevc3
dtrexc
dtrmm
dtrmv
idamax
ieeeck
iladlc
iladlr
ilaenv
iparmq
lsame
xerbla
dlamch
dlamc3
EOF
if [ "$bad" -ne 0 ]; then
  exit 1
fi
for sym in ira_dgemm_ ira_dgesvd_ ira_dgeev_ ira_xerbla_; do
  if ! grep -qx "$sym" "$work/defs"; then
    echo "error: $lib does not define $sym (lap_local was not linked)" >&2
    exit 1
  fi
done
echo "symbols: $lib defines no unprefixed BLAS/LAPACK names"
if [ "$symbols_only" -eq 1 ]; then
  exit 0
fi

cc -c "$src" -o "$work/link_external_blas.o"
if cc -o "$work/without_blas" "$work/link_external_blas.o" -L"$libdir" -lira >"$work/without.log" 2>&1; then
  echo "error: dgemm_ resolved without a BLAS library" >&2
  exit 1
fi
if ! grep -q 'dgemm_' "$work/without.log"; then
  echo "error: link without a BLAS failed, but not on dgemm_" >&2
  cat "$work/without.log" >&2
  exit 1
fi
echo "link_external_blas: -lira alone leaves dgemm_ undefined"

cc -o "$work/with_blas" "$work/link_external_blas.o" -L"$libdir" -lira -lopenblas
LD_LIBRARY_PATH="$libdir${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" "$work/with_blas"
echo "link_external_blas: -lira -lopenblas runs"
