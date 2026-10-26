#!/bin/sh
# Link tests/link_external_blas.c against one libira.
# The hide is in force when dgemm_ is not a dynamic definition: the
# object then fails to link with -lira alone and links with a BLAS.
set -eu
if [ $# -ne 1 ]; then
  echo "usage: $0 /path/to/libira.so" >&2
  exit 2
fi
lib=$1
here=$(CDPATH= cd -- "$(dirname "$0")" && pwd)
src=$here/link_external_blas.c
libdir=$(CDPATH= cd -- "$(dirname "$lib")" && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

if nm -D "$lib" | awk '$2 ~ /^[TWD]$/ && $3 ~ /^(dgemm_|dgemv_|dgesvd_|xerbla_)$/ { found = 1 } END { exit !found }'; then
  echo "error: $lib exports reference BLAS/LAPACK symbols" >&2
  exit 1
fi

cc -c "$src" -o "$work/link_external_blas.o"
if cc -o "$work/without_blas" "$work/link_external_blas.o" -L"$libdir" -lira >"$work/without.log" 2>&1; then
  if ! readelf -d "$lib" | grep -Eq 'NEEDED.*blas|NEEDED.*lapack'; then
    echo "error: dgemm_ resolved from $lib, which does not need a BLAS" >&2
    exit 1
  fi
  echo "link_external_blas: $lib does not export reference BLAS; the executable took dgemm_ from libira's NEEDED BLAS"
else
  if ! grep -q 'dgemm_' "$work/without.log"; then
    echo "error: link without a BLAS failed, but not on dgemm_" >&2
    cat "$work/without.log" >&2
    exit 1
  fi
  echo "link_external_blas: $lib hides reference BLAS; -lira alone leaves dgemm_ undefined"
fi

cc -o "$work/with_blas" "$work/link_external_blas.o" -L"$libdir" -lira -lopenblas
LD_LIBRARY_PATH="$libdir${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" "$work/with_blas"
echo "link_external_blas: -lira -lopenblas runs"
