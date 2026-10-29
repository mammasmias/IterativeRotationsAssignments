#!/bin/sh
# Build libira with make and with CMake, then check the shared library
# exports no unprefixed BLAS/LAPACK names and links with OpenBLAS.
set -eu
root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
tar -C "$root" \
  --exclude .git \
  --exclude build \
  --exclude lib \
  --exclude include \
  --exclude src/Obj \
  -cf - . | tar -C "$work" -xf -

echo "make: all"
make -C "$work/src" all
"$work/tests/link_external_blas.sh" "$work/lib/libira.so"

echo "cmake: shared library and tests"
cmake -S "$work" -B "$work/build" \
  -DBUILD_SHARED_LIBS=ON \
  -DIRA_BUILD_TESTS=ON \
  -DIRA_ENABLE_NATIVE_OPTIMIZATION=ON
cmake --build "$work/build"
so=$(find "$work/build" -name 'libira.so' -o -name 'libira.so.*' | head -n 1)
if [ -z "$so" ]; then
  echo "cmake produced no libira.so" >&2
  exit 1
fi
"$work/tests/link_external_blas.sh" "$so"
( cd "$work/build" && ctest --output-on-failure )
echo "check_builds: ok"
