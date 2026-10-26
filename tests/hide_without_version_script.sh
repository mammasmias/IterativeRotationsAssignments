#!/bin/sh
# A compiler wrapper that rejects --version-script must still hide
# dgemm_, dgemv_, dgesvd_, and xerbla_. The same tree's test-link-blas
# must relink when libira.so was built with a different LIBLAPACK.
set -eu
root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
git -C "$root" archive HEAD | tar -x -C "$work"

cat > "$work/wrap.sh" << 'EOF'
#!/bin/sh
real=${IRA_REAL_COMPILER:?}
for arg in "$@"; do
  case $arg in
    *version-script*)
      echo "rejecting $arg" >&2
      exit 1
      ;;
  esac
done
exec "$real" "$@"
EOF
chmod +x "$work/wrap.sh"

echo "make: version-script path"
if ! make -C "$work/src" test-link-blas >"$work/make-vs.log" 2>&1; then
  cat "$work/make-vs.log" >&2
  exit 1
fi
"$work/tests/link_external_blas.sh" "$work/lib/libira.so" >/dev/null
echo "make: identical bundled link is a no-op"
make -C "$work/src" shlib LIBLAPACK= >"$work/make-noop.log" 2>&1
if grep -q ' -shared ' "$work/make-noop.log"; then
  echo "identical shlib relinked" >&2
  cat "$work/make-noop.log" >&2
  exit 1
fi

make -C "$work/src" clean >/dev/null
echo "make: rejecting wrapper"
log=$work/make-wrapper.log
if ! env IRA_REAL_COMPILER=gcc make -C "$work/src" test-link-blas CC="$work/wrap.sh" >"$log" 2>&1; then
  cat "$log" >&2
  exit 1
fi
if ! grep -q 'IRA warning: linker rejected --version-script' "$log"; then
  echo "make wrapper did not warn" >&2
  cat "$log" >&2
  exit 1
fi
if grep -q 'Nothing to be done' "$log"; then
  echo "fresh test-link-blas skipped the link" >&2
  exit 1
fi
"$work/tests/link_external_blas.sh" "$work/lib/libira.so" >/dev/null

echo "make: relink after an external BLAS"
env IRA_REAL_COMPILER=gcc make -C "$work/src" shlib CC="$work/wrap.sh" LIBLAPACK=-lopenblas >"$work/make-openblas.log" 2>&1
if ! grep -q 'IRA warning: linker rejected --version-script' "$work/make-openblas.log"; then
  echo "external BLAS link did not warn" >&2
  exit 1
fi
env IRA_REAL_COMPILER=gcc make -C "$work/src" test-link-blas CC="$work/wrap.sh" >"$work/make-relink.log" 2>&1
if grep -q 'Nothing to be done' "$work/make-relink.log"; then
  echo "test-link-blas left the OpenBLAS library in place" >&2
  cat "$work/make-relink.log" >&2
  exit 1
fi
if readelf -d "$work/lib/libira.so" | grep -Eq 'NEEDED.*blas|NEEDED.*lapack'; then
  echo "relinked library still needs an external BLAS" >&2
  exit 1
fi
"$work/tests/link_external_blas.sh" "$work/lib/libira.so" >/dev/null

echo "cmake: rejecting wrapper"
cat > "$work/wrap-cc.sh" << 'EOF'
#!/bin/sh
export IRA_REAL_COMPILER=/usr/bin/gcc
exec "$(CDPATH= cd -- "$(dirname "$0")" && pwd)/wrap.sh" "$@"
EOF
cat > "$work/wrap-fc.sh" << 'EOF'
#!/bin/sh
export IRA_REAL_COMPILER=/usr/bin/gfortran
exec "$(CDPATH= cd -- "$(dirname "$0")" && pwd)/wrap.sh" "$@"
EOF
chmod +x "$work/wrap-cc.sh" "$work/wrap-fc.sh"
rm -rf "$work/build"
if ! cmake -S "$work" -B "$work/build" \
  -DBUILD_SHARED_LIBS=ON \
  -DCMAKE_DISABLE_FIND_PACKAGE_LAPACK=ON \
  -DIRA_BUILD_TESTS=ON \
  -DCMAKE_C_COMPILER="$work/wrap-cc.sh" \
  -DCMAKE_Fortran_COMPILER="$work/wrap-fc.sh" \
  >"$work/cmake-cfg.log" 2>&1; then
  cat "$work/cmake-cfg.log" >&2
  exit 1
fi
if ! grep -q 'IRA warning: linker rejected --version-script' "$work/cmake-cfg.log"; then
  echo "cmake did not warn" >&2
  cat "$work/cmake-cfg.log" >&2
  exit 1
fi
if ! cmake --build "$work/build" >"$work/cmake-build.log" 2>&1; then
  cat "$work/cmake-build.log" >&2
  exit 1
fi
so=$(find "$work/build" -name 'libira.so' -o -name 'libira.so.*' | head -n 1)
if [ -z "$so" ]; then
  echo "cmake produced no libira.so" >&2
  exit 1
fi
"$work/tests/link_external_blas.sh" "$so" >/dev/null
if ! (cd "$work/build" && ctest --output-on-failure >"$work/ctest.log" 2>&1); then
  cat "$work/ctest.log" >&2
  exit 1
fi
if ! grep -q 'link_external_blas' "$work/ctest.log"; then
  echo "ctest did not run link_external_blas" >&2
  cat "$work/ctest.log" >&2
  exit 1
fi
echo "hide_without_version_script: ok"
