#
# Fortran compiler
#
F90 = gfortran

#
# archiver
#
AR = ar

#
# lapack flag, optional:
# if not provided, local copy will be compiled (not optimized).
#
LIBLAPACK?= -llapack

#=======================
#
# compilation flags
#
FFLAGS =  -fcheck=bounds -ffree-line-length-512 -fPIC -cpp -Ofast -march=native -ffast-math -funroll-loops
#FFLAGS = -fcheck=bounds -ffree-line-length-512 -fPIC -cpp -O0
#FFLAGS += -g -DDEBUG -pedantic -std=f2008

# absolute path of this Makefile
SRC := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
OBJ = $(SRC)/Obj
MOD = $(abspath $(SRC)/../include)
LIBOUT = $(abspath $(SRC)/../lib)

# Hide reference BLAS/LAPACK symbols in the shared library (default).
# IRA_HIDE_BUNDLED_BLAS=0 keeps the bundled export. Callers that used
# dgemm_ from that export must link a BLAS themselves when it is hidden.
# GNU ld takes src/libira.map as --version-script. A linker that rejects
# that flag must not leave dgemm_ exported. gfortran -fvisibility=hidden
# still emits those procedures as global (gfortran 13.3), and
# objcopy --localize-symbol rewrites .symtab while .dynsym stays global
# (binutils 2.42). The fallback archives lap.f and passes --exclude-libs,
# which drops the dynamic exports and keeps calls inside the library.
# Rejecting both flags stops the bundled shared link.
IRA_HIDE_BUNDLED_BLAS ?= 1
CC ?= cc
IRA_VERSION_SCRIPT_FLAG :=
IRA_EXCLUDE_LIBS_FLAG :=
ifneq ($(IRA_HIDE_BUNDLED_BLAS),0)
IRA_VERSION_SCRIPT_FLAG := $(shell echo 'int ira_vs(void){return 0;}' | $(CC) -shared -x c - -o /dev/null -Wl,--version-script=$(SRC)/libira.map >/dev/null 2>&1 && echo -Wl,--version-script=$(SRC)/libira.map)
ifeq ($(strip $(IRA_VERSION_SCRIPT_FLAG)),)
IRA_EXCLUDE_LIBS_FLAG := $(shell echo 'int ira_vs(void){return 0;}' | $(CC) -shared -x c - -o /dev/null -Wl,--exclude-libs,liblap_ref.a >/dev/null 2>&1 && echo -Wl,--exclude-libs,liblap_ref.a)
endif
endif



# Produced filenames
#
IRALIB = libira.a
SHLIB  = libira.so

# version strings
sem_version  := $(shell head -n 1 $(SRC)/../version.txt)
date_version := $(shell tail -n 1 $(SRC)/../version.txt)
IRA_VERSION := -DSTRVERSION=\"${sem_version}\" -DDATEVERSION=${date_version}

# object files
#
# utils
futils := \
	$(SRC)/ira_precision.f90 \
	$(SRC)/timer.f90 \
	$(SRC)/m_ira_error.f90 \
	$(SRC)/m_ira_tools.f90 \
	$(SRC)/m_sofi_tools.f90 \
	$(SRC)/m_ira_version.f90 \
	$(SRC)/sorting_module.f90 \
	$(SRC)/m_set_candidate.f90 \
	$(SRC)/m_cshda.f90 \
	$(SRC)/m_ira_routines.f90 \
	$(SRC)/m_sofi_routines.f90 \
	$(SRC)/ira_mod.f90 \
	$(SRC)/ira_external.f90 \
	$(SRC)/sofi_external.f90

outils:= $(patsubst $(SRC)%,$(OBJ)%,$(futils:%.f90=%.o))

ifndef LIBLAPACK
outils:=$(OBJ)/lap.o $(outils)
endif



flib := $(SRC)/library_ira.f90 \
	$(SRC)/library_sofi.f90
olib := $(patsubst $(SRC)%,$(OBJ)%,$(flib:%.f90=%.o))

# Version-script links lap.o into the shared library. The fallback links
# the archive: --exclude-libs only applies to an archive on the link line.
# Recursive so olib above is visible.
SH_OBJECTS = $(outils) $(olib)
SH_EXTRA =
ifeq ($(strip $(IRA_EXCLUDE_LIBS_FLAG)),)
else
ifndef LIBLAPACK
SH_OBJECTS = $(filter-out $(OBJ)/lap.o,$(outils)) $(olib)
SH_EXTRA = $(IRA_EXCLUDE_LIBS_FLAG) $(OBJ)/liblap_ref.a
endif
endif


default: msg
all : lib shlib

# rules
$(OBJ)/%.o: $(SRC)/%.f90
	$(F90) $(FFLAGS) $(IRA_VERSION) -J$(MOD) -I$(MOD) -c $^ -o $@

# local lapack file if needed (keep -O0 to be sure)
$(OBJ)/lap.o: $(SRC)/../lap_local/lap.f
	$(F90) -J$(MOD) -c -fPIC -O0 -funroll-loops -g $^ -o $@

$(OBJ)/liblap_ref.a: $(OBJ)/lap.o
	rm -f $@
	$(AR) rcs $@ $(OBJ)/lap.o



lib: $(OBJ) $(MOD) $(LIBOUT) $(IRALIB)
shlib: $(OBJ) $(MOD) $(LIBOUT) $(SHLIB)

#
# object folder
#
$(OBJ):
	@if [ ! -d $(OBJ) ]; then mkdir $(OBJ) ; fi
$(LIBOUT):
	@if [ ! -d $(LIBOUT) ]; then mkdir $(LIBOUT) ; fi
$(MOD):
	@if [ ! -d $(MOD) ]; then mkdir $(MOD) ; fi


#====
# pack the library
#====
$(IRALIB): $(outils) $(olib)
	$(AR) -rcv $(LIBOUT)/$@ $^
	@ln -sf $(LIBOUT)/$@ .

#====
# pack the shared library
#====
# LIBLAPACK and the hide flags are part of the link. A newer libira.so
# is not left in place when only those values change: the stamp is
# rewritten, and this target is linked again.
.PHONY: force
force:

$(OBJ)/.shlib-flags: force | $(OBJ)
	@mkdir -p $(OBJ)
	@printf '%s\n' 'liblapack=$(LIBLAPACK) hide=$(IRA_HIDE_BUNDLED_BLAS) vs=$(IRA_VERSION_SCRIPT_FLAG) ex=$(IRA_EXCLUDE_LIBS_FLAG)' > $@.tmp
	@if cmp -s $@.tmp $@; then rm -f $@.tmp; else mv $@.tmp $@; fi

$(SHLIB): $(SH_OBJECTS) $(SRC)/libira.map $(OBJ)/.shlib-flags $(if $(SH_EXTRA),$(OBJ)/liblap_ref.a,)
	@if [ "$(IRA_HIDE_BUNDLED_BLAS)" != "0" ] && [ -z "$(IRA_VERSION_SCRIPT_FLAG)" ]; then \
		if [ -z "$(LIBLAPACK)" ] && [ -z "$(IRA_EXCLUDE_LIBS_FLAG)" ]; then \
			echo "IRA error: linker rejected --version-script and --exclude-libs. Bundled dgemm_ would stay exported. Set IRA_HIDE_BUNDLED_BLAS=0 to keep that export." >&2; \
			exit 1; \
		fi; \
		if [ -z "$(LIBLAPACK)" ]; then \
			echo "IRA warning: linker rejected --version-script. Hiding bundled reference BLAS with --exclude-libs instead of exporting dgemm_." >&2; \
		else \
			echo "IRA warning: linker rejected --version-script. This link names an external BLAS, so bundled lap.f is not in $(SHLIB)." >&2; \
		fi; \
	fi
	$(F90) -o $(LIBOUT)/$@ -shared -J$(MOD) -I$(MOD) $(SH_OBJECTS) $(FFLAGS) $(LIBLAPACK) $(SH_EXTRA) \
		-Wl,-soname,libira.so $(IRA_VERSION_SCRIPT_FLAG)
	@ln -sf $(LIBOUT)/$@ .
	@if [ "$(IRA_HIDE_BUNDLED_BLAS)" != "0" ]; then \
		if nm -D $(LIBOUT)/$@ | awk '$$2 ~ /^[TWD]$$/ && $$3 ~ /^(dgemm_|dgemv_|dgesvd_|xerbla_)$$/ { found = 1 } END { exit !found }'; then \
			echo "error: $(SHLIB) exports reference BLAS/LAPACK symbols"; exit 1; \
		fi; \
	fi

# Shared library with the bundled reference BLAS, then a link of a
# caller that must name its own BLAS. `shlib` is a prerequisite, and
# the recipe links that library again with an empty LIBLAPACK. The
# link-flag stamp is what makes the second link happen when ../lib
# already holds a library built against another BLAS. The result is
# left in ../lib.
.PHONY: test-link-blas
test-link-blas: shlib $(SRC)/libira.map $(SRC)/../tests/link_external_blas.sh
	$(MAKE) shlib LIBLAPACK=
	$(SRC)/../tests/link_external_blas.sh $(LIBOUT)/$(SHLIB)

msg:
	@echo 'Use one of the following:'
	@echo ''
	@echo ' make lib     -> compiles the $(IRALIB) library'
	@echo ' make shlib   -> compiles the $(SHLIB) shared library'
	@echo ' make all     -> all of the above'
	@echo ''
	@echo ' make clean   -> remove all files generated by make'
	@echo ''
	@echo 'Quitting.'




clean:
	rm -rf $(IRALIB) $(SHLIB)
	rm -rf $(OBJ) $(MOD) $(LIBOUT)


