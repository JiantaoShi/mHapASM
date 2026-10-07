# Build mhapasm and mhapconvert against htslib.
#   make                          # htslib found automatically (Homebrew, conda, /usr/local)
#   make HTSLIB=/path/to/prefix   # prefix holding include/htslib and lib/libhts.*
#   make static HTSLIB=...        # link libhts.a (needs the static zlib/bz2/lzma/deflate libs)

CC      ?= cc
CFLAGS  ?= -O2 -g -Wall -Wextra -Wno-unused-parameter
HTSLIB  ?= $(firstword $(foreach d,$(CONDA_PREFIX) /opt/homebrew/opt/htslib /usr/local/opt/htslib /usr/local /usr,$(if $(wildcard $(d)/include/htslib/sam.h),$(d))))

ifeq ($(HTSLIB),)
$(error htslib not found; run 'make HTSLIB=<prefix>')
endif

CPPFLAGS += -I$(HTSLIB)/include
LDFLAGS  += -L$(HTSLIB)/lib -Wl,-rpath,$(HTSLIB)/lib
LDLIBS   += -lhts -lz -lm -lpthread
PROGS     = mhapasm mhapconvert

all: $(PROGS)

$(PROGS): %: src/%.c src/bsread.c src/bsread.h
	$(CC) -std=gnu99 $(CFLAGS) $(CPPFLAGS) -o $@ src/$@.c src/bsread.c $(LDFLAGS) $(LDLIBS)

static: src/mhapasm.c src/mhapconvert.c src/bsread.c src/bsread.h
	for p in $(PROGS); do \
	    $(CC) -std=gnu99 $(CFLAGS) $(CPPFLAGS) -o $$p src/$$p.c src/bsread.c $(HTSLIB)/lib/libhts.a \
	        -lz -lbz2 -llzma -ldeflate -lcurl -lm -lpthread || exit 1; \
	done

test: $(PROGS)
	bash test/run_tests.sh

clean:
	rm -rf $(PROGS) $(addsuffix .dSYM,$(PROGS))

.PHONY: all static test clean
