# gts-cxx: plain make, no configuration step.
#
#   make            the library, lib/libgts-cxx.so, and the test programs
#   make lib        the library only
#   make install    the headers and the library into PREFIX (/usr/local)
#   make check      the unit tests, under AddressSanitizer
#   make ref-check  runs the comparison cases through the C reference only
#   make compare    runs the cases through both and compares the results
#
# The library has no dependencies. Only the reference program needs the C
# GTS library and glib, found through pkg-config. It casts callbacks to GTS's
# generic function types, as GTS itself does, so the warning about such casts
# is turned off for it.

PREFIX   ?= /usr/local

CXX      ?= g++
CXXFLAGS ?= -O2 -g
CXXFLAGS += -std=c++23 -Wall -Wextra
CPPFLAGS += -Iinclude

GTS_CFLAGS := $(shell pkg-config --cflags gts)
GTS_LIBS   := $(shell pkg-config --libs gts)

LIB     := lib/libgts-cxx.so
SRCS    := $(wildcard src/*.cxx)
OBJS    := $(SRCS:src/%.cxx=build/%.o)
HEADERS := $(wildcard include/gts/*.h)

REF  := test/compare/gts-ref
CXXT := test/compare/gts-cxx
UNIT := test/unit/test-core

.PHONY: all lib install check ref ref-check compare clean

all: $(LIB) $(REF) $(CXXT) $(UNIT)

lib: $(LIB)

install: $(LIB)
	install -d $(DESTDIR)$(PREFIX)/include/gts $(DESTDIR)$(PREFIX)/lib
	install -m 644 $(HEADERS) $(DESTDIR)$(PREFIX)/include/gts
	install -m 755 $(LIB) $(DESTDIR)$(PREFIX)/lib

build/%.o: src/%.cxx $(HEADERS)
	@mkdir -p build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -fPIC -c -o $@ $<

$(LIB): $(OBJS)
	@mkdir -p lib
	$(CXX) $(CXXFLAGS) -shared -o $@ $^

ref: $(REF)

$(REF): test/compare/gts-ref.cxx
	$(CXX) $(CXXFLAGS) -Wno-cast-function-type $(GTS_CFLAGS) -o $@ $< $(GTS_LIBS)

$(CXXT): test/compare/gts-cxx.cxx $(LIB) $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< -Llib -lgts-cxx -Wl,-rpath,'$$ORIGIN/../../lib'

# The unit test compiles the library sources itself, under AddressSanitizer.
$(UNIT): test/unit/test-core.cxx $(SRCS) $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -O1 -fsanitize=address -fno-omit-frame-pointer -o $@ $< $(SRCS)

check: $(UNIT)
	$(UNIT)

ref-check: $(REF)
	test/compare/run-ref.sh

compare: $(REF) $(CXXT)
	test/compare/compare.sh

clean:
	rm -rf build lib $(REF) $(CXXT) $(UNIT)
	rm -rf test/compare/out
