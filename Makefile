TARGET := svtrek
SRCS := $(wildcard *.c)
OBJS := $(SRCS:.c=.o)

# directories
CURRENT_DIR := $(shell pwd)
DEPS_DIR    := $(CURRENT_DIR)/deps
BIN_DIR     := $(CURRENT_DIR)/bin

# compiler
GXX := gcc
CXXFLAGS = -Wall -Wextra -O3 -pthread
TIME := /usr/bin/time -v

ifdef DEBUG
    CXXFLAGS += -DDEBUG
endif

# abPOA
ABPOA_DIR      := $(DEPS_DIR)/abPOA
ABPOA_CXXFLAGS := -I$(ABPOA_DIR)/include
ABPOA_LDFLAGS  := $(ABPOA_DIR)/lib/libabpoa.a -lm -lz

# htslib
#   HTSLIB=auto    use system htslib if pkg-config finds it, else local  (default)
#   HTSLIB=system  force the system htslib reported by pkg-config
#   HTSLIB=local   force the copy cloned + built under deps/htslib
HTSLIB     ?= auto
HTSLIB_DIR := $(DEPS_DIR)/htslib
HTSLIB_URL := https://github.com/samtools/htslib.git

ifeq ($(HTSLIB),system)
    USE_SYS_HTS := yes
else ifeq ($(HTSLIB),local)
    USE_SYS_HTS := no
else
    USE_SYS_HTS := $(shell pkg-config --exists htslib 2>/dev/null && echo yes || echo no)
endif

ifeq ($(USE_SYS_HTS),yes)
    HTSLIB_CXXFLAGS := $(shell pkg-config --cflags htslib)
    HTSLIB_LDFLAGS  := $(shell pkg-config --libs htslib)
else
    HTSLIB_CXXFLAGS := -I$(HTSLIB_DIR)/include
    HTSLIB_LDFLAGS  := -L$(HTSLIB_DIR)/lib -lhts -Wl,-rpath,$(HTSLIB_DIR)/lib -lz
endif

# build
$(TARGET): $(OBJS)
	$(GXX) $(CXXFLAGS) -o $@ $^ $(HTSLIB_LDFLAGS) $(ABPOA_LDFLAGS)
	rm -f *.o

%.o: %.c
	$(GXX) $(CXXFLAGS) $(HTSLIB_CXXFLAGS) $(ABPOA_CXXFLAGS) -c $< -o $@

# dependencies
.PHONY: install deps deps-abpoa deps-htslib clean distclean

install: deps

deps: deps-abpoa deps-htslib

# abPOA
deps-abpoa:
	@if [ ! -f "$(ABPOA_DIR)/Makefile" ]; then \
		echo ">> initialising abPOA submodule ($(ABPOA_DIR))"; \
		git submodule update --init --recursive -- deps/abPOA; \
	fi
	@echo ">> building abPOA"
	$(MAKE) -C $(ABPOA_DIR)

# htslib
deps-htslib:
ifeq ($(USE_SYS_HTS),yes)
	@echo ">> using system htslib ($$(pkg-config --modversion htslib))"
else
	@if [ ! -f "$(HTSLIB_DIR)/lib/libhts.a" ] && [ ! -f "$(HTSLIB_DIR)/lib/libhts.so" ] && [ ! -f "$(HTSLIB_DIR)/lib/libhts.dylib" ]; then \
		if [ ! -f "$(HTSLIB_DIR)/configure.ac" ] && [ ! -f "$(HTSLIB_DIR)/configure" ]; then \
			echo ">> system htslib not found; cloning into $(HTSLIB_DIR)"; \
			git clone --recurse-submodules $(HTSLIB_URL) $(HTSLIB_DIR); \
		fi; \
		echo ">> building local htslib"; \
		cd $(HTSLIB_DIR) && autoreconf -i && ./configure && $(MAKE) && $(MAKE) prefix=$(HTSLIB_DIR) install; \
	else \
		echo ">> local htslib already built ($(HTSLIB_DIR))"; \
	fi
endif

clean:
	rm -rf *.o $(TARGET)

# remove the locally cloned/built htslib
distclean: clean
	rm -rf $(HTSLIB_DIR)