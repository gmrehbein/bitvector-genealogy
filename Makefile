# Top-level Makefile for BitVector Genealogy Project
# Provides targets to compile both the C++ and Rust versions.

.PHONY: all build-cpp build-rs run-cpp run-rs

TEST_DATA := test.data
PARENT := parent.txt
ZSCORE ?= 3.0

# 'all' generates data, compiles both C++ (cpp) and Rust (rs) binaries, then runs and compares
all: generate build-cpp build-rs run-cpp run-rs compare

# 1) generate the data
generate:
	@echo "Generating test data and true parent mapping..."
	$(PYTHON) python/populate.py
	@echo "Wrote $(TEST_DATA) & $(PARENT)"

# Build the C++ version using Meson (configurable compiler)
# Allow overriding C and C++ compilers; defaults to system default (on macOS, clang)
# To compile with gcc: make build-cpp CXX=g++
CC ?= cc
CXX ?= clang++

# Build the C++ version using Meson
build-cpp:
	@echo "Building C++ version..."
	cd cpp && \
	CC=$(CC) CXX=$(CXX) meson setup builddir --buildtype=release && \
	meson compile -C builddir

# Build the Rust version using Cargo
build-rs:
	@echo "Building Rust version..."
	cd rs && \
	cargo build --release

# Run the C++ version
run-cpp: build-cpp
	@echo "Running C++ version..."
	./cpp/builddir/bvg-cpp --zscore $(ZSCORE) $(TEST_DATA) > out_cpp.txt

# Run the Rust version
run-rs: build-rs
	@echo "Running Rust version..."
	./rs/target/release/bvg-rs -i $(TEST_DATA) -z $(ZSCORE) -o out_rs.txt

# Compare both outputs to the true parent.txt
compare:
	@echo "Comparing C++ output..."
	$(PYTHON) python/compare_parents.py $(PARENT) out_cpp.txt
	@echo "Comparing Rust output..."
	$(PYTHON) python/compare_parents.py $(PARENT) out_rs.txt

# Clean build artifacts and data
clean:
	rm -rf cpp/builddir rs/target
	rm -f out_cpp.txt out_rs.txt $(TEST_DATA) $(PARENT)

# allow overriding ZS
zs:
	@echo "Current z-score = $(ZS)"

# Benchmark execution time of C++ and Rust versions
benchmark: build-cpp build-rs
	@echo "Benchmarking C++ version..."
	@time -p ./cpp/builddir/bvg-cpp --zscore $(ZSCORE) $(TEST_DATA) >/dev/null
	@echo "Benchmarking Rust version..."
	@time -p ./rs/target/release/bvg-rs -i $(TEST_DATA) -z $(ZSCORE) -o /dev/null