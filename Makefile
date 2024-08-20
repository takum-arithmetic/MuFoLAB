# See LICENSE file for copyright and license details
# Takum Linear Algebra Benchmarks
.POSIX:
.SUFFIXES:
.SUFFIXES: .format .jl .output

include config.mk

COMMON =\
	src/Crutches\
	src/Experiments\
	src/format\
	src/QR\
	src/LU\
	src/TestMatrices\
	src/TestMatricesGenerator\

EXPERIMENT =\
	src/solve_lu\
	src/solve_qr\

GENERATOR =\
	src/generate_sparse_test_matrices\

all: $(EXPERIMENT:=.output)

src/generate_sparse_test_matrices.output: src/generate_sparse_test_matrices.jl src/TestMatrices.jl config.mk Makefile
src/solve_lu.output: src/solve_lu.jl src/Experiments.jl src/LU.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output config.mk Makefile
src/solve_qr.output: src/solve_qr.jl src/Experiments.jl src/QR.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output config.mk Makefile

.jl.format:
	@# work around JuliaFormatter not supporting tabs for indentation
	@# by unexpanding a very wide 16-blank-indent
	$(JULIA) $(JULIA_FLAGS) -- "src/format.jl" "$<" && unexpand -t 16 "$<" > "$<.temp" && mv -f "$<.temp" "$<" && touch "$@"

.jl.output:
	@# experiments print a list of output files, store it an output witness
	$(JULIA) $(JULIA_FLAGS) -- "$<" $(JULIA_SCRIPT_FLAGS) > "$@.temp" && mv -f "$@.temp" "$@"

clean:
	@# use the output witnesses to clean up the output files, except
	@# those from the generators
	for w in $(EXPERIMENT:=.output); do if [ -f "$$w" ]; then xargs rm -f < "$$w"; fi; done
	rm -f $(EXPERIMENT:=.output) $(EXPERIMENT:=.output.temp)
	rm -f $(COMMON:=.format) $(EXPERIMENT:=.format)

clean-generated:
	@# remove the generated files using the output witnesses
	for w in $(GENERATOR:=.output); do if [ -f "$$w" ]; then xargs rm -f < "$$w"; fi; done

format: $(COMMON:=.format) $(EXPERIMENT:=.format) $(GENERATOR:=.format)

.PHONY: all clean format
