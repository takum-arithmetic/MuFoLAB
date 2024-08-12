# See LICENSE file for copyright and license details
# Takum Linear Algebra Benchmarks
.POSIX:
.SUFFIXES:
.SUFFIXES: .format .jl .output

include config.mk

COMMON =\
	src/format\
	src/QR\
	src/TestMatrices\
	src/TestMatricesGenerator\
	src/Utilities\

EXPERIMENT =\
	src/solve_direct\
	src/solve_qr\

GENERATOR =\
	src/generate_sparse_test_matrices\

all: $(EXPERIMENT:=.output)

src/generate_sparse_test_matrices.output: src/generate_sparse_test_matrices.jl src/TestMatrices.jl
src/solve_direct.output: src/solve_direct.jl src/Utilities.jl
src/solve_qr.output: src/solve_qr.jl src/Utilities.jl src/QR.jl

.jl.format:
	@# work around JuliaFormatter not supporting tabs for indentation
	@# by unexpanding a very wide 16-blank-indent
	$(JULIA) $(JULIAFLAGS) "src/format.jl" "$<" && unexpand -t 16 "$<" > "$<.temp" && mv -f "$<.temp" "$<" && touch "$@"

.jl.output:
	@# experiments print a list of output files, store it an output witness
	$(JULIA) $(JULIAFLAGS) "$<" > "$@.temp" && mv -f "$@.temp" "$@"

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
