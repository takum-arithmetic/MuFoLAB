# See LICENSE file for copyright and license details
# Takum Linear Algebra Benchmarks
.POSIX:
.SUFFIXES:
.SUFFIXES: .format .jl .output .output_sorted

include config.mk

COMMON =\
	src/Crutches\
	src/Experiments\
	src/Float128Conversions\
	src/format\
	src/QR\
	src/LU\
	src/sort_csv\
	src/TestMatrices\
	src/TestMatricesGenerator\

EXPERIMENT =\
	src/solve_gmres_ilu\
	src/solve_lu\
	src/solve_mpir_float_08_16_32\
	src/solve_mpir_posit_08_16_32\
	src/solve_mpir_takum_08_16_32\
	src/solve_mpir_float_16_16_32\
	src/solve_mpir_posit_16_16_32\
	src/solve_mpir_takum_16_16_32\
	src/solve_mpir_float_16_32_32\
	src/solve_mpir_posit_16_32_32\
	src/solve_mpir_takum_16_32_32\
	src/solve_mpir_float_16_32_64\
	src/solve_mpir_posit_16_32_64\
	src/solve_mpir_takum_16_32_64\
	src/solve_qr\

GENERATOR =\
	src/generate_full_test_matrices\
	src/generate_sparse_test_matrices\

all: $(EXPERIMENT:=.output_sorted)

src/generate_sparse_test_matrices.output: src/generate_sparse_test_matrices.jl src/TestMatrices.jl config.mk Makefile
src/generate_full_test_matrices.output: src/generate_full_test_matrices.jl src/TestMatrices.jl config.mk Makefile

src/solve_gmres_ilu.output: src/solve_gmres_ilu.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output src/generate_full_test_matrices.output config.mk Makefile
src/solve_lu.output: src/solve_lu.jl src/Experiments.jl src/Float128Conversions.jl src/LU.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output src/generate_full_test_matrices.output config.mk Makefile
src/solve_mpir_float_08_16_32.output: src/solve_mpir_float_08_16_32.jl src/Experiments.jl src/Float128Conversions.jl src/LU.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output src/generate_full_test_matrices.output config.mk Makefile
src/solve_mpir_posit_08_16_32.output: src/solve_mpir_posit_08_16_32.jl src/Experiments.jl src/Float128Conversions.jl src/LU.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output src/generate_full_test_matrices.output config.mk Makefile
src/solve_mpir_takum_08_16_32.output: src/solve_mpir_takum_08_16_32.jl src/Experiments.jl src/Float128Conversions.jl src/LU.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output src/generate_full_test_matrices.output config.mk Makefile
src/solve_mpir_float_16_16_32.output: src/solve_mpir_float_16_16_32.jl src/Experiments.jl src/Float128Conversions.jl src/LU.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output src/generate_full_test_matrices.output config.mk Makefile
src/solve_mpir_posit_16_16_32.output: src/solve_mpir_posit_16_16_32.jl src/Experiments.jl src/Float128Conversions.jl src/LU.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output src/generate_full_test_matrices.output config.mk Makefile
src/solve_mpir_takum_16_16_32.output: src/solve_mpir_takum_16_16_32.jl src/Experiments.jl src/Float128Conversions.jl src/LU.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output src/generate_full_test_matrices.output config.mk Makefile
src/solve_mpir_float_16_32_32.output: src/solve_mpir_float_16_32_32.jl src/Experiments.jl src/Float128Conversions.jl src/LU.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output src/generate_full_test_matrices.output config.mk Makefile
src/solve_mpir_posit_16_32_32.output: src/solve_mpir_posit_16_32_32.jl src/Experiments.jl src/Float128Conversions.jl src/LU.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output src/generate_full_test_matrices.output config.mk Makefile
src/solve_mpir_takum_16_32_32.output: src/solve_mpir_takum_16_32_32.jl src/Experiments.jl src/Float128Conversions.jl src/LU.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output src/generate_full_test_matrices.output config.mk Makefile
src/solve_mpir_float_16_32_64.output: src/solve_mpir_float_16_32_64.jl src/Experiments.jl src/Float128Conversions.jl src/LU.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output src/generate_full_test_matrices.output config.mk Makefile
src/solve_mpir_posit_16_32_64.output: src/solve_mpir_posit_16_32_64.jl src/Experiments.jl src/Float128Conversions.jl src/LU.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output src/generate_full_test_matrices.output config.mk Makefile
src/solve_mpir_takum_16_32_64.output: src/solve_mpir_takum_16_32_64.jl src/Experiments.jl src/Float128Conversions.jl src/LU.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output src/generate_full_test_matrices.output config.mk Makefile
src/solve_qr.output: src/solve_qr.jl src/Experiments.jl src/Float128Conversions.jl src/QR.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output src/generate_full_test_matrices.output config.mk Makefile

src/solve_gmres_ilu.output_sorted: src/solve_gmres_ilu.output src/sort_csv.jl
src/solve_lu.output_sorted: src/solve_lu.output src/sort_csv.jl
src/solve_mpir_float_08_16_32.output_sorted: src/solve_mpir_float_08_16_32.output src/sort_csv.jl
src/solve_mpir_posit_08_16_32.output_sorted: src/solve_mpir_posit_08_16_32.output src/sort_csv.jl
src/solve_mpir_takum_08_16_32.output_sorted: src/solve_mpir_takum_08_16_32.output src/sort_csv.jl
src/solve_mpir_float_16_16_32.output_sorted: src/solve_mpir_float_16_16_32.output src/sort_csv.jl
src/solve_mpir_posit_16_16_32.output_sorted: src/solve_mpir_posit_16_16_32.output src/sort_csv.jl
src/solve_mpir_takum_16_16_32.output_sorted: src/solve_mpir_takum_16_16_32.output src/sort_csv.jl
src/solve_mpir_float_16_32_64.output_sorted: src/solve_mpir_float_16_32_64.output src/sort_csv.jl
src/solve_mpir_posit_16_32_64.output_sorted: src/solve_mpir_posit_16_32_64.output src/sort_csv.jl
src/solve_mpir_takum_16_32_64.output_sorted: src/solve_mpir_takum_16_32_64.output src/sort_csv.jl
src/solve_mpir_float_16_32_32.output_sorted: src/solve_mpir_float_16_32_32.output src/sort_csv.jl
src/solve_mpir_posit_16_32_32.output_sorted: src/solve_mpir_posit_16_32_32.output src/sort_csv.jl
src/solve_mpir_takum_16_32_32.output_sorted: src/solve_mpir_takum_16_32_32.output src/sort_csv.jl
src/solve_qr.output_sorted: src/solve_qr.output src/sort_csv.jl

.jl.format:
	@# work around JuliaFormatter not supporting tabs for indentation
	@# by unexpanding a very wide 16-blank-indent
	$(JULIA) $(JULIA_FLAGS) -- "src/format.jl" "$<" && unexpand -t 16 "$<" > "$<.temp" && mv -f "$<.temp" "$<" && touch "$@"

.jl.output:
	@# experiments print a list of output files, store it an output witness
	$(JULIA) $(JULIA_FLAGS) -- "$<" $(JULIA_SCRIPT_FLAGS) > "$@.temp" && mv -f "$@.temp" "$@"

.output.output_sorted:
	@# use the output witness files and process each .csv file contained
	@# into a .sorted.csv file, outputting another witness file
	@# (.output_sorted) containing the file names
	$(JULIA) $(JULIA_FLAGS) -- "src/sort_csv.jl" "$<" > "$@.temp" && mv -f "$@.temp" "$@"

clean:
	@# use the output witnesses to clean up the output files, except
	@# those from the generators
	for w in $(EXPERIMENT:=.output); do if [ -f "$$w" ]; then xargs rm -f < "$$w"; fi; done
	for w in $(EXPERIMENT:=.output_sorted); do if [ -f "$$w" ]; then xargs rm -f < "$$w"; fi; done
	for d in $(EXPERIMENT); do if [ -d "`basename "$$d"`" ]; then rmdir "out/`basename "$$d"`"; fi; done
	rm -f $(EXPERIMENT:=.output) $(EXPERIMENT:=.output.temp) $(EXPERIMENT:=.output_sorted) $(EXPERIMENT:=.output_sorted.temp)
	rm -f $(COMMON:=.format) $(EXPERIMENT:=.format)

clean-generated:
	@# remove the generated files using the output witnesses
	for w in $(GENERATOR:=.output); do if [ -f "$$w" ]; then xargs rm -f < "$$w"; fi; done

format: $(COMMON:=.format) $(EXPERIMENT:=.format) $(GENERATOR:=.format)

.PHONY: all clean format
