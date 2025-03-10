# See LICENSE file for copyright and license details
# MuFoLAB - Multi-Format Linear Algebra Benchmarks
.POSIX:
.SUFFIXES:
.SUFFIXES: .format .jl .output .output_sorted .sh

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
	src/eigen_graph_biological_08\
	src/eigen_graph_biological_16\
	src/eigen_graph_biological_32\
	src/eigen_graph_biological_64\
	src/eigen_graph_social_08\
	src/eigen_graph_social_16\
	src/eigen_graph_social_32\
	src/eigen_graph_social_64\
	src/eigen_graph_infrastructure_08\
	src/eigen_graph_infrastructure_16\
	src/eigen_graph_infrastructure_32\
	src/eigen_graph_infrastructure_64\
	src/eigen_graph_misc_08\
	src/eigen_graph_misc_16\
	src/eigen_graph_misc_32\
	src/eigen_graph_misc_64\
	src/eigen_smallest_08\
	src/eigen_smallest_16\
	src/eigen_smallest_32\
	src/eigen_smallest_64\
	src/conversion\
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
	src/generate_graph_test_matrices\
	src/generate_graphs\
	src/generate_sparse_test_matrices\
	src/generate_stochastic_test_matrices\
	src/sparse_matrix_condition_numbers\

all: $(EXPERIMENT:=.output_sorted)

src/generate_full_test_matrices.output: src/generate_full_test_matrices.jl src/TestMatrices.jl config.mk Makefile
src/generate_graphs.output: src/generate_graphs.sh config.mk Makefile
src/generate_graph_test_matrices.output: src/generate_graph_test_matrices.jl src/TestMatrices.jl src/generate_graphs.output config.mk Makefile
src/generate_sparse_test_matrices.output: src/generate_sparse_test_matrices.jl src/TestMatrices.jl config.mk Makefile
src/generate_stochastic_test_matrices.output: src/generate_stochastic_test_matrices.jl src/TestMatrices.jl config.mk Makefile
src/sparse_matrix_condition_numbers.output: src/sparse_matrix_condition_numbers.jl src/TestMatrices.jl config.mk Makefile

src/eigen_graph_biological_08.output: src/eigen_graph_biological_08.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_graph_test_matrices.output config.mk Makefile
src/eigen_graph_biological_16.output: src/eigen_graph_biological_16.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_graph_test_matrices.output config.mk Makefile
src/eigen_graph_biological_32.output: src/eigen_graph_biological_32.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_graph_test_matrices.output config.mk Makefile
src/eigen_graph_biological_64.output: src/eigen_graph_biological_64.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_graph_test_matrices.output config.mk Makefile
src/eigen_graph_social_08.output: src/eigen_graph_social_08.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_graph_test_matrices.output config.mk Makefile
src/eigen_graph_social_16.output: src/eigen_graph_social_16.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_graph_test_matrices.output config.mk Makefile
src/eigen_graph_social_32.output: src/eigen_graph_social_32.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_graph_test_matrices.output config.mk Makefile
src/eigen_graph_social_64.output: src/eigen_graph_social_64.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_graph_test_matrices.output config.mk Makefile
src/eigen_graph_infrastructure_08.output: src/eigen_graph_infrastructure_08.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_graph_test_matrices.output config.mk Makefile
src/eigen_graph_infrastructure_16.output: src/eigen_graph_infrastructure_16.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_graph_test_matrices.output config.mk Makefile
src/eigen_graph_infrastructure_32.output: src/eigen_graph_infrastructure_32.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_graph_test_matrices.output config.mk Makefile
src/eigen_graph_infrastructure_64.output: src/eigen_graph_infrastructure_64.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_graph_test_matrices.output config.mk Makefile
src/eigen_graph_misc_08.output: src/eigen_graph_misc_08.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_graph_test_matrices.output config.mk Makefile
src/eigen_graph_misc_16.output: src/eigen_graph_misc_16.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_graph_test_matrices.output config.mk Makefile
src/eigen_graph_misc_32.output: src/eigen_graph_misc_32.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_graph_test_matrices.output config.mk Makefile
src/eigen_graph_misc_64.output: src/eigen_graph_misc_64.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_graph_test_matrices.output config.mk Makefile
src/eigen_smallest_08.output: src/eigen_smallest_08.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output src/generate_full_test_matrices.output config.mk Makefile
src/eigen_smallest_16.output: src/eigen_smallest_16.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output src/generate_full_test_matrices.output config.mk Makefile
src/eigen_smallest_32.output: src/eigen_smallest_32.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output src/generate_full_test_matrices.output config.mk Makefile
src/eigen_smallest_64.output: src/eigen_smallest_64.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output src/generate_full_test_matrices.output config.mk Makefile
src/conversion.output: src/conversion.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output config.mk Makefile
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

src/eigen_graph_biological_08.output_sorted: src/eigen_graph_biological_08.output src/sort_csv.jl
src/eigen_graph_biological_16.output_sorted: src/eigen_graph_biological_16.output src/sort_csv.jl
src/eigen_graph_biological_32.output_sorted: src/eigen_graph_biological_32.output src/sort_csv.jl
src/eigen_graph_biological_64.output_sorted: src/eigen_graph_biological_64.output src/sort_csv.jl
src/eigen_graph_social_08.output_sorted: src/eigen_graph_social_08.output src/sort_csv.jl
src/eigen_graph_social_16.output_sorted: src/eigen_graph_social_16.output src/sort_csv.jl
src/eigen_graph_social_32.output_sorted: src/eigen_graph_social_32.output src/sort_csv.jl
src/eigen_graph_social_64.output_sorted: src/eigen_graph_social_64.output src/sort_csv.jl
src/eigen_graph_infrastructure_08.output_sorted: src/eigen_graph_infrastructure_08.output src/sort_csv.jl
src/eigen_graph_infrastructure_16.output_sorted: src/eigen_graph_infrastructure_16.output src/sort_csv.jl
src/eigen_graph_infrastructure_32.output_sorted: src/eigen_graph_infrastructure_32.output src/sort_csv.jl
src/eigen_graph_infrastructure_64.output_sorted: src/eigen_graph_infrastructure_64.output src/sort_csv.jl
src/eigen_graph_misc_08.output_sorted: src/eigen_graph_misc_08.output src/sort_csv.jl
src/eigen_graph_misc_16.output_sorted: src/eigen_graph_misc_16.output src/sort_csv.jl
src/eigen_graph_misc_32.output_sorted: src/eigen_graph_misc_32.output src/sort_csv.jl
src/eigen_graph_misc_64.output_sorted: src/eigen_graph_misc_64.output src/sort_csv.jl
src/eigen_smallest_08.output_sorted: src/eigen_smallest_08.output src/sort_csv.jl
src/eigen_smallest_16.output_sorted: src/eigen_smallest_16.output src/sort_csv.jl
src/eigen_smallest_32.output_sorted: src/eigen_smallest_32.output src/sort_csv.jl
src/eigen_smallest_64.output_sorted: src/eigen_smallest_64.output src/sort_csv.jl
src/conversion.output_sorted: src/conversion.output src/sort_csv.jl
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

.sh.format:
	@# no-op

.jl.output:
	@# experiments print a list of output files, store it an output witness
	$(JULIA) $(JULIA_FLAGS) -- "$<" $(JULIA_SCRIPT_FLAGS) > "$@.temp" && mv -f "$@.temp" "$@"

.sh.output:
	@# shell scripts print a list of output files, store it an output witness
	$(SH) "$<" > "$@.temp" && mv -f "$@.temp" "$@"

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
