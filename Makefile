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

EXPERIMENT_EIGEN =\
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
	src/eigen_general_08\
	src/eigen_general_16\
	src/eigen_general_32\
	src/eigen_general_64\
	src/eigen_properties_graph_biological\
	src/eigen_properties_graph_social\
	src/eigen_properties_graph_infrastructure\
	src/eigen_properties_graph_misc\
	src/eigen_properties_general\

EXPERIMENT_SOLVE =\
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


EXPERIMENT =\
	src/convert\
	$(EXPERIMENT_EIGEN)\
	$(EXPERIMENT_SOLVE)\

GENERATOR =\
	src/generate_full_test_matrices\
	src/generate_graph_test_matrices\
	src/generate_graphs\
	src/generate_sparse_test_matrices\
	src/generate_stochastic_test_matrices\
	src/sparse_matrix_condition_numbers\

EIGEN_PLOTS =\
	plots/eigen_general/eigen_general\
	plots/eigen_graph_biological/eigen_graph_biological\
	plots/eigen_graph_infrastructure/eigen_graph_infrastructure\
	plots/eigen_graph_social/eigen_graph_social\
	plots/eigen_graph_misc/eigen_graph_misc\
	plots/eigen_properties/eigen_properties\

all: $(EXPERIMENT:=.output_sorted)
eigen: $(EXPERIMENT_EIGEN:=.output_sorted) plots/eigen.pdf
solve: $(EXPERIMENT_SOLVE:=.output_sorted)

#src/generate_full_test_matrices.output: src/generate_full_test_matrices.jl src/TestMatrices.jl config.mk Makefile
#src/generate_graphs.output: src/generate_graphs.sh config.mk Makefile
#src/generate_graph_test_matrices.output: src/generate_graph_test_matrices.jl src/TestMatrices.jl src/generate_graphs.output config.mk Makefile
#src/generate_sparse_test_matrices.output: src/generate_sparse_test_matrices.jl src/TestMatrices.jl config.mk Makefile
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
src/eigen_general_08.output: src/eigen_general_08.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output config.mk Makefile
src/eigen_general_16.output: src/eigen_general_16.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output config.mk Makefile
src/eigen_general_32.output: src/eigen_general_32.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output config.mk Makefile
src/eigen_general_64.output: src/eigen_general_64.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output config.mk Makefile
src/eigen_properties_general.output: src/eigen_properties_general.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output config.mk Makefile
src/eigen_properties_graph_biological.output: src/eigen_properties_graph_biological.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_graph_test_matrices.output config.mk Makefile
src/eigen_properties_graph_social.output: src/eigen_properties_graph_social.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_graph_test_matrices.output config.mk Makefile
src/eigen_properties_graph_infrastructure.output: src/eigen_properties_graph_infrastructure.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_graph_test_matrices.output config.mk Makefile
src/eigen_properties_graph_misc.output: src/eigen_properties_graph_misc.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_graph_test_matrices.output config.mk Makefile
src/convert.output: src/convert.jl src/Experiments.jl src/Float128Conversions.jl src/TestMatrices.jl src/generate_sparse_test_matrices.output config.mk Makefile
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
src/eigen_general_08.output_sorted: src/eigen_general_08.output src/sort_csv.jl
src/eigen_general_16.output_sorted: src/eigen_general_16.output src/sort_csv.jl
src/eigen_general_32.output_sorted: src/eigen_general_32.output src/sort_csv.jl
src/eigen_general_64.output_sorted: src/eigen_general_64.output src/sort_csv.jl
src/eigen_properties_general.output_sorted: src/eigen_properties_general.output src/sort_csv.jl
src/eigen_properties_graph_biological.output_sorted: src/eigen_properties_graph_biological.output src/sort_csv.jl
src/eigen_properties_graph_social.output_sorted: src/eigen_properties_graph_social.output src/sort_csv.jl
src/eigen_properties_graph_infrastructure.output_sorted: src/eigen_properties_graph_infrastructure.output src/sort_csv.jl
src/eigen_properties_graph_misc.output_sorted: src/eigen_properties_graph_misc.output src/sort_csv.jl
src/convert.output_sorted: src/convert.output src/sort_csv.jl
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

plots/eigen_general/eigen_general.pdf: plots/eigen_general/eigen_general.tex src/eigen_general_08.output_sorted src/eigen_general_16.output_sorted src/eigen_general_32.output_sorted src/eigen_general_64.output_sorted
plots/eigen_graph_biological/eigen_graph_biological.pdf: plots/eigen_graph_biological/eigen_graph_biological.tex src/eigen_graph_biological_08.output_sorted src/eigen_graph_biological_16.output_sorted src/eigen_graph_biological_32.output_sorted src/eigen_graph_biological_64.output_sorted
plots/eigen_graph_infrastructure/eigen_graph_infrastructure.pdf: plots/eigen_graph_infrastructure/eigen_graph_infrastructure.tex src/eigen_graph_infrastructure_08.output_sorted src/eigen_graph_infrastructure_16.output_sorted src/eigen_graph_infrastructure_32.output_sorted src/eigen_graph_infrastructure_64.output_sorted
plots/eigen_graph_social/eigen_graph_social.pdf: plots/eigen_graph_social/eigen_graph_social.tex src/eigen_graph_social_08.output_sorted src/eigen_graph_social_16.output_sorted src/eigen_graph_social_32.output_sorted src/eigen_graph_social_64.output_sorted
plots/eigen_graph_misc/eigen_graph_misc.pdf: plots/eigen_graph_misc/eigen_graph_misc.tex src/eigen_graph_misc_08.output_sorted src/eigen_graph_misc_16.output_sorted src/eigen_graph_misc_32.output_sorted src/eigen_graph_misc_64.output_sorted
plots/eigen_properties/eigen_properties.pdf: plots/eigen_properties/eigen_properties.tex src/eigen_properties_general.output_sorted src/eigen_properties_graph_biological.output_sorted src/eigen_properties_graph_infrastructure.output_sorted src/eigen_properties_graph_social.output_sorted src/eigen_properties_graph_misc.output_sorted

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

$(EIGEN_PLOTS:=.pdf):
	$(LATEXMK) -pdf -cd -shell-escape $(@:.pdf=.tex)

plots/eigen.pdf: plots/eigen.tex $(EIGEN_PLOTS:=.pdf)
	$(LATEXMK) -pdf -cd plots/eigen.tex

clean:
	@# use the output witnesses to clean up the output files, except
	@# those from the generators
	for w in $(EXPERIMENT:=.output); do if [ -f "$$w" ]; then xargs rm -f < "$$w"; fi; done
	for w in $(EXPERIMENT:=.output_sorted); do if [ -f "$$w" ]; then xargs rm -f < "$$w"; fi; done
	for d in $(EXPERIMENT); do if [ -d "`basename "$$d"`" ]; then rmdir "out/`basename "$$d"`"; fi; done
	rm -f $(EXPERIMENT:=.output) $(EXPERIMENT:=.output.temp) $(EXPERIMENT:=.output_sorted) $(EXPERIMENT:=.output_sorted.temp)
	rm -f $(COMMON:=.format) $(EXPERIMENT:=.format)
	rm -f plots/eigen_*/08* plots/eigen_*/16* plots/eigen_*/32* plots/eigen_*/64* plots/eigen_*/*.auxlock
	$(LATEXMK) -C -cd $(EIGEN_PLOTS:=.tex) plots/eigen.tex
	@# remove empty folders in out/ recursively
	find out/ -empty -type d -delete

clean-generated:
	@# remove the generated files using the output witnesses
	for w in $(GENERATOR:=.output); do if [ -f "$$w" ]; then xargs rm -f < "$$w"; fi; done
	@# remove empty folders in out/ recursively
	find out/ -empty -type d -delete

format: $(COMMON:=.format) $(EXPERIMENT:=.format) $(GENERATOR:=.format)

.PHONY: all clean format
