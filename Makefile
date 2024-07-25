# See LICENSE file for copyright and license details
# Takum Linear Algebra Benchmarks
.POSIX:
.SUFFIXES:

include config.mk

COMMON =\
	src/Utilities\

EXPERIMENT =\
	src/solve_direct\

	#src/solve_mixed_iterative_refinement\
	#src/solve_squeeze\

all: $(EXPERIMENT:=.w)

src/solve_direct.w: src/solve_direct.jl src/Utilities.jl

# each experiment program prints the files it generated to stdout, we store
# their names in a witness file
$(EXPERIMENT:=.w):
	$(JULIA) $(JULIAFLAGS) $(EXPERIMENT:=.jl) > $(EXPERIMENT:=.w)

# go over each witness file, and if it exists, remove all files listed in it.
# Then remove the witness file.
clean:
	for w in $(EXPERIMENT:=.w); do if [ -f "$$w" ]; then xargs rm -f < "$$w"; fi; done
	rm -f $(EXPERIMENT:=.w)

# use JuliaFormatter to automatically format the code. Given JuliaFormatter
# does not support tabs for indentation we let it run with an obnoxiously
# large blank-indentation value of 16 to then automatically detect and
# convert them to tabs using unexpand(1). We set the value so high to
# avoid converting blanks properly used for alignment.
# Along the way we increase the margin using a 1-level-indent as the
# reference case based on a soft target of 85 characters per row. Here
# we assume the default case of a tab being 8 blanks wide
format:
	for s in $(COMMON:=.jl) $(EXPERIMENT:=.jl); do printf "using JuliaFormatter; format_file(\"$$s\"; indent = 16, margin = 85 + 2 * (16 - 8), always_for_in = true, whitespace_typedefs = true, whitespace_ops_in_indices = true, remove_extra_newlines = true, pipe_to_function_call = true, short_to_long_function_def = true, always_use_return = true, align_struct_field = true, align_conditional = true, align_assignment = true, align_pair_arrow = true, align_matrix = true, conditional_to_if = true, normalize_line_endings = \"unix\", trailing_comma = true, indent_submodule = true, separate_kwargs_with_semicolon = true, short_circuit_to_if = true); exit(0);\n" | $(JULIA) $(JULIAFLAGS); unexpand -t 16 "$$s" > "$$s.tmp"; mv -f "$$s.tmp" "$$s"; done

.PHONY: all clean format
