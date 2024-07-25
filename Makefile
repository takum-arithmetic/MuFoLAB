# See LICENSE file for copyright and license details
# Takum Linear Algebra Benchmarks
.POSIX:
.SUFFIXES:
.SUFFIXES: .format .jl .output

include config.mk

COMMON =\
	src/format\
	src/Utilities\

EXPERIMENT =\
	src/solve_direct\

	#src/solve_mixed_iterative_refinement\
	#src/solve_squeeze\

src/solve_direct.output: src/solve_direct.jl src/Utilities.jl

all: $(EXPERIMENT:=.output)

# use JuliaFormatter to automatically format the code. Given JuliaFormatter
# does not support tabs for indentation we let it run with an obnoxiously
# large blank-indentation value of 16 to then automatically detect and
# convert them to tabs using unexpand(1). We set the value so high to
# avoid converting blanks properly used for alignment.
# Along the way we increase the margin using a 1-level-indent as the
# reference case based on a soft target of 85 characters per row. Here
# we assume the default case of a tab being 8 blanks wide
.jl.format:
	$(JULIA) $(JULIAFLAGS) "src/format.jl" "$<" && unexpand -t 16 "$<" > "$<.temp" && mv -f "$<.temp" "$<" && touch "$@"

# each experiment program prints the files it generated to stdout; we store
# their names in a temporary output witness file and rename it to the final
# output when the command was successful.
.jl.output:
	$(JULIA) $(JULIAFLAGS) "$<" > "$@.temp" && mv -f "$@.temp" "$@"

# besides removing all witnesses (and possible temporary ones), use the
# output witnesses to clean up all output files
clean:
	for w in $(EXPERIMENT:=.output); do if [ -f "$$w" ]; then xargs rm -f < "$$w"; fi; done
	rm -f $(EXPERIMENT:=.output) $(EXPERIMENT:=.output.temp)
	rm -f $(COMMON:=.format) $(EXPERIMENT:=.format)

format: $(COMMON:=.format) $(EXPERIMENT:=.format)

.PHONY: all clean format
