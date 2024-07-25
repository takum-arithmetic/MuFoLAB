# See LICENSE file for copyright and license details.
using JuliaFormatter

# We set the indent to 16 to have 16 blanks for each indentation level,
# which is easily detected by unexpand(1). We don't want to set it
# lower as we do not want to unexpand any blanks used for alignment,
# which is desired.
#
# Each line should contain at most 85 characters. We use a heuristic
# for the common case of 2 indentation levels to expand it by the
# difference of the common tab (8 characters) and the set indent value
# such that our large value for indent does not yield short lines.
format_file(
	ARGS[1];
	indent = 16,
	margin = 85 + 2 * (16 - 8),
	always_for_in = true,
	whitespace_typedefs = true,
	whitespace_ops_in_indices = true,
	remove_extra_newlines = true,
	pipe_to_function_call = true,
	short_to_long_function_def = true,
	always_use_return = true,
	align_struct_field = true,
	align_conditional = true,
	align_assignment = true,
	align_pair_arrow = true,
	align_matrix = true,
	conditional_to_if = true,
	normalize_line_endings = "unix",
	trailing_comma = true,
	indent_submodule = true,
	separate_kwargs_with_semicolon = true,
	short_circuit_to_if = true,
)
