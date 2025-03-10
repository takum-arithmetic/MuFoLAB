# See LICENSE file for copyright and license details.
using CSV
using DataFrames

@kwdef struct ReplacementValues
	zero::Float64
	negative_infinity::Float64
	positive_infinity::Float64
end

function ReplacementValues(type_name::String, df::DataFrame)

	# A ratio of the decadic dynamic range of the value group, how
	# far should the replacement values, each, be from the 'true'
	# dataset in the plot?
	plot_distance_ratio = 0.15

	# We check in which group the type_name is
	groups = [
		["Float8_4", "Float8_5", "Posit8", "Takum8", "LinearTakum8"],
		["BFloat16", "Float16", "Posit16", "Takum16", "LinearTakum16"],
		["Float32", "Posit32", "Takum32", "LinearTakum32"],
		["Float64", "Posit64", "Takum64", "LinearTakum64"],
	]

	type_group = nothing
	for group in groups
		if type_name in group
			type_group = group
		end
	end
	if type_group === nothing
		throw(ErrorException("Type '$(type_name)' matches no group"))
	end

	# Filter out those type names in the type_group which are not
	# present in the DataFrame
	filter!(c -> c in names(df), type_group)

	# Extract the type columns belonging to the type group from the
	# DataFrame as a matrix
	type_group_columns = Matrix(df[!, type_group])

	# Remove all zero and non-finite values and write the rest in a vector
	type_group_values = type_group_columns[isfinite.(
		type_group_columns
	) .&& .!iszero.(type_group_columns)]

	# Return early if the collection is empty
	if isempty(type_group_values)
		return ReplacementValues(;
			zero = 0.0,
			negative_infinity = 0.0,
			positive_infinity = 0.0,
		)
	end

	# Sort the vector
	type_group_values = sort(type_group_values)

	# Get the decadic dynamic range of the values
	dynamic_range = log10(maximum(type_group_values) / minimum(type_group_values))

	# The plot distance is the ratio applied to the full dynamic range
	plot_distance = plot_distance_ratio * dynamic_range

	# We obtain the replacement values by looking at the decadic
	# logarithms of the value maxima and minima, adding the distance
	# and rounding up or down to the next full integer.
	return ReplacementValues(;
		zero = exp10(
			floor(
				log10(minimum(type_group_values)) -
				plot_distance,
			),
		),
		negative_infinity = exp10(
			ceil(
				log10(maximum(type_group_values)) +
				plot_distance,
			),
		),
		positive_infinity = exp10(
			ceil(
				log10(maximum(type_group_values)) +
				2 * plot_distance,
			),
		),
	)
end

function sort_csv(input_file_name::String)
	# Read the CSV as a DataFrame
	df = CSV.read(input_file_name, DataFrame)

	# Discard the first column containing the matrix names
	df = df[:, 2:end]

	# Replace the values 0.0 (not directly plottable in a logarithmic
	# plot), -Inf (indicating singularity) and Inf
	# (indicating under-/overflow) with replacement values for easier
	# plotting.
	type_names = names(df)

	# We make a copy of df where the replacement values are not yet
	# entered
	df_original = copy(df)

	for i in 1:ncol(df)
		if contains(input_file_name, "mpir")
			replacement_values = ReplacementValues(;
				zero = 1e-1,
				negative_infinity = 1e3,
				positive_infinity = 1e4,
			)
		else
			replacement_values = ReplacementValues(
				type_names[i],
				df_original,
			)
		end

		df[(df[:, i] .== 0.0), i]  .= replacement_values.zero
		df[(df[:, i] .== Inf), i]  .= replacement_values.positive_infinity
		df[(df[:, i] .== -Inf), i] .= replacement_values.negative_infinity
	end

	# Sort each column
	for i in 1:ncol(df)
		df[:, i] = sort(df[:, i])
	end

	# Add a column to the left of the DataFrame containing the percent
	# values
	insertcols!(df, 1, :percent => (0:(nrow(df) - 1)) ./ (nrow(df) - 1))

	# Write DataFrame as a CSV with the input file name, but with
	# the .sorted.csv extension
	output_file_name = chopsuffix(input_file_name, ".csv") * ".sorted.csv"

	CSV.write(output_file_name, df)

	# Write the name of the output file to standard output
	return println(output_file_name)
end

begin
	if length(ARGS) != 1
		throw(ArgumentError("This program takes exactly one argument"))
	end

	for file_name in eachline(ARGS[1])
		sort_csv(file_name)
	end
end
