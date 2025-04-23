#!/bin/sh

MAXIMUM_FILESIZE="500K"

# print the to-be-generated file paths as a witness to standard output before we start
printf "out/graph_all_ids\nout/graph_files\n"

# get the network page and scrape the contained graph download links.
# Remove the common URL preffix and .zip suffix, yielding a list of graph IDs of the
# form category/name. We only do that if the out/graph_all_ids file does not
# exist
if [ ! -e "out/graph_all_ids" ]; then
	curl -sS https://networkrepository.com/networks.php | grep -oE 'https://nrvis.com/./download/data/[^"]+' | sed -e 's|^https://nrvis.com/./download/data/||' -e 's|\.zip$||' > out/graph_all_ids
fi

# create the output directory
mkdir -p "out/graphs"

# do not clear the out/graph_files file, as we only append further

# loop over the graph IDs
while IFS= read -r GRAPH_ID; do
	# get the graph category and graph name from the graph ID, abusing
	# dirname and basename a bit
	GRAPH_CATEGORY=$(dirname "$GRAPH_ID")
	GRAPH_NAME=$(basename "$GRAPH_ID")

	# check if the graph file already exists (the suffix can be anything,
	# like .mtx or .egde, so we match for any file), which prompts
	# us to skip
	if ls "out/graphs/$GRAPH_CATEGORY/$GRAPH_NAME".* >/dev/null 2>&1; then
		continue
	fi

	# check if the graph is blacklisted (malformed, complex, 0-indexing)
	BLACKLIST="\
		dynamic/ia-hospital-ward-proximity-attr \
		ia/ia-hospital-ward-proximity-attr \
		labeled/escorts \
		misc/inf-contiguous-usa \
		misc/cavity15 \
		misc/dwg961a \
		misc/dwg961b \
		misc/eco-everglades \
		misc/eco-florida \
		misc/eco-mangwet \
		misc/eco-stmarks \
		misc/M80PI_n \
		misc/mhd1280a \
		misc/mhd1280b \
		misc/mplate \
		misc/qc324 \
		rand/ba_1k_100k \
		rand/ba_1k_10k \
		rand/ba_1k_150k \
		rand/ba_1k_15k \
		rand/ba_1k_20k \
		rand/ba_1k_25k \
		rand/ba_1k_2k \
		rand/ba_1k_30k \
		rand/ba_1k_40k \
		rand/ba_1k_4k \
		rand/ba_1k_50k \
		rand/ba_1k_6k \
		rand/ba_1k_8k \
		rand/er_graph_1k_100k \
		rand/er_graph_1k_10k \
		rand/er_graph_1k_12k \
		rand/er_graph_1k_14k \
		rand/er_graph_1k_200k \
		rand/er_graph_1k_20k \
		rand/er_graph_1k_25k \
		rand/er_graph_1k_30k \
		rand/er_graph_1k_35k \
		rand/er_graph_1k_40k \
		rand/er_graph_1k_4k \
		rand/er_graph_1k_50k \
		rand/er_graph_1k_60k \
		rand/er_graph_1k_6k \
		rand/er_graph_1k_8k \
		rand/geo1k_100k \
		rand/geo1k_10k \
		rand/geo1k_12k \
		rand/geo1k_14k \
		rand/geo1k_150k \
		rand/geo1k_20k \
		rand/geo1k_25k \
		rand/geo1k_30k \
		rand/geo1k_35k \
		rand/geo1k_40k \
		rand/geo1k_4k \
		rand/geo1k_50k \
		rand/geo1k_6k \
		rand/geo1k_8k \
		rec/rec-movielens-user-tag-10m \
		soc/soc-firm-hi-tech \
		tscc/scc_rt_islam \
		"

	BLACKLISTED=false
	for BLACKLISTED_GRAPH_ID in $BLACKLIST; do
		if [ "$GRAPH_ID" = "$BLACKLISTED_GRAPH_ID" ]; then
			BLACKLISTED=true
		fi
	done
	if [ "$BLACKLISTED" = true ]; then
		continue
	fi

	# generate the download URL from the graph ID
	GRAPH_URL="https://nrvis.com/download/data/$GRAPH_CATEGORY/$GRAPH_NAME.zip"

	# generate a temporary file suffix
	TMPFILE=$(mktemp)

	# fetch the graph zip file into the temporary file path, enforcing a
	# size limit. If we breach the limit, we simply remove the temporary
	# file and continue. We also fail silently on HTTP errors, as some
	# zip files don't exist.
	# This approach saves us from sending two requests per file
	# (a HEAD request to get the size and then the download itself).
	# It's a bit greedy, but makes more sense.
	if curl -fs --max-filesize $MAXIMUM_FILESIZE -o "$TMPFILE" "$GRAPH_URL"; then
		# the file size is within the limit

		# first check if the bloody zip file is even well-formed,
		# given some in the archive are NOT!
		if ! unzip -t "$TMPFILE" >/dev/null 2>&1; then
			# clean up and continue, not even mentioning
			# this crap
			rm -f "$TMPFILE"
			continue
		fi

		# the zip file contains a readme.html and the graph file
		# itself, named like the GRAPH_NAME and an arbitrary suffix
		# (mtx, edgelist, etc.). We only want the latter, so we
		# list the zip contents, obtain the filename and specifically
		# unpack it. This is a bit hacky given unzip(1) sucks.
		GRAPH_FILE=$(unzip -l "$TMPFILE" | awk '{print $NF}' | grep "$GRAPH_NAME\." | head -n 1)

		if [ -z "$GRAPH_FILE" ]; then
			# print error and clean up
			printf "Error: The archive $GRAPH_ID does not contain any file with basename $GRAPH_NAME\n" >&2
			rm -f "$TMPFILE"
			exit 1
		fi

		# create a directory for the graph category
		mkdir -p "out/graphs/$GRAPH_CATEGORY"

		# extract the graph file
		unzip -j "$TMPFILE" "$GRAPH_FILE" -d "out/graphs/$GRAPH_CATEGORY" >/dev/null 2>&1

		# get the graph path
		GRAPH_PATH="out/graphs/$GRAPH_CATEGORY/$(basename $GRAPH_FILE)"

		# set the file permissions on the graph file
		chmod 644 "$GRAPH_PATH"

		# apply some general post processing
		case "$GRAPH_FILE" in
		*.mtx)
			# ensure two percent signs before MatrixMarket, as the
			# Julia MatrixMarket package is pretty strict about it
			sed -i '1s/^%MatrixMarket/%%MatrixMarket/' "$GRAPH_PATH"

			# if a line begins with "%[0-9]" or "% [0-9]", remove the
			# comment, as then this file has commented out the matrix
			# dimensions, which is wrong.
			sed -i 's/^% \?\([0-9]\)/\1/' "$GRAPH_PATH"
			;;
		esac

		# apply some file-specific post processing
		case "$GRAPH_FILE" in
		*DSJC500-5.mtx)
			# the entry count is doubled
			sed -i '2s/^500 500 125248$/500 500 62624/' "$GRAPH_PATH"
			;;
		*p-hat500-1.mtx)
			# the dimensions are messed up
			sed -i '2s/^  500$/500 500 31569/' "$GRAPH_PATH"
			;;
		*mark3jac140.mtx)
			# the dimensions are messed up
			sed -i '2s/^64089 64089 399735$/4557 4557 19848/' "$GRAPH_PATH"
			;;
		esac

		# clean up the temporary zip file
		rm -f "$TMPFILE"

		# print the graph path to standard output as a witness
		# we take the basename as GRAPH_FILE is a path, but unzip(1)
		# with -j only yields the file in the end
		printf "$GRAPH_PATH\n"

		# also append the graph path to out/graph_files
		printf "$GRAPH_PATH\n" >> "out/graph_files"
	else
		# the file is too large, remove the temporary file and
		# skip it
		rm -f "$TMPFILE"
		continue
	fi
done < out/graph_all_ids
