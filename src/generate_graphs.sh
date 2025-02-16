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
		unzip -j "$TMPFILE" "$GRAPH_FILE" -d "out/graphs/$GRAPH_CATEGORY"

		# clean up the temporary zip file
		rm -f "$TMPFILE"

		# print the graph path to standard output as a witness
		# we take the basename as GRAPH_FILE is a path, but unzip(1)
		# with -j only yields the file in the end
		printf "out/graphs/$GRAPH_CATEGORY/$(basename $GRAPH_FILE)\n"

		# also append the graph path to out/graph_files
		printf "out/graphs/$GRAPH_CATEGORY/$(basename $GRAPH_FILE)\n" >> "out/graph_files"	
	else
		# the file is too large, remove the temporary file and
		# skip it
		rm -f "$TMPFILE"
		continue
	fi
done < out/graph_all_ids
