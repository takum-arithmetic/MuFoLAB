#!/bin/sh

# Build docker image
docker build -t "mufolab" .

# Run docker container, mounting the code folder and using it as a working directory.
docker run -it --rm -v "$PWD":/mnt/mufolab -w "/mnt/mufolab" --entrypoint sh mufolab -c "julia --project -e 'using Pkg; Pkg.instantiate()'; exec sh"
