FROM julia

# Install make(1) and LaTeX TeX Live with latexmk(1) and clean up afterwards
RUN apt-get update && apt-get install -y --no-install-recommends make texlive-full latexmk && rm -rf /var/lib/apt/lists/*

# Set the shell prompt name to clearly indicate it as the docker shell
ENV PS1="[MuFoLAB docker] \$ "
