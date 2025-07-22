FROM --platform=linux/amd64 ubuntu:20.04

# Set non-interactive mode for apt
ENV DEBIAN_FRONTEND=noninteractive

# Set environment variables
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH:-}"
ENV PYTHONPATH="/opt/orchard/src:${PYTHONPATH:-}"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget curl gcc g++ make perl zlib1g-dev \
    build-essential autoconf automake libtool pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda with Python 3.9
WORKDIR /opt
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_23.5.2-0-Linux-x86_64.sh \
    && bash Miniconda3-py39_23.5.2-0-Linux-x86_64.sh -b -p /opt/anaconda \
    && rm -f Miniconda3-py39_23.5.2-0-Linux-x86_64.sh

ENV PATH="/opt/anaconda/bin:$PATH"

# Initialize conda and create environment
RUN /opt/anaconda/bin/conda init bash \
    && /bin/bash -c "source /opt/anaconda/bin/activate && conda create --name pipeline python=3.9"

# Install conda-forge packages (matching your environment)
RUN /bin/bash -c "source /opt/anaconda/bin/activate pipeline && conda install -c conda-forge \
    numpy=1.26.4 \
    scipy=1.13.1 \
    matplotlib=3.9.4 \
    pandas=2.3.0 \
    astropy=6.0.1 \
    scikit-learn=1.6.1 \
    photutils=1.11.0 \
    fitsio=1.2.6 \
    cfitsio=4.6.2 \
    wcslib=8.2.2 \
    tqdm=4.67.1 \
    pyyaml=6.0.2 \
    requests=2.32.4"

# Install pip packages (matching your environment)
RUN /bin/bash -c "source /opt/anaconda/bin/activate pipeline && pip install \
    twirl==0.4.2 \
    astroquery==0.4.10 \
    scikit-image==0.21.0 \
    tabulate==0.9.0 \
    progressbar2==4.5.0 \
    psutil==7.0.0 \
    python-dotenv==1.1.0 \
    tifffile==2024.8.30 \
    pexpect==4.9.0"

# Set conda environment as default
ENV CONDA_DEFAULT_ENV=pipeline
ENV PATH="/opt/anaconda/envs/pipeline/bin:$PATH"

# Install Autoconf, Automake for casutools
WORKDIR /usr/local/src
RUN wget http://ftp.gnu.org/gnu/autoconf/autoconf-2.69.tar.gz \
    && tar xf autoconf-2.69.tar.gz \
    && cd autoconf-2.69 \
    && ./configure --prefix=/usr \
    && make \
    && make install \
    && cd .. \
    && rm -rf autoconf-2.69 autoconf-2.69.tar.gz \
    && wget http://ftp.gnu.org/gnu/automake/automake-1.16.1.tar.gz \
    && tar xf automake-1.16.1.tar.gz \
    && cd automake-1.16.1 \
    && ./configure --prefix=/usr \
    && make \
    && make install \
    && cd .. \
    && rm -rf automake-1.16.1 automake-1.16.1.tar.gz

# Install cfitsio and wcslib from source
WORKDIR /usr/local/src
RUN wget https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-4.2.0.tar.gz \
    && tar xzf cfitsio-4.2.0.tar.gz \
    && cd cfitsio-4.2.0 \
    && ./configure --prefix=/usr/local \
    && make \
    && make install \
    && cd .. \
    && rm -rf cfitsio-4.2.0* \
    && wget ftp://ftp.atnf.csiro.au/pub/software/wcslib/wcslib.tar.bz2 \
    && tar xjf wcslib.tar.bz2 \
    && cd wcslib-8.4 \
    && ./configure --prefix=/usr/local --with-cfitsio=/usr/local \
    && make \
    && make install \
    && cd .. \
    && rm -rf wcslib-8.4*

# Install casutools 1.0.44pre
WORKDIR /opt
COPY docker/casutools-1.0.44pre.tar.gz ./
RUN tar xzf casutools-1.0.44pre.tar.gz \
    && cd casutools-1.0.44pre \
    && ./configure --prefix=/opt/casutools --with-cfitsio=/usr/local --with-wcs=/usr/local \
    && make \
    && make install \
    && cd .. \
    && rm -rf casutools-1.0.44pre*

# Set environment variables for casutools
ENV PATH="/opt/casutools/bin:$PATH"

# Set up pipeline environment
ENV ORCHARD_PATH="/opt/orchard/src"
WORKDIR /opt/orchard

# Copy pipeline code
COPY src/ ./src/

# Create pipeline data directory structure
RUN mkdir -p /data

# Set working directory to pipeline source
WORKDIR /opt/orchard/src

# Ensure pipeline environment is activated by default
RUN echo "source /opt/anaconda/bin/activate pipeline" >> ~/.bashrc
ENV BASH_ENV=~/.bashrc
SHELL ["/bin/bash", "-c"]

CMD ["/bin/bash"]