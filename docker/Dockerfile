FROM --platform=linux/amd64 centos:7.9.2009

# Set environment variables
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH:-}"
ENV PYTHONPATH="/opt/orchard/src:${PYTHONPATH:-}"

# Install EPEL repository first (needed for astrometry.net)
RUN yum install -y epel-release

# Install system dependencies with retry logic
RUN for i in 1 2 3; do \
        yum install -y \
            wget curl gcc gcc-c++ make perl zlib-devel \
            autoconf automake libtool pkgconfig \
            libX11 libXext libXrender libXtst libXft libXpm \
            astrometry.net \
        && break || sleep 5; \
    done && \
    yum clean all

# Install Miniconda with Python 3.9
WORKDIR /opt
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_23.5.2-0-Linux-x86_64.sh \
    && bash Miniconda3-py39_23.5.2-0-Linux-x86_64.sh -b -p /opt/anaconda \
    && rm -f Miniconda3-py39_23.5.2-0-Linux-x86_64.sh

ENV PATH="/opt/anaconda/bin:$PATH"

# Initialize conda and create environment
RUN /opt/anaconda/bin/conda init bash \
    && /bin/bash -c "source /opt/anaconda/bin/activate && conda create --name pipeline python=3.9"

# Copy requirements and install Python packages
COPY docker/requirements.txt /opt/requirements.txt
RUN /bin/bash -c "source /opt/anaconda/bin/activate && conda activate pipeline && pip install numpy==1.26.4"
RUN /bin/bash -c "source /opt/anaconda/bin/activate && conda activate pipeline && pip install --no-cache-dir -r /opt/requirements.txt"

# Set environment variables
ENV CONDA_DEFAULT_ENV=pipeline
ENV PATH="/opt/anaconda/envs/pipeline/bin:$PATH"

# Install Autoconf, Automake, and dependencies
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

# Verify installation
RUN aclocal --version && automake --version && autoconf --version

# Install cfitsio and wcslib
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

# Set library path
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH:-}"

# Install casutools 1.0.43
WORKDIR /opt
COPY docker/casutools-1.0.43.tar.gz ./
RUN tar xzf casutools-1.0.43.tar.gz \
    && cd casutools-1.0.43 \
    && ./configure --prefix=/opt/casutools --with-cfitsio=/usr/local --with-wcs=/usr/local \
    && make \
    && make install \
    && cd .. \
    && rm -rf casutools-1.0.43*

# Set environment variables for casutools
ENV PATH="/opt/casutools/bin:$PATH"

# Create astrometry.net data directory and copy index files
RUN mkdir -p /usr/share/astrometry
COPY docker/astrometry-data/ /usr/share/astrometry/

# Configure astrometry.net to use the index files
RUN echo "add_path /usr/share/astrometry" > /etc/astrometry.cfg \
    && echo "autoindex" >> /etc/astrometry.cfg

# Set up pipeline environment
ENV ORCHARD_PATH="/opt/orchard/src"
ENV PYTHONPATH="/opt/orchard/src:${PYTHONPATH:-}"
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
