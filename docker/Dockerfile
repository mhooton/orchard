FROM --platform=linux/amd64 python:3.9

WORKDIR /usr/src/app

# Copy requirements and install Python packages
COPY docker/requirements.txt ./
RUN pip install numpy
RUN pip install --no-cache-dir -r requirements.txt

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
RUN wget https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-3.47.tar.gz \
    && tar xzf cfitsio-3.47.tar.gz \
    && cd cfitsio-3.47 \
    && ./configure --prefix=/usr/local \
    && make \
    && make install \
    && cd .. \
    && rm -rf cfitsio-3.47*

RUN wget ftp://ftp.atnf.csiro.au/pub/software/wcslib/wcslib.tar.bz2 \
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
    && make install

# Set environment variables
ENV PATH="/opt/casutools/bin:$PATH"
ENV ORCHARD_PATH="/opt/orchard/src"
ENV PYTHONPATH="/opt/orchard/src"

# Set up pipeline environment
WORKDIR /opt/orchard

# Copy pipeline code
COPY src/ ./src/

# Create pipeline data directory structure
RUN mkdir -p /data

# Set working directory
WORKDIR /opt/orchard/src

CMD ["/bin/bash"]
