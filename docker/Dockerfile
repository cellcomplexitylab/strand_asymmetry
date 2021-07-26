# Build and run as shown below.
# docker build -t ff .
# docker run --rm -v $(pwd):/tmp -u$(id -u):$(id -g) ff make

FROM ubuntu:focal

RUN apt-get update -y
RUN apt-get install -y wget

# Install R.
RUN DEBIAN_FRONTEND=noninteractive TZ="America/New_York" apt-get install -y \
	r-base-core=3.6.3-2 r-base-dev=3.6.3-2 \
	libfreetype6-dev libcurl4-openssl-dev libxml2-dev

# Install Python 2.7
RUN apt-get install -y python2.7 python2.7-dev && apt-get clean && rm -rf /var/lib/apt/lists/*
RUN ln /usr/bin/python2.7 /usr/bin/python

## Install CRAN R packages.
RUN R -e "install.packages('ggplot2')"
RUN R -e "install.packages('RCircos')"
# Install BioConductor R packages.
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install()"
RUN R -e "BiocManager::install('GenomicRanges')"
# Install home-made R packages.
COPY HMMt.tar.gz /
RUN R CMD INSTALL --clean HMMt.tar.gz && rm -r HMMt.tar.gz

# Install GEM mapper.
RUN wget -nv -O GEM.tar.bz2 https://sourceforge.net/projects/gemlibrary/files/gem-library/Binary%20pre-release%203/GEM-binaries-Linux-x86_64-core_i3-20130406-045632.tbz2/download && \
    tar -jxvf GEM.tar.bz2 && \
    cp GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin/gem-mapper /usr/local/bin/ && \
    cp GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin/gem-indexer /usr/local/bin/ && \
    cp GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin/gem-indexer_* /usr/local/bin/ && \
    rm -r GEM-binaries-Linux-x86_64-core_i3-20130406-045632 && \
    rm GEM.tar.bz2

# Install starcode.
RUN cd /etc/ && \
    wget https://github.com/gui11aume/starcode/archive/refs/heads/master.zip -O starcode-master.zip && \
    unzip starcode-master.zip && \
    cd starcode-master && \
    make && \
    ln -s $(pwd)/starcode /usr/local/bin/ && \
    rm ../starcode-master.zip

# Install seeq.
COPY seeq.tar.gz /
RUN tar -zxvf seeq.tar.gz -C /etc/ && \
    cd /etc/seeq && \
    make && \
    ln -s $(pwd)/seeq /usr/local/bin/ && \
    make lib && \
    python setup.py install && \
    rm /seeq.tar.gz

# Use /tmp as shared volume.
WORKDIR /tmp
