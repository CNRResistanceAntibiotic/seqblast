###############################################
# Dockerfile to build seqblast container image
# Based on Ubuntu 16.04
# Build with:
#   sudo docker build -t seqblast .
###############################################

    # Use ubuntu 16.04 base image
    FROM ubuntu:16.04

    # File author/maintainer info
    MAINTAINER Aurélien BIRER <abirer@chu-clermontferrand.fr>

    # set non-interactive mode
    ENV DEBIAN_FRONTEND noninteractive




    #Set dependencies
    RUN apt-get update
    RUN apt-get install --no-install-recommends -y \
            file \
            unzip \
            python \
            python-pip \
            python-dev \
            python-setuptools \
            build-essential \
            ca-certificates \
            git \
            wget \
            less \
            software-properties-common \
            libdatetime-perl \
            libxml-simple-perl \
            libdigest-md5-perl \
            default-jre \
            bioperl

    RUN wget -O - http://cpanmin.us | perl - --self-upgrade
    RUN cpanm Bio::Perl

    RUN add-apt-repository ppa:webupd8team/java -y
    RUN apt-get update && \
            echo oracle-java8-installer shared/accepted-oracle-license-v1-1 select true | /usr/bin/debconf-set-selections

    RUN apt install -y oracle-java8-installer


    RUN  pip install --upgrade pip && \
         pip install biopython==1.67 && \
         pip install python-docx


    WORKDIR /usr/local/

    #Install a5
    RUN   wget https://downloads.sourceforge.net/project/ngopt/a5_miseq_linux_20160825.tar.gz && \
          gunzip a5_miseq_linux_20160825.tar.gz && \
          tar -xvf a5_miseq_linux_20160825.tar && \
          export PATH=$PATH:bin

    #Install quast
    RUN   wget https://downloads.sourceforge.net/project/quast/quast-4.4.tar.gz && \
          gunzip quast-4.4.tar.gz && \
          tar -xvf quast-4.4.tar

    #Install Trimmomatic
    RUN    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip && \
           unzip Trimmomatic-0.36.zip && \
           mv /usr/local/Trimmomatic-0.36/trimmomatic-0.36.jar /usr/local/Trimmomatic-0.36/trimmomatic.jar

    #Install SPAdes
    RUN   wget http://cab.spbu.ru/files/release3.10.0/SPAdes-3.10.0-Linux.tar.gz && \
          tar -xvf SPAdes-3.10.0-Linux.tar.gz && \
          mv /usr/local/SPAdes-3.10.0-Linux /usr/local/SPAdes-3.10.0

    #Install Blast local
    RUN   wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz && \
          gunzip ncbi-blast-2.6.0+-x64-linux.tar.gz && \
          tar -xvf ncbi-blast-2.6.0+-x64-linux.tar
    
    # the install ca-certificates and adding "slCAinfo = /etc/ssl/certs/ca-certificates.crt" to .gitconfig
    # fixed tg cloning via git with the error:
    ## fatal: unable to access 'https://github.com/vysheng/tg.git/': Problem with the SSL CA cert (path? access rights?)
    RUN echo "[http]\n\tsslVerify = true\n\tslCAinfo = /etc/ssl/certs/ca-certificates.crt\n" >> ~/.gitconfig
  
    
    RUN git clone https://github.com/tseemann/prokka.git  &&\
        cd prokka  &&\
        git checkout v1.11 &&\
        ./bin/prokka --setupdb &&\
        mv ../prokka/ ../prokka-1.11/
        


    WORKDIR /tmp/
 
    RUN    git clone https://github.com/CNRResistanceAntibiotic/seqblast.git && \
           mv seqblast/src/ /usr/local/seqblast && \
	   ls /usr/local/seqblast && \
           ln /usr/local/seqblast/seqDetector.py /usr/bin/seqDetector && \
           ln /usr/local/seqblast/seqAssembler_batch.py /usr/bin/seqAssembler && \
           ln /usr/local/seqblast/trimmer.py /usr/bin/trimmer.py && \
           ln /usr/local/seqblast/a5.py /usr/bin/a5.py && \
           ln /usr/local/seqblast/spades.py /usr/bin/spades.py && \
           ln /usr/local/seqblast/runProkka.py /usr/bin/runProkka.py && \
           ln /usr/local/seqblast/seqMLST.py /usr/bin/seqMLST.py && \
           ln /usr/local/seqblast/runBlast.py /usr/bin/runBlast.py && \
           ln /usr/local/seqblast/parseBlast.py /usr/bin/parseBlast.py && \
           ln /usr/local/seqblast/cnrdbTofasta.py /usr/bin/cnrdbTofasta.py && \
           ln /usr/local/seqblast/dbtools.py /usr/bin/dbtools.py && \
           ln /usr/local/seqblast/annotGBK.py /usr/bin/annotGBK.py && \
           ln /usr/local/seqblast/mergeResults.py /usr/bin/mergeResults.py && \

/usr/bin/mergeResults.py

           cp -r /usr/local/seqblast/db_mlst /usr/bin && \
           cp /usr/local/seqblast/adapter_trimmomatic/adapter.fasta /usr/local/Trimmomatic-0.36/adapters/
    
    RUN    apt autoremove --purge --yes && \
           apt clean
