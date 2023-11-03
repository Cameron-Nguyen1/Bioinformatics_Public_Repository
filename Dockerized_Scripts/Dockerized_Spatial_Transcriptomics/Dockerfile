FROM r-base:4.2.2

RUN apt-get update && apt-get install -y \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libfreetype-dev \
    libhdf5-dev

RUN mkdir -p /opt/software/setup/R
RUN mkdir -p /usr/local/src/spatial
RUN mkdir -p /usr/local/work/

ADD install_deps.r /opt/software/setup/R
ADD spatial.r /usr/local/src/spatial

#REPLACE MYGITHUBPATTOKEN WITH YOUR PERSONAL ACCESS TOKEN FROM GITHUB
RUN sed -i s/"R_LIBS_USER"/"GITHUB_PAT=ghp_MYGITHUBPATTOKEN\nR_LIBS_USER"/ /etc/R/Renviron

RUN Rscript /opt/software/setup/R/install_deps.r

WORKDIR /usr/local/work/

ENTRYPOINT ["/bin/Rscript","/usr/local/src/spatial/spatial.r"]

#docker run \
# -v /mnt/d/Spatial_A1_D1:/usr/local/work/ \
# -it spatial --slides "A1,D1" --outputdir MyDir