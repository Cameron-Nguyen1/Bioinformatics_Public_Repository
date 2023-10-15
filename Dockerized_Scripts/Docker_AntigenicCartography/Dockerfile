FROM r-base:4.2.2

RUN apt-get update && apt-get install -y \
    libssl-dev \
    libfontconfig1-dev \
    libcurl4-openssl-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libxml2-dev \
    libfreetype-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev

RUN mkdir -p /opt/software/setup/R
RUN mkdir -p /usr/local/src/agcart
RUN mkdir -p /usr/local/work/
ADD install.R /opt/software/setup/R
ADD agcart.R /usr/local/src/agcart
RUN Rscript /opt/software/setup/R/install.R

WORKDIR /usr/local/work/

ENTRYPOINT ["/bin/Rscript","/usr/local/src/agcart/agcart.R"]

#docker run \
#    -v /mnt/c/users/came/Desktop/test_data:/usr/local/work/ \
#    -it agcart --Input='TEST_DATA.csv' --XY_Lim="-10,10,-10,10" --Prefix="mystudy" --Out="MyDir" --Point_Sizes="5,2" --Transparency=".8,1" --Antigen_Overprint="FALSE"
