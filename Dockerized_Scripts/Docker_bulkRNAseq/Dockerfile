FROM r-base:4.2.2

RUN mkdir -p /opt/software/setup/R
RUN mkdir -p /usr/local/src/brnaSeq
RUN mkdir -p /usr/local/work/

ADD ks6_docker_install.r /opt/software/setup/R
ADD ks6_docker.r /usr/local/src/brnaSeq

WORKDIR /usr/local/work/

RUN Rscript /opt/software/setup/R/ks6_docker_install.r

ENTRYPOINT ["/bin/Rscript","/usr/local/src/brnaSeq/ks6_docker.r"]

#docker run \
#    -v /path/to/data/:/usr/local/work/ \
#    -it bseq --studydesign Study_Design.txt --database EnsDb.Mmusculus.v79 \
#    --expression "infection=Group1-Group2" \
#    --ex_uh "(Group2_1+Group2_2+Group2_3+Group2_4+Group2_5)/5" \
#    --ex_eh "(Group1_1+Group1_2+Group1_3+Group1_4+Group1_5)/5" \
#    --species "Mus musculus" \
#    --out ${base}"/Output"