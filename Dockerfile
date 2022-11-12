FROM rocker/rstudio:4.2.2
RUN apt-get update && apt-get install -y  gdal-bin git-core imagemagick libcurl4-openssl-dev libgdal-dev libgeos-dev libgeos++-dev libgit2-dev libicu-dev libpng-dev libxt-dev libproj-dev libssl-dev libxml2-dev make pandoc pandoc-citeproc zlib1g-dev libbz2-dev && rm -rf /var/lib/apt/lists/*
RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/
RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site
RUN R -e 'install.packages("remotes")'
RUN R -e 'install.packages("BiocManager")'
RUN Rscript -e 'BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", update = FALSE)'
RUN Rscript -e 'BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene",update = FALSE)'
RUN Rscript -e 'BiocManager::install("TxDb.Hsapiens.UCSC.hg18.knownGene",update = FALSE)'
RUN Rscript -e 'remotes::install_version("shinyWidgets",upgrade="never", version = "0.7.3")'
RUN Rscript -e 'remotes::install_version("shinyjs",upgrade="never", version = "2.1.0")'
RUN Rscript -e 'remotes::install_version("shinycssloaders",upgrade="never", version = "1.0.0")'
RUN Rscript -e 'BiocManager::install("Rsamtools",update = FALSE)'
RUN Rscript -e 'BiocManager::install("minfi", update = FALSE)'
RUN Rscript -e 'BiocManager::install("IlluminaHumanMethylationEPICmanifest", update = FALSE)'
RUN Rscript -e 'BiocManager::install("IlluminaHumanMethylation450kmanifest", update = FALSE)'
RUN Rscript -e 'BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19", update = FALSE)'
RUN Rscript -e 'BiocManager::install("Homo.sapiens", update = FALSE)'
RUN Rscript -e 'BiocManager::install("Gviz", update = FALSE)'
RUN Rscript -e 'BiocManager::install("epimutacionsData", update = FALSE)'
RUN Rscript -e 'BiocManager::install("epimutacions", update = FALSE)'
RUN Rscript -e 'BiocManager::install("ensembldb", update = FALSE)'
RUN Rscript -e 'remotes::install_version("cowplot",upgrade="never", version = "1.1.1")'
RUN Rscript -e 'BiocManager::install("AnnotationHub", update = FALSE)'
RUN Rscript -e 'BiocManager::install("AnnotationDbi", update = FALSE)'
RUN Rscript -e 'remotes::install_version("testthat",upgrade="never", version = "3.1.4")'
RUN Rscript -e 'remotes::install_version("golem",upgrade="never", version = "0.3.2")'
RUN Rscript -e 'remotes::install_version("shinydashboard",upgrade="never", version = "0.7.2")'
RUN Rscript -e 'BiocManager::install("ExperimentHub", update = FALSE)'
RUN Rscript -e 'remotes::install_version("reshape2",upgrade="never", version = "1.4.4")'
RUN Rscript -e 'remotes::install_version("ggrepel",upgrade="never", version = "0.9.1")'
RUN Rscript -e 'remotes::install_github("isglobal-brge/epimutacionsShiny")'

EXPOSE 8787

RUN mkdir scripts
COPY users.sh /scripts
WORKDIR /scripts
RUN chmod +x users.sh
RUN ./users.sh
RUN rm ./users.sh