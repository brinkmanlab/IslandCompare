FROM quay.io/bgruening/galaxy:19.01

USER $GALAXY_USER

RUN mkdir -p /home/galaxy/local_tools && curl -L -s https://github.com/brinkmanlab/galaxy-tools/archive/master.tar.gz | tar xzf - --strip-components=1 -C /home/galaxy/local_tools

COPY viz/ $GALAXY_ROOT/config/plugins/visualizations/islandcompare

COPY colombo-v4.0-0.tar.bz2 $GALAXY_CONDA_PREFIX/colombo-v4.0-0.tar.bz2

RUN . $GALAXY_CONDA_PREFIX/bin/activate \
    && conda create --name __colombo@4.0 $GALAXY_CONDA_PREFIX/colombo-v4.0-0.tar.bz2 \
    && conda activate __colombo@4.0 \
    && conda update --all --yes --quiet \
    && conda deactivate \
    && . $GALAXY_CONDA_PREFIX/bin/deactivate

RUN curl -L -s https://patch-diff.githubusercontent.com/raw/galaxyproject/galaxy/pull/7543.diff | patch -f --verbose -p1 -d $GALAXY_ROOT \
    && curl -L -s https://patch-diff.githubusercontent.com/raw/galaxyproject/galaxy/pull/7435.diff | patch -f --verbose -p1 -d $GALAXY_ROOT \
    && curl -L -s https://patch-diff.githubusercontent.com/raw/galaxyproject/galaxy/pull/7458.diff | patch -f --verbose -p1 -d $GALAXY_ROOT \
    && curl -L -s https://patch-diff.githubusercontent.com/raw/galaxyproject/galaxy/pull/7453.diff | patch -f --verbose -p1 -d $GALAXY_ROOT \
    && make -C $GALAXY_ROOT client

USER root
