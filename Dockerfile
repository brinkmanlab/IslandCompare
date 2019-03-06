FROM quay.io/bgruening/galaxy:dev

USER $GALAXY_USER

RUN mkdir -p /home/galaxy/local_tools && curl -L -s https://github.com/brinkmanlab/galaxy-tools/archive/master.tar.gz | tar xzf - --strip-components=1 -C /home/galaxy/local_tools

COPY ./vis/ $GALAXY_ROOT/config/plugins/visualizations/vis

COPY ./colombo-v4.0-0.tar.bz2 $GALAXY_CONDA_PREFIX/colombo-v4.0-0.tar.bz2

RUN source $GALAXY_CONDA_PREFIX/bin/activate \
    && conda create --name __colombo@4.0 $GALAXY_CONDA_PREFIX/colombo-v4.0-0.tar.bz2 \
    && conda activate __colombo@4.0 \
    && conda update --all --yes --quiet \
    && conda deactivate \
    && source $GALAXY_CONDA_PREFIX/bin/deactivate
    
