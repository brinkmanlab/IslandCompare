FROM quay.io/bgruening/galaxy:19.01

USER $GALAXY_USER

RUN mkdir -p /home/galaxy/local_tools && curl -L -s https://github.com/brinkmanlab/galaxy-tools/archive/master.tar.gz | tar xzf - --strip-components=1 -C /home/galaxy/local_tools

COPY viz/ $GALAXY_ROOT/config/plugins/visualizations/viz

COPY colombo-v4.0-0.tar.bz2 $GALAXY_CONDA_PREFIX/colombo-v4.0-0.tar.bz2

RUN . $GALAXY_CONDA_PREFIX/bin/activate \
    && conda create --name __colombo@4.0 $GALAXY_CONDA_PREFIX/colombo-v4.0-0.tar.bz2 \
    && conda activate __colombo@4.0 \
    && conda update --all --yes --quiet \
    && conda deactivate \
    && . $GALAXY_CONDA_PREFIX/bin/deactivate

RUN rm $GALAXY_ROOT/client/galaxy/scripts/mvc/workflow/workflow-terminals.js \
    && rm $GALAXY_ROOT/client/galaxy/scripts/mvc/workflow/workflow-view-terminals.js \
    && wget -O $GALAXY_ROOT/client/galaxy/scripts/mvc/workflow/workflow-terminals.js https://raw.githubusercontent.com/mvdbeek/galaxy/6c22f13e3ad7773337f84e72ab21d295f934ee9e/client/galaxy/scripts/mvc/workflow/workflow-terminals.js \
    && wget -O $GALAXY_ROOT/client/galaxy/scripts/mvc/workflow/workflow-view-terminals.js https://raw.githubusercontent.com/mvdbeek/galaxy/6c22f13e3ad7773337f84e72ab21d295f934ee9e/client/galaxy/scripts/mvc/workflow/workflow-view-terminals.js

USER root
