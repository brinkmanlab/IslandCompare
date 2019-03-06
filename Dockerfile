FROM quay.io/bgruening/galaxy:dev

USER $GALAXY_USER

RUN curl -L -s https://github.com/brinkmanlab/galaxy-tools/archive/master.tar.gz | tar xzf - --strip-components=1 -C /local_tools

ADD ./vis $GALAXY_ROOT/config/plugins/visualizations/vis
