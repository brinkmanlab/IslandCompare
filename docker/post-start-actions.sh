# Docs: https://ephemeris.readthedocs.io/en/latest/commands.html
# Install shed tools
galaxy-wait -v
shed-tools install -v -g "http://localhost:80" -a admin -t /home/galaxy/workflows/IslandCompare.tools.yml
workflow-install -v -g "http://localhost:80" -a admin --publish_workflows -w /home/galaxy/workflows/IslandCompare.ga

# Patch shed tools
wget https://raw.githubusercontent.com/bgruening/galaxytools/bc9dc2fc8a806b312d679466d1ea581f3ded35ee/tools/text_processing/text_processing/awk.xml -O /galaxy-central/database/shed_tools/toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/a6f147a050a2/text_processing/awk.xml

# Patch tool_deps
wget https://raw.githubusercontent.com/bioconda/bioconda-recipes/1bc6236d5260f55c8f76302b11e044f49bebf73d/recipes/mauve/MauveCM -O $GALAXY_CONDA_PREFIX/envs/__mauve@_uv_/share/mauve-2.4.0.r4736-0/MauveCM
chmod +x $GALAXY_CONDA_PREFIX/envs/__mauve@_uv_/share/mauve-2.4.0.r4736-0/MauveCM
ln -s $GALAXY_CONDA_PREFIX/envs/__mauve@_uv_/share/mauve-2.4.0.r4736-0/MauveCM $GALAXY_CONDA_PREFIX/envs/__mauve@_uv_/bin/MauveCM
