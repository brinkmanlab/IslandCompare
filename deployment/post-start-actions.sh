# This must be ran AFTER tool dependencies are fully installed via admin panel or future ephemeris tool-deps install tool
# Docs: https://ephemeris.readthedocs.io/en/latest/commands.html
galaxy-wait -v
python -m ephemeris.install_tool_deps -v -g "http://localhost:80" -a admin -t $LOCAL_FILES/local_tools/tool_conf.xml
for file in $LOCAL_FILES/workflows/*.ga
do
    # Generate tool requirements
    workflow-to-tools -w "$file" -o "$file.yml" -l "IslandCompare"
    
    # Install shed tools
    shed-tools install -v -g "http://localhost:80" -a admin -t "$file.yml"
    
    # Install workflow
    workflow-install -v -g "http://localhost:80" -a admin --publish_workflows -w "$file"
done

# Patch shed tools
curl -L -s https://raw.githubusercontent.com/bgruening/galaxytools/352f2ecfb6684e174262b4681249c9faec37127b/tools/text_processing/text_processing/awk.xml -o /galaxy-central/database/shed_tools/toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/a6f147a050a2/text_processing/awk.xml

# Patch tool_deps
curl -L -s https://raw.githubusercontent.com/bioconda/bioconda-recipes/1bc6236d5260f55c8f76302b11e044f49bebf73d/recipes/mauve/MauveCM -o $GALAXY_CONDA_PREFIX/envs/__mauve@_uv_/share/mauve-2.4.0.r4736-0/MauveCM \
&& chmod +x $GALAXY_CONDA_PREFIX/envs/__mauve@_uv_/share/mauve-2.4.0.r4736-0/MauveCM \
&& ln -s $GALAXY_CONDA_PREFIX/envs/__mauve@_uv_/share/mauve-2.4.0.r4736-0/MauveCM $GALAXY_CONDA_PREFIX/envs/__mauve@_uv_/bin/MauveCM
