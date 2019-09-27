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
