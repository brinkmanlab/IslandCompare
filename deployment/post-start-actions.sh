# This must be ran AFTER tool dependencies are fully installed via admin panel or future ephemeris tool-deps install tool
# Docs: https://ephemeris.readthedocs.io/en/latest/commands.html
SERVER="${SERVER:-http://localhost:80}"
KEY="${KEY:-admin}"
WORKFLOWS="${WORKFLOWS:-$LOCAL_FILES/workflows}"
TOOLCONF="${TOOLCONF:-$LOCAL_FILES/local_tools/tool_conf.xml}"

galaxy-wait -v -g "$SERVER"
install_tool_deps -v -g "$SERVER" -a $KEY -t "$TOOLCONF"
for file in $WORKFLOWS/*.ga
do
    # Generate tool requirements
    workflow-to-tools -w "$file" -o "`basename $file`.yml" -l "IslandCompare"
    
    # Install shed tools
    shed-tools install -v -g "$SERVER" -a $KEY -t "`basename $file`.yml"
    
    # Install workflow
    workflow-install -v -g "$SERVER" -a $KEY --publish_workflows -w "$file"
done
