# Docs: https://ephemeris.readthedocs.io/en/latest/commands.html
# Install shed tools
galaxy-wait -v
shed-tools install -v -g "http://localhost:80" -a admin -t /home/galaxy/workflows/IslandCompare.tools.yml
workflow-install -v -g "http://localhost:80" -a admin --publish_workflows -w /home/galaxy/workflows/IslandCompare.ga

# Patch shed tools
wget https://raw.githubusercontent.com/bgruening/galaxytools/bc9dc2fc8a806b312d679466d1ea581f3ded35ee/tools/text_processing/text_processing/awk.xml -O /galaxy-central/database/shed_tools/toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/a6f147a050a2/text_processing/awk.xml
