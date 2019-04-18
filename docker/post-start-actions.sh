# Install shed tools
galaxy-wait -v
workflow-install -a admin -w /home/galaxy/workflows/IslandCompare.ga

# Patch shed tools
wget https://raw.githubusercontent.com/bgruening/galaxytools/bc9dc2fc8a806b312d679466d1ea581f3ded35ee/tools/text_processing/text_processing/awk.xml -o /shed_tools/toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/a6f147a050a2/text_processing/awk.xml
