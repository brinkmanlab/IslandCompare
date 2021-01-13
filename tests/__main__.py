import unittest
import httpimport

with httpimport.remote_repo("islandcompare", "https://raw.githubusercontent.com/brinkmanlab/islandcompare-cli/master/islandcompare.py"):
    import islandcompare

# Basic Pseudomonas

# Basic Burkholeria

# Draft without newick

# Draft with newick

# MCM Lockup
# ERR388703.gbk -r NZ_LN999987.1

# Empty BLAST failure
# VC14135 and VC13213

if __name__ == '__main__':
    unittest.main()