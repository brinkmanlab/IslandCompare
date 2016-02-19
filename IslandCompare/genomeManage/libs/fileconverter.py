from Bio import SeqIO

# Code modified version of snippet from
# http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank2fasta/
def convertGbkToFna(gbkpath,faapath):
    input_handle  = open(gbkpath, "r")
    output_handle = open(faapath, "w")

    for seq_record in SeqIO.parse(input_handle, "genbank") :
        output_handle.write(">%s %s\n%s\n" % (
            seq_record.id,
            seq_record.description,
            seq_record.seq))

    output_handle.close()
    input_handle.close()
    return faapath
