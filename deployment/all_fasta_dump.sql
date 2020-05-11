select
    genomeproject.assembly_accession as unique_build_id,
    genomeproject.assembly_accession as dbkey,
    concat(REPLACE(r.definition, ', complete genome.', ''), ' [', genomeproject.assembly_accession, ']') as label,
    concat(genomeproject.gpv_directory, '/', genomeproject.filename, '_genomic.fna') as file_path
from genomeproject join replicon r on genomeproject.gpv_id = r.gpv_id and r.rep_type = 'chromosome'
where
    genomeproject.version_id = 101
    and genomeproject.file_types is not null
    and genomeproject.file_types rlike '.fna'
    and r.rep_type = 'chromosome'