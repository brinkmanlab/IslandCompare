select
    concat(r.rep_accnum, '.', r.rep_version) as unique_build_id,
    concat(r.rep_accnum, '_', r.rep_version) as dbkey,
    concat(REPLACE(r.definition, ', complete genome.', ''), ' [', concat(r.rep_accnum, '.', r.rep_version), ']') as label,
    concat(genomeproject.gpv_directory, '/', genomeproject.filename, '_genomic.fna') as file_path
from genomeproject join replicon r on genomeproject.gpv_id = r.gpv_id and r.rep_type = 'chromosome'
where
    genomeproject.version_id = 101
    and genomeproject.file_types is not null
    and genomeproject.file_types rlike '.fna'
    and r.rep_type = 'chromosome'