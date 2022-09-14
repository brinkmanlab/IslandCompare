resource "galaxy_history" "data_managers" {
  name = "data_managers"
}

resource "galaxy_job" "rgi" {
  tool_id = { for tool in galaxy_repository.islandcompare["toolshed.g2.bx.psu.edu/repos/card/rgi"].tools: tool.tool_id => tool.tool_guid }["rgi_database_builder"]
  history_id = galaxy_history.data_managers.id
  params = {
    #"url" = "https://card.mcmaster.ca/download/0/broadstreet-v3.1.0.tar.bz2"
    "url" = "https://card.mcmaster.ca/latest/data"
    "name" = "latest"
  }
}

resource "galaxy_job" "microbedb" {
  tool_id = { for tool in galaxy_repository.microbedb.tools: tool.tool_id => tool.tool_guid }["microbedb_all_fasta"]
  history_id = galaxy_history.data_managers.id
  params = {
    "db" = var.microbedb_path
    "builds" = true
  }
}

resource "galaxy_job" "salmonella_gis" {
  tool_id = { for tool in galaxy_repository.data_manager_fetch_genome_dbkeys_all_fasta.tools: tool.tool_id => tool.tool_guid }["data_manager_fetch_genome_all_fasta_dbkey"]
  history_id = galaxy_history.data_managers.id
  params = {
    "dbkey_source|dbkey_source_selector" = "new"
    "dbkey_source|dbkey" = "salmonella_gis"
    "dbkey_source|dbkey_name" = "Salmonella Curated GIs [SPI1-10;SGI1-4;fels1-2;gifsy1-2]"
    "sequence_name" = "salmonella_gis"
    "sequence_id" = "salmonella_gis"
    "reference_source|reference_source_selector" = "url"
    "reference_source|user_url" = "https://github.com/brinkmanlab/IslandCompare/blob/master/workflow/workflows/data/salmonella_gis.fna"
    "sorting|sort_selector" = "as_is"
  }
}
