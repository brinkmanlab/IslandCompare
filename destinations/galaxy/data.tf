resource "galaxy_history" "data_managers" {
  name = "data_managers"
}

resource "galaxy_job" "rgi" {
  tool_id = { for tool in galaxy_repository.islandcompare["toolshed.g2.bx.psu.edu/repos/card/rgi"].tools: tool.tool_id => tool.tool_guid }["rgi_database_builder"]
  history_id = galaxy_history.data_managers.id
  params = {
    #"url" = "https://card.mcmaster.ca/download/0/broadstreet-v3.1.0.tar.bz2" TODO change to lastest after https://github.com/arpcard/rgi_wrapper/pull/13
    "url" = "https://card.mcmaster.ca/download/0/broadstreet-v3.0.7.tar.gz"
    "name" = "latest"
  }
}

resource "galaxy_job" "microbedb" {
  tool_id = { for tool in galaxy_repository.microbedb.tools: tool.tool_id => tool.tool_guid }["microbedb_all_fasta"]
  history_id = galaxy_history.data_managers.id
  params = {
    "path" = var.microbedb_path
    "builds" = true
  }
}