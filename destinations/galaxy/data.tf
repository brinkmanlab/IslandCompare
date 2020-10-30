resource "galaxy_history" "data_managers" {
  name = "data_managers"
}

resource "galaxy_job" "microbedb" {
  depends_on = [galaxy_repository.microbedb]
  tool_id = "toolshed.g2.bx.psu.edu/repos/brinkmanlab/microbedb/microbedb_all_fasta/1.0"
  #tool_id = galaxy_repository.microbedb.tools["load_fasta"].id
  history_id = galaxy_history.data_managers.id
  params = {
    "source" = "path"
    "source|path" = "${var.data_dir}/microbedb/microbe.sqlite"
  }
}

resource "galaxy_job" "rgi" {
  tool_id = "toolshed.g2.bx.psu.edu/repos/card/rgi/rgi_database_builder/1.0.0"
  #tool_id = [get rgi from galaxy_repository.islandcompare[]].tools["rgi_database_builder"].id
  history_id = galaxy_history.data_managers.id
  params = {
    "url" = "https://card.mcmaster.ca/download/1/software-v5.1.1.tar.bz2"
    "name" = "5.1.1"
  }
}