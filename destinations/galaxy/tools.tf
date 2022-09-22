data "galaxy_workflow_repositories" "islandcompare" {
  json = local.workflow
}

resource "galaxy_repository" "islandcompare" {
  for_each = { for repo in data.galaxy_workflow_repositories.islandcompare.repositories: "${repo.tool_shed}/repos/${repo.owner}/${repo.name}" => repo }
  tool_shed = each.value.tool_shed
  owner = each.value.owner
  name = each.value.name
  changeset_revision = each.value.changeset_revision

  lifecycle {
    ignore_changes = [changeset_revision,install_repository_dependencies,install_resolver_dependencies,install_tool_dependencies,remove_from_disk,sub_repositories,tools]
  }
}

resource "galaxy_repository" "microbedb" {
  tool_shed = "toolshed.g2.bx.psu.edu"
  owner = "brinkmanlab"
  name = "microbedb"
  changeset_revision = "2f6ef3a184df"
  lifecycle {
    ignore_changes = [changeset_revision,install_repository_dependencies,install_resolver_dependencies,install_tool_dependencies,remove_from_disk,sub_repositories,tools]
  }
}

resource "galaxy_repository" "data_manager_fetch_genome_dbkeys_all_fasta" {
  tool_shed = "toolshed.g2.bx.psu.edu"
  owner = "devteam"
  name = "data_manager_fetch_genome_dbkeys_all_fasta"
  changeset_revision = "4d3eff1bc421"
  lifecycle {
    ignore_changes = [changeset_revision,install_repository_dependencies,install_resolver_dependencies,install_tool_dependencies,remove_from_disk,sub_repositories,tools]
  }
}
