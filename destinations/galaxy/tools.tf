data "galaxy_workflow_repositories" "islandcompare" {
  json = local.workflow
}

resource "galaxy_repository" "islandcompare" {
  for_each = { for repo in data.galaxy_workflow_repositories.islandcompare.repositories: "${repo.tool_shed}/repos/${repo.owner}/${repo.name}" => repo }
  tool_shed = each.value.tool_shed
  owner = each.value.owner
  name = each.value.name
  changeset_revision = each.value.changeset_revision
}

resource "galaxy_repository" "microbedb" {
  tool_shed = "toolshed.g2.bx.psu.edu"
  owner = "brinkmanlab"
  name = "microbedb"
  changeset_revision = "7c2f036e13a7"
}