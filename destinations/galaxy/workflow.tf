locals {
  workflow = file("${path.module}/../../workflow/workflows/IslandCompare.ga")
}

resource "galaxy_stored_workflow" "islandcompare" {
  json = local.workflow
  published = true
}