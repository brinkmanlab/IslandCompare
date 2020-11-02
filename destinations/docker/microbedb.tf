resource "docker_container" "microbedb" {
  image = "cvmfs/service"
  name = "microbedb"
  restart = "unless-stopped"

  env = [
    "CVMFS_HTTP_PROXY=stratum-0.brinkmanlab.ca",
    "CVMFS_REPOSITORIES=microbedb.brinkmanlab.ca",
  ]

  capabilities {
    add = ["SYS_ADMIN"]
  }
  devices {
    host_path = "/dev/fuse"
  }
}