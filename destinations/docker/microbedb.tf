resource "docker_container" "microbedb" {
  image = "cvmfs/service"
  name = "microbedb"
  restart = "unless-stopped"

  entrypoint = ["cvmfs2", "-d", "-f",  "-o", "allow_other", "-o", "config=/etc/cvmfs/default.d/19-brinkman.conf", "microbedb.brinkmanlab.ca", "/cvmfs/microbedb.brinkmanlab.ca"]

  security_opts = ["apparmor=unconfined"]
  capabilities {
    add = ["SYS_ADMIN"]
  }
  devices {
    host_path = "/dev/fuse"
  }
  mounts {
    type = "bind"
    target = "/cvmfs/microbedb.brinkmanlab.ca"
    source = var.microbedb_mount_path
    bind_options {
      propagation = "shared"
    }
    read_only = true
  }
  mounts {
    type = "bind"
    target = "/etc/cvmfs/default.d/19-brinkman.conf"
    source = "${abspath(path.module)}/cvmfs.config"
    read_only = true
  }
  mounts {
    type = "bind"
    target = "/etc/cvmfs/keys/microbedb.brinkmanlab.ca.pub"
    source = "${abspath(path.module)}/microbedb.brinkmanlab.ca.pub"
    read_only = true
  }
}