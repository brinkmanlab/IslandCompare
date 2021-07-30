# TODO replace with microbedb module git@github.com:brinkmanlab/MicrobeDB.git

resource "docker_image" "microbedb" {
  count = var.enable_CVMFS ? 1 : 0
  name = "cvmfs/service"
}

resource "docker_container" "microbedb" {
  count = var.enable_CVMFS ? 1 : 0
  image = docker_image.microbedb.0.latest
  name = "microbedb"
  restart = "unless-stopped"

  entrypoint = ["cvmfs2", "-d", "-f",  "-o", "allow_other", "-o", "config=/etc/cvmfs/default.d/19-brinkman.conf", "microbedb.brinkmanlab.ca", "/cvmfs/microbedb.brinkmanlab.ca"]

  #security_opts = ["apparmor=unconfined"]
  capabilities {
    add = ["SYS_ADMIN"]
  }
  devices {
    container_path = "/dev/fuse"
    host_path = "/dev/fuse"
    permissions    = "rwm"
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
    source = var.microbedb_key_path
    read_only = true
  }
}