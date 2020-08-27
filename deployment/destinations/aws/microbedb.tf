resource "aws_iam_user" "S3Reader" {
  name = "S3Reader"
  path = "/${local.instance}/"
}

data "aws_iam_policy_document" "S3Reader" {
  statement {
    sid    = "S3Reader"
    effect = "Allow"
    actions = [
      "s3:GetObject",
      "s3:ListBucket"
    ]

    resources = ["*"]
  }
}

resource "aws_iam_policy" "S3Reader" {
  name        = "S3Reader"
  path        = "/${local.instance}/"
  description = "Autoscaler cluster policy"

  policy = data.aws_iam_policy_document.S3Reader.json
}

resource "aws_iam_user_policy_attachment" "S3Reader" {
  user       = aws_iam_user.S3Reader.name
  policy_arn = aws_iam_policy.S3Reader.arn
}

resource "aws_iam_access_key" "S3Reader" {
  user = aws_iam_user.S3Reader.name
}

data "aws_region" "current" {}

resource "kubernetes_job" "microbedb" {
  metadata {
    generate_name = "init-microbedb-"
    namespace     = local.instance
  }
  spec {
    template {
      metadata {}
      spec {
        security_context {
          run_as_user = 1000
          run_as_group = 1000
          fs_group = 1000
        }
        container {
          name  = "init-microbedb"
          image = "rclone/rclone"
          command = ["sh", "-c"]
          args = [<<-EOF
            rclone config create aws s3 provider AWS env_auth false access_key_id '${aws_iam_access_key.S3Reader.id}' secret_access_key '${aws_iam_access_key.S3Reader.secret}' region '${data.aws_region.current.name}' >> /dev/null &&
            rclone copy aws:/microbedb ${var.data_dir}/microbedb -v &&
            sed -i 's/\/data\/NCBI_genomes\/MicrobeDBv2\//\/${var.data_dir}\/microbedb\//' /${var.data_dir}/microbedb/all_fasta.loc
            EOF
          ]
          volume_mount {
            mount_path = var.data_dir
            name       = "data"
          }
        }
        node_selector = {
          WorkClass = "service"
        }
        volume {
          name = "data"
          persistent_volume_claim {
            claim_name = var.user_data_volume_name
          }
        }
        restart_policy = "Never"
      }
    }
    backoff_limit = 1
  }
  wait_for_completion = true
  timeouts {
    create = "10m"
  }
}