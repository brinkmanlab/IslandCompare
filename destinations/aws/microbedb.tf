data "aws_iam_policy_document" "s3reader_assume_role_with_oidc" {
  statement {
    effect = "Allow"

    actions = ["sts:AssumeRoleWithWebIdentity"]

    principals {
      type        = "Federated"
      identifiers = [var.eks.oidc_provider_arn]
    }

    condition {
      test     = "StringEquals"
      variable = "${trimprefix(var.eks.cluster_oidc_issuer_url, "https://")}:sub"
      values   = ["system:serviceaccount:${local.namespace.metadata.0.name}:s3reader"]
    }

    condition {
      test     = "StringEquals"
      variable = "${trimprefix(var.eks.cluster_oidc_issuer_url, "https://")}:aud"
      values   = ["sts.amazonaws.com"]
    }
  }
}

resource "aws_iam_role" "S3Reader" {
  name_prefix = "S3Reader"
  assume_role_policy = data.aws_iam_policy_document.s3reader_assume_role_with_oidc.json
}

resource "aws_iam_role_policy_attachment" "S3Reader" {
  role       = aws_iam_role.S3Reader.name
  policy_arn = "arn:aws:iam::aws:policy/AmazonS3ReadOnlyAccess"
}

resource "kubernetes_service_account" "S3Reader" {
  metadata {
    name      = "s3reader"
    namespace = local.namespace.metadata.0.name
    annotations = {
      "eks.amazonaws.com/role-arn" = aws_iam_role.S3Reader.arn
    }
  }
}

resource "kubernetes_job" "microbedb" {
  metadata {
    generate_name = "init-microbedb-"
    namespace     = local.instance
  }
  spec {
    template {
      metadata {}
      spec {
        service_account_name = kubernetes_service_account.S3Reader.metadata.0.name
        automount_service_account_token = true
        security_context {
          run_as_user = var.uwsgi_uid
          run_as_group = var.uwsgi_gid
          fs_group = var.uwsgi_gid
        }
        container {
          name  = "init-microbedb"
          image = "rclone/rclone"
          command = ["sh", "-c"]
          args = [<<-EOF
            rclone config create aws s3 provider AWS env_auth true region '${data.aws_region.current.name}' --config /tmp/rclone.conf >> /dev/null &&
            rclone copy aws:/microbedb ${var.data_dir}/microbedb -v --config /tmp/rclone.conf
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
    create = "20m"
    update = "20m"
  }
}

data "null_data_source" "microbedb" {
  depends_on = [kubernetes_job.microbedb]
  inputs = {
    path = "${var.data_dir}/microbedb/microbedb.sqlite"
  }
}