from django.conf import settings
from django.db.models.signals import post_save
from django.dispatch import receiver
from rest_framework.authtoken.models import Token

# Create your models here.


@receiver(post_save, sender=settings.AUTH_USER_MODEL)
def create_auth_token(sender, instance=None, created=False, **kwargs):
    """
    Creates and stores an authentication token in the database when a user is created.
    :param sender:
    :param instance:
    :param created:
    :param kwargs:
    :return:
    """
    if created:
        Token.objects.create(user=instance)
