from rest_framework import serializers
from django.contrib.auth.models import User


class UserSerializer(serializers.ModelSerializer):
    """
    Serializes user account information.
    """
    class Meta:
        model = User
        fields = ('email', 'username', 'password')
        extra_kwargs = {
            'password': {'write_only': True},
        }

    def create(self, validated_data):
        """
        Creates a user account in the database from validated data
        :param validated_data:
        :return:
        """
        user = User(
            email=validated_data['email'],
            username=validated_data['username'],
        )
        user.set_password(validated_data['password'])
        user.save()
        return user

    def update(self, instance, validated_data):
        """
        Updates a user account's email address
        :param instance:
        :param validated_data:
        :return:
        """
        instance.email = validated_data['email']
        instance.save()
        return instance
