from genomes.models import Genome
from rest_framework import serializers


class GenomeSerializer(serializers.ModelSerializer):
    class Meta:
        model = Genome
        fields = ('id', 'name', 'gbk')
        read_only_fields = ('id',)

    def validate(self, data):
        if Genome.objects.filter(
            name=data['name'],
            owner=self.context['request'].user
        ).exists():
            raise serializers.ValidationError("Genome with name: {} already exists".format(data['name']))
        return data

    def create(self, validated_data):
        genome = Genome(
            name=validated_data['name'],
            gbk=validated_data['gbk'],
            owner=self.context['request'].user
        )
        genome.save()
        return genome

    def update(self, instance, validated_data):
        instance.name = validated_data['name']
        instance.gbk = validated_data['gbk']
        instance.save()
        return instance


class GenomeUploadSerializer(serializers.Serializer):
    genomes = serializers.FileField()

    def create(self, validated_data):
        return validated_data

    def to_internal_value(self, data):
        return data

    def to_representation(self, instance):
        return instance

    def update(self, instance, validated_data):
        return validated_data
