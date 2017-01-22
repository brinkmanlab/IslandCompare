from genomes.models import Genome
from rest_framework import serializers


class GenomeSerializer(serializers.ModelSerializer):
    class Meta:
        model = Genome
        fields = ('name', 'gbk')

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


