from rest_framework import serializers
from analysis.models import Analysis, AnalysisComponent, AnalysisType
from genomes.models import Genome
from io import StringIO
import csv
import os
from celery.result import AsyncResult


class AnalysisTypeSerializer(serializers.ModelSerializer):
    class Meta:
        model = AnalysisType
        fields = ('name',)
        read_only_fields = ('name',)


class AnalysisComponentSerializer(serializers.ModelSerializer):
    status = serializers.SerializerMethodField()
    type = AnalysisTypeSerializer(read_only=True)

    class Meta:
        model = AnalysisComponent
        fields = ('id', 'start_time', 'complete_time', 'status', 'type')
        read_only_fields = ('id', 'start_time', 'complete_time', 'status', 'type')

    def get_status(self, obj):
        result = AsyncResult(obj.celery_task_id)
        return result.status


class AnalysisSerializer(serializers.ModelSerializer):
    analysiscomponent_set = AnalysisComponentSerializer(many=True, read_only=True)

    class Meta:
        model = Analysis
        fields = ('id', 'name', 'genomes', 'submit_time', 'start_time',
                  'complete_time', 'celery_task_id', 'analysiscomponent_set')
        read_only_fields = ('id', 'genomes', 'submit_time', 'start_time',
                            'complete_time', 'celery_task_id', 'analysiscomponent_set')


class ValidGenomeField(serializers.Field):
    def to_representation(self, value):
        return value

    def to_internal_value(self, data):
        genome = Genome.objects.filter(id=data, owner=self.context['request'].user)
        if not genome.exists():
            raise serializers.ValidationError("Genome with id: {}, does not exist.".format(data))
        return genome.get()


class RunAnalysisSerializer(serializers.Serializer):
    name = serializers.CharField(max_length=100)
    genomes = serializers.ListField(
        child=ValidGenomeField()
    )

    def validate_genomes(self, value):
        if len(value) <= 1:
            raise serializers.ValidationError("Need 2 or more genomes to create an analysis. Only received {}"
                                              .format(len(value)))
        return value

    def create(self, validated_data):
        return validated_data

    def update(self, instance, validated_data):
        instance.genomes = validated_data['genomes']
        return instance


class ReportVisualizationOverviewSerializer(serializers.Serializer):
    def to_representation(self, instance):
        analysis = Analysis.objects.get(id=instance["analysis"])
        genomes = analysis.genomes

        output = dict()

        output["genomes"] = dict()
        for genome in genomes.all():
            output["genomes"][genome.id] = dict()
            output["genomes"][genome.id]["name"] = genome.name

        return output

    def to_internal_value(self, data):
        return

    def create(self, validated_data):
        return

    def update(self, instance, validated_data):
        return


class ReportCsvSerializer(serializers.BaseSerializer):
    def to_representation(self, instance):
        output = StringIO()
        writer = csv.writer(output)

        writer.writerow(["IslandPath GIs"])
        writer.writerow([])

        for key in instance["islandpath_gis"].keys():
            writer.writerow([os.path.basename(instance["gbk_paths"][key])])
            for row in instance["islandpath_gis"][key]:
                writer.writerow(row)
            writer.writerow([])

        contents = output.getvalue()
        output.close()

        return contents

    def to_internal_value(self, data):
        return

    def create(self, validated_data):
        return

    def update(self, instance, validated_data):
        return
