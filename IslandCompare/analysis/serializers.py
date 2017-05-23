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
    analysiscomponent_set = serializers.SerializerMethodField()

    class Meta:
        model = Analysis
        fields = ('id', 'name', 'genomes', 'submit_time', 'start_time',
                  'complete_time', 'celery_task_id', 'analysiscomponent_set')
        read_only_fields = ('id', 'genomes', 'submit_time', 'start_time',
                            'complete_time', 'celery_task_id', 'analysiscomponent_set')

    def get_analysiscomponent_set(self, obj):
        output_dict = dict()
        for analysis_component in obj.analysiscomponent_set.all():
            data = AnalysisComponentSerializer(analysis_component).data
            output_dict[analysis_component.type.name] = data
        return output_dict


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
    @staticmethod
    def get_spaced_colors(n):
        max_value = 16581375 #255**3
        interval = int(max_value / n)
        colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]

        return ['#%02x%02x%02x' % (int(i[:2], 16), int(i[2:4], 16), int(i[4:], 16)) for i in colors]

    def to_representation(self, instance):
        analysis = Analysis.objects.get(id=instance["analysis"])
        genomes = analysis.genomes

        output = dict()

        output["genomes"] = dict()
        for genome in genomes.all():
            output["genomes"][genome.id] = dict()
            output["genomes"][genome.id]["name"] = genome.name
            output["genomes"][genome.id]["length"] = instance["gbk_metadata"][str(genome.id)]['size']
            output["genomes"][genome.id]["genomic_islands"] = dict()
            output["genomes"][genome.id]["genomic_islands"]["sigi"] = [{'start': island[0], 'end': island[1]} for island in instance["sigi_gis"][str(genome.id)]]
            output["genomes"][genome.id]["genomic_islands"]["islandpath"] = [{'start': island[0], 'end':island[1]} for island in instance["islandpath_gis"][str(genome.id)]]
            output["genomes"][genome.id]["genomic_islands"]["merged"] = [{'start': island[0], 'end':island[1]} for island in instance["merge_gis"][str(genome.id)]]

        number_clusters = instance["cluster_gis"]["numberClusters"]
        color_index = self.get_spaced_colors(number_clusters)

        for genome_id in instance["gbk_paths"].keys():
            for gi_index in range(len(output["genomes"][int(genome_id)]["genomic_islands"]["merged"])):
                clusters = instance['cluster_gis'][str(genome_id)]
                cluster_index = int(clusters[str(gi_index)])
                output["genomes"][int(genome_id)]["genomic_islands"]["merged"][gi_index]['color'] = color_index[cluster_index]

        output["newick"] = instance["newick"]
        output["alignment"] = instance["alignment"]

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
        fieldnames = ['name', 'start', 'end', 'method', 'cluster_id']
        writer = csv.DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()

        method_keys = ["islandpath_gis", "sigi_gis"]

        for key in instance["gbk_paths"]:
            counter = 0
            for island in instance["merge_gis"][key]:
                genome = Genome.objects.get(id=key)
                writer.writerow({'name': genome.name,
                                 'start': island[0],
                                 'end': island[1],
                                 'method': 'merge_gis',
                                 'cluster_id': instance["cluster_gis"][str(genome.id)][str(counter)]
                                 })
                counter += 1

        for method in method_keys:
            for key in instance["gbk_paths"]:
                for island in instance[method][key]:
                    genome = Genome.objects.get(id=key)
                    writer.writerow({'name': genome.name,
                                     'start': island[0],
                                     'end': island[1],
                                     'method': method})

        contents = output.getvalue()
        output.close()

        return contents

    def to_internal_value(self, data):
        return

    def create(self, validated_data):
        return

    def update(self, instance, validated_data):
        return
