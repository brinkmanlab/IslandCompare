from rest_framework import serializers
from analysis.models import Analysis, AnalysisComponent, AnalysisType
from genomes.models import Genome
from io import StringIO
import csv
from celery.result import AsyncResult
from Bio import Phylo


class AnalysisTypeSerializer(serializers.ModelSerializer):
    """
    Serializer for analysis type
    """
    class Meta:
        model = AnalysisType
        fields = ('name',)
        read_only_fields = ('name',)


class AnalysisComponentSerializer(serializers.ModelSerializer):
    """
    Serializer for analysis components
    """
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
    """
    Serializer for analysis
    """
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
    """
    Validates user genome field to ensure that genome id exists for the user.
    """
    def to_representation(self, value):
        return value

    def to_internal_value(self, data):
        genome = Genome.objects.filter(id=data, owner=self.context['request'].user)
        if not genome.exists():
            raise serializers.ValidationError("Genome with id: {}, does not exist.".format(data))
        return genome.get()


class RunAnalysisSerializer(serializers.Serializer):
    """
    Serializer for a user form to run an analysis
    """
    name = serializers.CharField(max_length=100)
    genomes = serializers.ListField(
        child=ValidGenomeField()
    )
    newick = serializers.FileField(required=False)
    gi = serializers.FileField(required=False)

    def validate_genomes(self, value):
        if len(value) <= 1:
            raise serializers.ValidationError("Need 2 or more genomes to create an analysis. Received {}"
                                              .format(len(value)))
        return value

    def validate(self, data):
        if 'newick' in data.keys():
            selected_genomes = Genome.objects.filter(id__in=[genome.id for genome in data['genomes']])

            tree = Phylo.read(StringIO(data['newick'].read().decode('utf-8')), 'newick')
            terminals = tree.get_terminals()

            for leaf in terminals:
                selected_genome = selected_genomes.filter(owner__exact=self.context['request'].user,
                                                          name__exact=leaf.name)
                if not selected_genome.exists():
                    raise serializers.ValidationError("Genome with Name: {} Does not Exist".format(leaf))
        return data

    def create(self, validated_data):
        return validated_data

    def update(self, instance, validated_data):
        instance.genomes = validated_data['genomes']
        instance.newick = validated_data['newick']
        return instance


class ReportVisualizationOverviewSerializer(serializers.Serializer):
    """
    Serializer for returning the data needed by the user to build a visualization
    """
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
        gi_types = {"sigi_gis", "islandpath_gis", "merge_gis"}.intersection(instance)
        if "cluster_gis" in instance:
            number_clusters = instance["cluster_gis"]["numberClusters"]
            color_index = self.get_spaced_colors(number_clusters)
        for genome in genomes.all():
            output["genomes"][genome.id] = dict()
            output["genomes"][genome.id]["name"] = genome.name
            output["genomes"][genome.id]["length"] = instance["gbk_metadata"][str(genome.id)]['size']
            output["genomes"][genome.id]["amr_genes"] = [{'start': amr['orf_start'], 'end': amr['orf_end'], 'strand': amr['orf_strand']} for amr in instance["amr_genes"][str(genome.id)]]
            output["genomes"][genome.id]["genomic_islands"] = dict()
            if "user_gis" in instance:
                output["genomes"][genome.id]["genomic_islands"]["user"] = instance["user_gis"][str(genome.id)]
                print(output["genomes"][genome.id]["genomic_islands"]["user"])
            else:
                for gi_type in gi_types:
                    output["genomes"][genome.id]["genomic_islands"][gi_type[:-4].replace("merge", "merged")] = [{'start': island[0], 'end': island[1]} for island in instance[gi_type][str(genome.id)]]
            if "cluster_gis" in instance:
                for gi_index in range(len(instance["merge_gis"][str(genome.id)])):
                    clusters = instance['cluster_gis'][str(genome.id)]
                    cluster_index = int(clusters[str(gi_index)])
                    output["genomes"][genome.id]["genomic_islands"]["merged"][gi_index]['color'] = color_index[cluster_index]

        output["newick"] = instance["newick"]
        output["alignment"] = instance["alignment"]
        output["failed_components"] = instance["failed_components"]

        return output

    def to_internal_value(self, data):
        return

    def create(self, validated_data):
        return

    def update(self, instance, validated_data):
        return


class ReportCsvSerializer(serializers.BaseSerializer):
    """
    Serializer needed to return a csv file to the user
    """
    def to_representation(self, instance):
        output = StringIO()
        fieldnames = ['name', 'start', 'end', 'method', 'cluster_id']
        writer = csv.DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()

        method_keys = ["islandpath_gis", "sigi_gis"]

        if "merge_gis" in instance:
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
            if method in instance:
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

class ReportGeneCsvSerializer(serializers.BaseSerializer):
    """
    Serializer to return a csv of genes contained within GIs
    """
    def to_representation(self, instance):
        output = StringIO()
        fieldnames = ['genome', 'name', 'start', 'end', 'strand', 'gi_start', 'gi_end', 'gi_method']
        writer = csv.DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()

        method_keys = ["merge_gis", "islandpath_gis", "sigi_gis"]
        methods = [method for method in method_keys if method in instance]

        for key in instance["gbk_paths"]:
            genome = Genome.objects.get(id=key)
            for method in methods:
                for island in instance[method][key]:
                    for gi_gene in genome.gene_set.filter(start__gte=island[0]).filter(end__lte=island[1]):
                        writer.writerow({'genome': genome.name,
                                         'name': gi_gene.name,
                                         'start': gi_gene.start,
                                         'end': gi_gene.end,
                                         'strand': gi_gene.strand,
                                         'gi_start': island[0],
                                         'gi_end': island[1],
                                         'gi_method': method})

        contents = output.getvalue()
        output.close()

        return contents