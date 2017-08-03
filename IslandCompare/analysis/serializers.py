from rest_framework import serializers
from analysis.models import Analysis, AnalysisComponent, AnalysisType
from genomes.models import Genome
from genomes.serializers import GenomeSerializer
from io import StringIO
import csv, ast, re
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
        cluster_flag = True
        if 'newick' in data.keys():
            selected_genomes = Genome.objects.filter(id__in=[genome.id for genome in data['genomes']])

            tree = Phylo.read(StringIO(data['newick'].read().decode('utf-8')), 'newick')
            terminals = tree.get_terminals()

            for leaf in terminals:
                leaf.name = re.sub(r'(\.genbank|\.gbff)$', ".gbk", leaf.name)
                selected_genome = selected_genomes.filter(owner__exact=self.context['request'].user,
                                                          name__exact=leaf.name)
                if not selected_genome.exists():
                    raise serializers.ValidationError("Genome: {} not included in this analysis".format(leaf))
        if 'gi' in data.keys():
            selected_genomes = [genome.name for genome in data['genomes']]
            formatted_lines = []

            # Replace carriage returns to avoid decode errors
            user_gis = re.sub(b'\r\n?', b'\n', data["gi"].read())
            # Replace problematic file extensions
            user_gis = re.sub(b'.genbank|.gbff', b'.gbk', user_gis)
            for line in user_gis.decode("utf-8").split("\n"):
                if line is not "":
                    formatted_lines.append(line)
                    columns = line.split("\t")
                    if len(columns) == 4:
                        cluster_flag = False
                    elif len(columns) != 3:
                        raise serializers.ValidationError("Improperly formatted genomic islands file")
                    if columns[0] not in selected_genomes:
                        raise serializers.ValidationError("Genome: {} not included in this analysis".format(columns[0]))

            # Replace user supplied gi data with properly formatted version
            data["gi"] = "\n".join(formatted_lines)
        # Tells AnalysisRunView whether to add MashMclClusterPipelineComponent
        data["cluster_gis"] = cluster_flag
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

        output = dict()

        output["genomes"] = dict()
        if "user_gi" in instance["pipeline_components"]:
            gi_method = "user"
        else:
            gi_method = "merge"
        # numberClusters signifies that clustering of genomic islands was run
        if "numberClusters" in instance:
            clustering = True
            color_list = self.get_spaced_colors(instance["numberClusters"])
            clusters = ast.literal_eval(analysis.clusters)
        else:
            clustering = False
        for genome in analysis.genomes.all():
            output["genomes"][genome.id] = dict()
            output["genomes"][genome.id]["name"] = genome.name
            output["genomes"][genome.id]["length"] = instance["gbk_metadata"][str(genome.id)]['size']
            output["genomes"][genome.id]["amr_genes"] = [{'start': amr['orf_start'], 'end': amr['orf_end'], 'strand': amr['orf_strand']} for amr in instance["amr_genes"][str(genome.id)]]
            output["genomes"][genome.id]["genomic_islands"] = dict()
            if gi_method == "user":
                gis = genome.genomicisland_set.filter(method="user").filter(usergenomicisland__analysis=analysis)
                gi_dicts = [{'start': gi.start, 'end': gi.end, 'color': gi.usergenomicisland.color} for gi in gis]
                output["genomes"][genome.id]["genomic_islands"]["user"] = gi_dicts
            else:
                for gi_type in ["sigi", "islandpath", "merge"]:
                    gis = genome.genomicisland_set.filter(method=gi_type)
                    if gis.exists():
                        gi_dicts = [{'start': gi.start, 'end': gi.end} for gi in gis]
                        output["genomes"][genome.id]["genomic_islands"][gi_type] = gi_dicts
            if clustering and gi_method in output["genomes"][genome.id]["genomic_islands"]:
                for gi_index in range(len(output["genomes"][genome.id]["genomic_islands"][gi_method])):
                    cluster_index = int(clusters[str(genome.id)][str(gi_index)])
                    output["genomes"][genome.id]["genomic_islands"][gi_method][gi_index]["cluster"] = cluster_index
                    output["genomes"][genome.id]["genomic_islands"][gi_method][gi_index]["color"] = color_list[cluster_index]

        output["newick"] = instance["newick"]
        output["alignment"] = instance["alignment"]
        output["failed_components"] = instance["failed_components"]
        output["cluster_gis"] = clustering and gi_method # 'user', 'merge', or False

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
    def to_representation(self, analysis):
        output = StringIO()
        fieldnames = ['name', 'start', 'end', 'method', 'cluster_id']
        writer = csv.DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()

        method_keys = ["merge", "islandpath", "sigi"]

        for method in method_keys:
            for genome in analysis.genomes.all():
                for island in genome.genomicisland_set.filter(method=method):
                    row = {
                        'name': genome.name,
                        'start': island.start,
                        'end': island.end,
                        'method': method + "_gis"
                    }
                    if method == "merge":
                        cluster = island.genomicislandcluster_set.get(analysis=analysis)
                        row['cluster_id'] = cluster.number
                    writer.writerow(row)

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
    def to_representation(self, analysis):
        output = StringIO()
        fieldnames = ['genome', 'gene_name', 'locus_tag', 'product', 'start', 'end', 'strand', 'gi_start', 'gi_end', 'gi_method']
        writer = csv.DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()

        method_keys = ["merge", "islandpath", "sigi"]

        for genome in analysis.genomes.all():
            for method in method_keys:
                for island in genome.genomicisland_set.filter(method=method):
                    for gi_gene in genome.gene_set.filter(start__gte=island.start).filter(end__lte=island.end):
                        writer.writerow({'genome': genome.name,
                                         'gene_name': gi_gene.gene,
                                         'locus_tag': gi_gene.locus_tag,
                                         'product': gi_gene.product,
                                         'start': gi_gene.start,
                                         'end': gi_gene.end,
                                         'strand': gi_gene.strand,
                                         'gi_start': island.start,
                                         'gi_end': island.end,
                                         'gi_method': method + "_gis"})

        contents = output.getvalue()
        output.close()

        return contents

class AnalysisGenomicIslandSerializer(serializers.Serializer):
    """
    Serializer for GenomicIsland objects
    """
    method = serializers.CharField(max_length=10)
    start = serializers.IntegerField()
    end = serializers.IntegerField()
    genome = GenomeSerializer()
