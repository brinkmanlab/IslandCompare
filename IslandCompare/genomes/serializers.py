from genomes.models import Genome, Gene
from rest_framework import serializers
from Bio import SeqIO
from Bio.Seq import UnknownSeq
import logging
from io import StringIO

class GenomeSerializer(serializers.ModelSerializer):
    """
    Serializer for genomes submitted by the user.
    Ensures that genome file names are unique and that gbk files are valid
    """
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

    def validate_gbk(self, value):
        """
        Ensures that gbk files contain only a single record and the record contains both annotation and sequence
        :param value:
        :return:
        """
        # 'value' is of type <class 'django.core.files.uploadedfile.TemporaryUploadedFile'>
        # SeqIO.parse cannot read this as its contents are of type bytes rather than str
        # 'gbk' is created as a decoded 'value', to be passed to SeqIO.parse
        with StringIO(value.read().decode("utf-8")) as gbk:
            value.seek(0)
            gbk_records = list(SeqIO.parse(gbk, 'genbank'))
        if len(gbk_records) != 1:
            raise serializers.ValidationError("Genbank File contains {} records".format(len(gbk_records)))
        else:
            if type(gbk_records[0].seq) is UnknownSeq:
                raise serializers.ValidationError("Unable to read sequence from Genbank File")
            if len(gbk_records[0].features) <= 1:
                raise serializers.ValidationError("Features not included in Genbank File")
        return value

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
    """
    Serializer for obtaining a gbk file from the user
    """
    genomes = serializers.FileField()

    def create(self, validated_data):
        return validated_data

    def to_internal_value(self, data):
        return data

    def to_representation(self, instance):
        return instance

    def update(self, instance, validated_data):
        return validated_data


class GenomeGenesSerializer(serializers.Serializer):
    """
    Serializer for returning genes for a genome by parsing the genbank file
    """
    logger = logging.getLogger(__name__)
    start_cut_off = None
    end_cut_off = None

    def get_genes_from_gbk(self, filePath, start_cut_off=None, end_cut_off=None):
        if start_cut_off is not None and end_cut_off is None:
            self.logger.warning("Start cut off is set but end cut off is not. Not using specified cut offs")
        if start_cut_off is None and end_cut_off is not None:
            self.logger.warning("End cut off is set but start cut off is not. Not using specified cut offs")

        # Given a path to a gbk file, this will return all genes
        geneList = []
        for record in SeqIO.parse(open(filePath), "genbank"):
            for feature in record.features:
                geneInfo = {}
                if feature.type in ["gene", "rRNA", "tRNA"]:
                    # Bio.SeqIO returns 1 for (+) and  -1 for (-)
                    geneInfo['strand'] = feature.location.strand
                    geneInfo['start'] = feature.location.start
                    geneInfo['end'] = feature.location.end
                    geneInfo['type'] = feature.type
                    try:
                        geneInfo['note'] = feature.qualifiers['note']
                    except:
                        logging.info("No Notes Found For This Gene")
                    try:
                        geneInfo['name'] = feature.qualifiers['gene'][0]
                    except:
                        logging.info("No Name Found For This Gene")
                        try:
                            geneInfo['name'] = feature.qualifiers['locus_tag'][0]
                        except:
                            logging.info("No Locus Found For This Gene")
                    if start_cut_off is not None and end_cut_off is not None:
                        if int(geneInfo['start']) > start_cut_off and int(geneInfo['end']) < end_cut_off:
                            geneList.append(geneInfo)
                    else:
                        geneList.append(geneInfo)
            # Only gather data from the first genome in a gbk file
            break
        return geneList

    def create(self, validated_data):
        return validated_data

    def to_representation(self, instance):
        self.logger.debug("Start Cut Off Used is {}\nEnd Cut Off Used is {}".format(self.start_cut_off, self.end_cut_off))
        output = self.get_genes_from_gbk(instance.gbk.path, self.start_cut_off, self.end_cut_off)
        return {'genes': output}

    def to_internal_value(self, data):
        return data

    def update(self, instance, validated_data):
        return validated_data

class GeneSerializer(serializers.Serializer):
    """
    Serializer for Gene objects
    """
    type = serializers.CharField(max_length=4)
    gene = serializers.CharField(max_length=12)
    locus_tag = serializers.CharField(max_length=12)
    product = serializers.CharField(max_length=200)
    start = serializers.IntegerField()
    end = serializers.IntegerField()
    strand = serializers.IntegerField()