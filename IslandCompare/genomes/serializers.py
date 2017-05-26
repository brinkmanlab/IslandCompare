from genomes.models import Genome
from rest_framework import serializers
from Bio import SeqIO
import logging


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

    def validate_gbk(self, value):
        gbk_records = SeqIO.parse(value, 'genbank')
        if len(list(gbk_records)) > 1:
            raise serializers.ValidationError("Genbank File contains {} records".format(len(list(gbk_records))))
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
    logger = logging.getLogger(__name__)
    start_cut_off = None
    end_cut_off = None

    def get_genes_from_gbk(self, filePath, start_cut_off=None, end_cut_off=None):
        if start_cut_off is not None and end_cut_off is None:
            self.logger.warning("Start cut off is set but end cut off is not. Not using specified cut offs")
        if start_cut_off is None and end_cut_off is not None:
            self.logger.warning("End cut off is set but start cut off is not. Not using specified cut offs")

        # Given a path to a gbk file, this will return all CDS
        geneList = []
        for record in SeqIO.parse(open(filePath), "genbank"):
            for feature in record.features:
                geneInfo = {}
                if feature.type == 'gene' or feature.type == 'CDS':
                    # Bio.SeqIO returns 1 for (+) and  -1 for (-)
                    geneInfo['strand'] = feature.location.strand
                    geneInfo['start'] = feature.location.start
                    geneInfo['end'] = feature.location.end
                    try:
                        geneInfo['note'] = feature.qualifiers['note']
                    except:
                        logging.info("No Notes Found For This Gene")
                    try:
                        geneInfo['name'] = feature.qualifiers['gene'][0]
                    except:
                        logging.info("No Name Found For This Gene")
                        try:
                            geneInfo['name'] = feature.qualifiers['locus_tag']
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