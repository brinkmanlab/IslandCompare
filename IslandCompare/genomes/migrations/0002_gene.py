# -*- coding: utf-8 -*-
# Generated by Django 1.10.5 on 2017-07-13 21:23
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('genomes', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Gene',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=50)),
                ('start', models.IntegerField()),
                ('end', models.IntegerField()),
                ('strand', models.SmallIntegerField()),
                ('type', models.CharField(max_length=4)),
                ('genome', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='genomes.Genome')),
            ],
        ),
    ]