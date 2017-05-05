from django.test import TestCase
from rest_framework.test import APIRequestFactory, force_authenticate
from genomes.models import Genome
from genomes.views import GenomeListView
from django.contrib.auth.models import User
from django.core.files.uploadedfile import SimpleUploadedFile
import os
from rest_framework.reverse import reverse
from rest_framework.test import APIClient

# Create your tests here.


class ListGenomesTestCase(TestCase):
    test_username = "username"
    test_user = None

    test_name = "test_genome"
    test_gbk = None
    test_genome = None

    def setUp(self):
        self.factory = APIRequestFactory()

        self.test_user = User(username=self.test_username)
        self.test_user.save()

        self.test_genome = Genome(name=self.test_name,
                                  owner=self.test_user,
                                  gbk=self.test_gbk)
        self.test_genome.save()

    def test_authenticated_list_genomes(self):
        url = reverse('genome')

        request = self.factory.get(url)
        force_authenticate(request, user=self.test_user)
        response = GenomeListView.as_view()(request)

        self.assertEqual(200, response.status_code)
        self.assertEqual(self.test_name, response.data[0]['name'])
        self.assertEqual(self.test_gbk, response.data[0]['gbk'])

    def test_unauthenticated_list_genomes(self):
        url = reverse('genome')

        request = self.factory.get(url)
        response = GenomeListView.as_view()(request)

        self.assertEqual(403, response.status_code)


class CreateGenomeTestCase(TestCase):
    test_username = "username"
    test_user = None

    new_name = "test_genome"
    new_gbk_name = "test.gbk"
    new_gbk_contents = bytes("test", 'utf-8')
    new_gbk = SimpleUploadedFile(new_gbk_name, new_gbk_contents)
    new_genome = None

    def setUp(self):
        self.factory = APIRequestFactory()

        self.test_user = User(username=self.test_username)
        self.test_user.save()

    def test_authenticated_create_genome(self):
        url = reverse('genome')

        request = self.factory.post(url,
                                    {'name': self.new_name,
                                     'gbk': self.new_gbk})
        force_authenticate(request, user=self.test_user)
        response = GenomeListView.as_view()(request)

        self.assertEqual(201, response.status_code)
        self.assertTrue(Genome.objects.filter(name=self.new_name).exists())

        genome = Genome.objects.filter(name=self.new_name).get()
        self.assertEqual(self.new_gbk_contents, genome.gbk.read())

    def test_unauthenticated_create_genome(self):
        url = reverse('genome')

        request = self.factory.post(url,
                                    {'name': self.new_name,
                                     'gbk': self.new_gbk})
        response = GenomeListView.as_view()(request)

        self.assertEqual(403, response.status_code)

    def test_gbk_delete_on_model_delete(self):
        genome = Genome(name=self.new_name,
                        owner=self.test_user,
                        gbk=self.new_gbk)
        genome.save()
        gbk_path = genome.gbk.name
        genome.delete()

        self.assertFalse(os.path.isfile(gbk_path))

    def test_gbk_none_delete_on_model_delete(self):
        genome = Genome(name=self.new_name,
                        owner=self.test_user,
                        gbk=None)
        genome.save()
        genome.delete()

    def tearDown(self):
        for genome in Genome.objects.all():
            genome.delete()


class RetrieveGenomeTestCase(TestCase):
    test_username = "username"
    test_user = None

    test_name = "test_genome"
    test_gbk = None
    test_genome = None

    def setUp(self):
        self.factory = APIRequestFactory()

        self.test_user = User(username=self.test_username)
        self.test_user.save()

        self.test_genome = Genome(name=self.test_name,
                                  owner=self.test_user,
                                  gbk=self.test_gbk)
        self.test_genome.save()

    def test_authenticated_retrieve_genome(self):
        url = reverse('genome_details', kwargs={'pk': self.test_genome.id})

        client = APIClient()
        client.force_authenticate(user=self.test_user)
        response = client.get(url)

        self.assertEqual(200, response.status_code)
        self.assertEqual(self.test_name, response.data['name'])
        self.assertEqual(self.test_gbk, response.data['gbk'])

    def test_unauthenticated_retrieve_genome(self):
        url = reverse('genome_details', kwargs={'pk': self.test_genome.id})

        client = APIClient()
        response = client.get(url)

        self.assertEqual(403, response.status_code)

    def test_incorrect_user_retrieve_genome(self):
        url = reverse('genome_details', kwargs={'pk': self.test_genome.id})

        incorrect_user = User(username="test_user")
        incorrect_user.save()

        client = APIClient()
        client.force_authenticate(user=incorrect_user)
        response = client.get(url)

        self.assertEqual(404, response.status_code)

    def test_retrieve_genome_does_not_exist(self):
        url = reverse('genome_details', kwargs={'pk': self.test_genome.id+1})

        client = APIClient()
        client.force_authenticate(user=self.test_user)
        response = client.get(url)

        self.assertEqual(404, response.status_code)


class UpdateGenomeTestCase(TestCase):
    test_username = "username"
    test_user = None

    test_name = "test_genome"
    test_gbk = None
    test_genome = None

    def setUp(self):
        self.factory = APIRequestFactory()

        self.test_user = User(username=self.test_username)
        self.test_user.save()

        self.test_genome = Genome(name=self.test_name,
                                  owner=self.test_user,
                                  gbk=self.test_gbk)
        self.test_genome.save()

    def test_authenticated_update_genome(self):
        url = reverse('genome_details', kwargs={'pk': self.test_genome.id})

        updated_name = "updated_genome"
        updated_gbk_name = "test.gbk"
        updated_gbk_content = bytes("test", 'utf-8')
        updated_gbk = SimpleUploadedFile(updated_gbk_name, updated_gbk_content)

        client = APIClient()
        client.force_authenticate(user=self.test_user)
        response = client.put(url,
                              {'name': updated_name,
                               'gbk': updated_gbk})

        genome = Genome.objects.get(id=self.test_genome.id)

        self.assertEqual(200, response.status_code)
        self.assertEqual(updated_name, genome.name)
        self.assertEqual(updated_gbk_content, genome.gbk.read())

    def test_unauthenticated_update_genome(self):
        url = reverse('genome_details', kwargs={'pk': self.test_genome.id})

        updated_name = "updated_genome"
        updated_gbk_name = "test.gbk"
        updated_gbk_content = bytes("test", 'utf-8')
        updated_gbk = SimpleUploadedFile(updated_gbk_name, updated_gbk_content)

        client = APIClient()
        response = client.put(url,
                              {'name': updated_name,
                               'gbk': updated_gbk})

        self.assertEqual(403, response.status_code)

    def test_incorrect_user_update_genome(self):
        url = reverse('genome_details', kwargs={'pk': self.test_genome.id})

        incorrect_user = User(username="incorrect_user")
        incorrect_user.save()

        updated_name = "updated_genome"
        updated_gbk_name = "test.gbk"
        updated_gbk_content = bytes("test", 'utf-8')
        updated_gbk = SimpleUploadedFile(updated_gbk_name, updated_gbk_content)

        client = APIClient()
        client.force_authenticate(user=incorrect_user)
        response = client.put(url,
                              {'name': updated_name,
                               'gbk': updated_gbk})

        self.assertEqual(404, response.status_code)

    def test_update_genome_does_not_exist(self):
        url = reverse('genome_details', kwargs={'pk': self.test_genome.id+1})

        updated_name = "updated_genome"
        updated_gbk_name = "test.gbk"
        updated_gbk_content = bytes("test", 'utf-8')
        updated_gbk = SimpleUploadedFile(updated_gbk_name, updated_gbk_content)

        client = APIClient()
        client.force_authenticate(user=self.test_user)
        response = client.put(url,
                              {'name': updated_name,
                               'gbk': updated_gbk})

        self.assertEqual(404, response.status_code)

    def tearDown(self):
        for genome in Genome.objects.all():
            genome.delete()


class DeleteGenomeTestCase(TestCase):
    test_username = "username"
    test_user = None

    test_name = "test_genome"
    test_gbk_name = "test.gbk"
    test_gbk_contents = bytes("test", 'utf-8')
    test_gbk = SimpleUploadedFile(test_gbk_name, test_gbk_contents)
    test_genome = None

    def setUp(self):
        self.factory = APIRequestFactory()

        self.test_user = User(username=self.test_username)
        self.test_user.save()

        self.test_genome = Genome(name=self.test_name,
                                  owner=self.test_user,
                                  gbk=self.test_gbk)
        self.test_genome.save()

    def test_authenticated_delete_genome(self):
        url = reverse('genome_details', kwargs={'pk': self.test_genome.id})

        client = APIClient()
        client.force_authenticate(user=self.test_user)
        response = client.delete(url)

        self.assertEqual(204, response.status_code)
        self.assertFalse(Genome.objects.filter(id=self.test_genome.id).exists())

    def test_unauthenticated_delete_genome(self):
        url = reverse('genome_details', kwargs={'pk': self.test_genome.id})

        client = APIClient()
        response = client.delete(url)

        self.assertEqual(403, response.status_code)

    def test_incorrect_user_delete_genome(self):
        url = reverse('genome_details', kwargs={'pk': self.test_genome.id})

        incorrect_user = User(username="incorrect_user")
        incorrect_user.save()

        client = APIClient()
        client.force_authenticate(user=incorrect_user)
        response = client.delete(url)

        self.assertEqual(404, response.status_code)

    def test_delete_genome_does_not_exist(self):
        url = reverse('genome_details', kwargs={'pk': self.test_genome.id+1})

        client = APIClient()
        client.force_authenticate(user=self.test_user)
        response = client.delete(url)

        self.assertEqual(404, response.status_code)

    def tearDown(self):
        for genome in Genome.objects.all():
            genome.delete()
