from django.test import TestCase
from rest_framework.test import APIRequestFactory, force_authenticate
from accounts.views import UserRegistrationView, UserRetrieveUpdateView
from django.contrib.auth.models import User

# Create your tests here.


class RegistrationTestCase(TestCase):
    new_email = "email@test.ca"
    new_username = "username"
    new_password = "password"
    existing_email = "exist@test.ca"
    existing_username = "existing"
    existing_password = "existing_password"

    def setUp(self):
        self.factory = APIRequestFactory()

        user = User(username=self.existing_username,
                    email=self.existing_email)
        user.set_password(self.existing_password)
        user.save()

    def test_success_registration(self):
        request = self.factory.post('/accounts/registration',
                                    {'email': self.new_email,
                                     'username': self.new_username,
                                     'password': self.new_password})
        response = UserRegistrationView.as_view()(request)

        self.assertTrue(User.objects.filter(username=self.new_username,
                                            email=self.new_email).exists())

        self.assertEqual(201, response.status_code)

    def test_duplicate_registration(self):
        request = self.factory.post('/accounts/registration',
                                    {'email': self.existing_email,
                                     'username': self.existing_username,
                                     'password': self.existing_password})
        response = UserRegistrationView.as_view()(request)

        self.assertEqual(400, response.status_code)


class UserRetrieveTestCase(TestCase):
    test_email = "email@test.ca"
    test_username = "username"
    test_password = "password"
    test_user = None

    def setUp(self):
        self.factory = APIRequestFactory()

        self.test_user = User(username=self.test_username,
                              email=self.test_email)
        self.test_user.set_password(self.test_password)
        self.test_user.save()

    def test_authenticated_retrieve(self):
        request = self.factory.get('/accounts/')
        force_authenticate(request, user=self.test_user)
        response = UserRetrieveUpdateView.as_view()(request)

        self.assertEqual(200, response.status_code)
        self.assertEqual(self.test_email, response.data['email'])
        self.assertEqual(self.test_username, response.data['username'])

    def test_unauthenticated_retrieve(self):
        request = self.factory.get('/accounts/')
        response = UserRetrieveUpdateView.as_view()(request)

        self.assertEqual(403, response.status_code)


class UserUpdateTestCase(TestCase):
    test_email = "email@test.ca"
    test_username = "username"
    test_password = "password"
    test_user = None

    def setUp(self):
        self.factory = APIRequestFactory()

        self.test_user = User(username=self.test_username,
                              email=self.test_email)
        self.test_user.set_password(self.test_password)
        self.test_user.save()

    def test_authenticated_update(self):
        updated_email = "update@test.ca"

        request = self.factory.put('/accounts/',
                                   {'username': self.test_username,
                                    'email': updated_email,
                                    'password': self.test_password})
        force_authenticate(request, user=self.test_user)
        response = UserRetrieveUpdateView.as_view()(request)

        self.assertEqual(200, response.status_code)
        self.assertEqual(updated_email, self.test_user.email)

    def test_unauthenticated_update(self):
        unauthenticated_email = "unauthenticated@test.ca"

        request = self.factory.put('/accounts/',
                                   {'username': self.test_username,
                                    'email': unauthenticated_email,
                                    'password': self.test_password})
        response = UserRetrieveUpdateView.as_view()(request)

        self.assertEqual(403, response.status_code)
        self.assertEqual(self.test_email, self.test_user.email)
