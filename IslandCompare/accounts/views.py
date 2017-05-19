from rest_framework import generics
from rest_framework.permissions import AllowAny, IsAuthenticated
from accounts.serializers import UserSerializer
from rest_framework.authtoken import views

# Create your views here.


class UserRegistrationView(generics.CreateAPIView):
    """
    Registers a new account in the system.
    """
    permission_classes = [AllowAny]
    serializer_class = UserSerializer


class UserRetrieveUpdateView(generics.RetrieveUpdateAPIView):
    """
    Retrieves or Updates the user's existing account in the system.
    """
    permission_classes = [IsAuthenticated]
    serializer_class = UserSerializer

    def get_object(self):
        return self.request.user


class UserTokenAuthView(views.ObtainAuthToken):
    """
    Retrieves a users authentication token
    """