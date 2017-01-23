from django.conf.urls import url
from accounts import views

urlpatterns = [
    url(r'^$', views.UserRetrieveUpdateView.as_view(), name="account"),
    url(r'^register/$', views.UserRegistrationView.as_view(), name="register"),
    url(r'^auth-token/', views.UserTokenAuthView.as_view(), name="auth"),
]
