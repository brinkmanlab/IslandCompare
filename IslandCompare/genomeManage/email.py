from django.core.mail import send_mail
import logging
from django.conf import settings

def sendAnalysisCompleteEmail(userEmail, jobId):
    if settings.SEND_EMAIL:
        send_mail('IslandCompare - Job Complete', 'Your Job has completed',
                  settings.EMAIL_SENDER, [userEmail], fail_silently=False)

    logging.info("Sent Email To : "+userEmail)