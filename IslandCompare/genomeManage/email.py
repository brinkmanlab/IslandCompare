from django.core.mail import send_mail
import logging

EMAILSENDER = 'from@example.com'
SENDEMAIL = False

def sendAnalysisCompleteEmail(userEmail, jobId):
    if SENDEMAIL:
        send_mail('IslandCompare - Job Complete', 'Your Job has completed',
                  EMAILSENDER, [userEmail], fail_silently=False)

    logging.info("Sent Email To : "+userEmail)