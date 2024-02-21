import random
import smtplib
from email.mime.text import MIMEText

from backend.load_config import *


def send_mail(user, pwd, sender, receiver, content, title):
    mail_host = "smtp.qq.com"
    message = MIMEText(content, "plain", "utf-8")
    message["From"] = sender
    # message["To"] = receiver
    message["To"] = ",".join(receiver)
    message["Subject"] = title
    try:
        smtpObj = smtplib.SMTP()
        smtpObj.connect(mail_host, 25)  # 25 为 SMTP 端口号
        smtpObj.login(user, pwd)
        smtpObj.sendmail(user, receiver, message.as_string())
        print("邮件发送成功")
    except Exception as e:
        print("邮件发送失败")
        raise e


def gen_random_code():
    code = random.sample(list(range(10, 101)), 6)
    code = list(map(lambda x: str(x), code))
    code = "".join(code)
    return code


def send_mail_deepmass(code, mail_receiver):
    email_content = "感谢您使用deepmass ，您的验证码为：%s" % code
    email_title = "DEEPMASS"

    send_mail(
        mail_user, mail_pwd, mail_sender, mail_receiver, email_content, email_title
    )


if __name__ == "__main__":
    # Send_Email(111111)
    code = gen_random_code()
