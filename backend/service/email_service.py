import random
import smtplib
from email.mime.text import MIMEText

from backend.dao.captcha_dao import CaptchaDAO
from backend.load_config import GLOBAL_CONFIG

MAIL_USER = GLOBAL_CONFIG["email"]["mail_user"]
MAIL_PWD = GLOBAL_CONFIG["email"]["mail_pwd"]
MAIL_SENDER = GLOBAL_CONFIG["email"]["mail_sender"]
MAIL_HOST = GLOBAL_CONFIG["email"]["host"]
MAIL_PORT = GLOBAL_CONFIG["email"]["port"]


class EmailSenderService:
    def __init__(self):
        super().__init__()
        self.dao = CaptchaDAO()

    def send_mail(self, receiver,code):
        email_content = "感谢您使用deepmass ，您的验证码为：%s" % code
        email_title = "DEEPMASS"

        message = MIMEText(email_content, "plain", "utf-8")
        message["From"] = MAIL_SENDER
        message["To"] = ",".join(receiver)
        message["Subject"] = email_title
        try:
            smtpObj = smtplib.SMTP()
            smtpObj.connect(MAIL_HOST, 25)  # 25 为 SMTP 端口号
            smtpObj.login(MAIL_USER, MAIL_PWD)
            smtpObj.sendmail(MAIL_USER, receiver, message.as_string())
            print("邮件发送成功")
        except Exception as e:
            print("邮件发送失败")
            raise e

    def send_captcha(self,email):
        code = self.gen_random_code()
        self.dao.insert_log(code, email)
        self.send_mail(email,code)
        self.dao.commit()

    def gen_random_code(self):
        code = random.sample(list(range(10, 101)), 6)
        code = list(map(lambda x: str(x), code))
        code = "".join(code)
        return code
