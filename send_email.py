import random
import smtplib
from email.mime.text import MIMEText


def sendMail(user, pwd, sender, receiver, content, title):
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


def Send_Email(code, mail_receiver):
    mail_user = "3571267146@qq.com"
    mail_pwd = "gpzhdwzjwgmrchgd"
    mail_sender = "3571267146@qq.com"
    # shoujian = input('请输入收件人：')
    # csr = input('请输入抄送人：')
    # receivers = csr.split(' ')
    # mail_receiver = receivers + [shoujian]

    email_content = "感谢您使用deepmass ，您的验证码为：%s" % code
    email_title = "DEEPMASS"

    sendMail(
        mail_user, mail_pwd, mail_sender, mail_receiver, email_content, email_title
    )


if __name__ == "__main__":
    # Send_Email(111111)
    code = gen_random_code()
    Send_Email(code, "3571267146@qq.com")
