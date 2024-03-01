import hashlib

from backend.dao.user_dao import UserDAO
from backend.service.code_service import CaptchaService
from backend.service.login_log_service import LoginLogService


class UserService:
    def __init__(self):
        self.dao = UserDAO()

    def auth_login(self, email, password):
        # 对密码取MD5
        password = hashlib.new("md5", password.encode("utf-8")).hexdigest()
        # 验证密码
        login_flag = self.dao.login(email, password)
        # 登录成功，留下记录
        if login_flag:
            LoginLogService().insert_login_log(email)
        return login_flag

    def user_register(self, email, password, name, captcha):
        validate_flag = CaptchaService().validate(email, captcha)
        if validate_flag:
            return {"msg": "验证码错误或已过期"}
        if self.dao.query_email_exist(email):
            return {"msg": "用户已经注册"}
        self.dao.add_user(email, password, name)
        return {"msg": "注册成功"}
