from datetime import datetime

from backend.dao.captcha_dao import CaptchaDAO


class CaptchaService:
    def __init__(self):
        super().__init__()
        self.dao = CaptchaDAO()

    def validate(self, email: str, code: str) -> bool:
        # 查询验证码
        res_code = self.dao.query_captcha_code(email)

        now = datetime.now().timestamp()
        self.dao.delete_captcha(email)
        # 对比数据库时间戳是否超时，若超时则返回
        if (
            res_code is not None
            and code == res_code.verify_code
            and now - res_code.verify_code < 60 * 10
        ):
            return True
        else:
            return False
