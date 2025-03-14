from datetime import datetime

from backend.dao.captcha_dao import CaptchaDAO


class CaptchaService:
    def __init__(self):
        super().__init__()
        self.dao = CaptchaDAO()

    def validate(self, email: str, code: str) -> bool:
        self.dao.session.expire_all()  # 强制刷新
        res_code = self.dao.query_captcha_code(email)
        now = datetime.now().timestamp()
        # print(f"当前时间：{now}")
        if res_code:
            pass
            # print(f"验证码记录：{res_code.verify_code}, 生成时间：{res_code.verify_time}, 差值：{now - res_code.verify_time}")
        else:
            print("没有查询到验证码记录")
        if res_code is None:
            return True
        if code != res_code.verify_code:
            return True
        if now - res_code.verify_time >= 60 * 10:
            return True
        # 验证通过后删除验证码记录
        self.dao.delete_captcha(email)
        return False

