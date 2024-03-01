# 导入SQLAlchemy模块
# 导入SQLAlchemy模块
from datetime import datetime

from sqlalchemy import select

from backend.dao.basedao import BaseDao
from backend.entity.captcha import Code


# 定义UserDao类，继承自BaseDao
class CaptchaDAO(BaseDao):
    def __init__(self):
        super().__init__()

    def insert_log(self, email, code):
        # 把验证码存入数据，附带当前时间，email、code、time，实际列名称以数据库为准
        verify_time = datetime.now().timestamp()

        ed_code = Code(contact_info=email, verify_code=code, verify_time=verify_time)
        self.session.query(Code).filter(Code.contact_info == email).delete()
        self.session.add(ed_code)

    def query_captcha_code(self, email):
        res_code = (
            self.session.execute(select(Code).where(Code.contact_info.in_([email])))
            .scalars()
            .one_or_none()
        )
        return res_code



    def delete_captcha(self, email):
        self.session.query(Code).filter_by(contact_info=email).delete()
        self.session.commit()

    def commit(self):
        self.session.commit()
