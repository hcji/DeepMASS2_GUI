# 导入SQLAlchemy模块
import time

from backend.dao.basedao import BaseDao
from backend.entity.login import Login


# 定义UserDao类，继承自BaseDao
class LoginLogDAO(BaseDao):
    def __init__(self):
        super().__init__()

    def insert_log(self, email):
        login_time = time.time()
        log = Login(contact_info=email, login_time=login_time)
        self.session.add(log)
        self.session.commit()
