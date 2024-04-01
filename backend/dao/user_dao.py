# 导入SQLAlchemy模块
import uuid

from sqlalchemy import select

from backend.dao.basedao import BaseDao
from backend.entity.user import User


# 定义UserDao类，继承自BaseDao
class UserDAO(BaseDao):
    def __init__(self):
        super().__init__()

    def login(self, email, password):
        # 在数据库中查询username与passwd是否匹配，username设置为unique
        results = self.session.execute(
            select(User.passwd).where(User.contact_info.in_([email]))
        )
        res = results.scalars().one_or_none()
        if res == password:
            return True
        else:
            return False

    def query_email_exist(self, email):
        obj = self.session.execute(select(User).where(User.contact_info.in_([email])))
        results = obj.scalars().all()
        return len(results) > 0

    def add_user(self, email, password, name):
        user = User(id=uuid.uuid4().hex, name=name, contact_info=email, passwd=password)
        self.session.add(user)
        self.session.commit()
