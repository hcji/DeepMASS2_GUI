# 导入SQLAlchemy模块
import hashlib
import uuid

from sqlalchemy import select

from backend.dao.basedao import BaseDao
from backend.entity.user import User


# 定义UserDao类，继承自BaseDao
class UserDAO(BaseDao):
    def __init__(self):
        super().__init__()

    def login(self, email, password):
        # 对密码取sha256
        password = hashlib.new("sha256", password.encode("utf-8")).hexdigest()

        # 在数据库中查询username与passwd是否匹配，username设置为unique
        results = self.session.execute(
            select(User.passwd).where(User.contact_info.in_([email]))
        )
        res = results.scalars().one_or_none()
        # print(f"login - 输入密码哈希：{password}, 数据库密码：{res}")

        if res == password:
            return True
        else:
            return False

    def query_email_exist(self, email):
        obj = self.session.execute(select(User).where(User.contact_info.in_([email])))
        results = obj.scalars().all()
        return len(results) > 0

    def add_user(self, email, password, name):
        # 对密码取sha256
        password = hashlib.new("sha256", password.encode("utf-8")).hexdigest()
        user = User(id=uuid.uuid4().hex, name=name, contact_info=email, passwd=password)
        self.session.add(user)
        self.session.commit()

    def update_password(self, email, new_password):
        new_password_hashed = hashlib.new("sha256", new_password.encode("utf-8")).hexdigest()
        # 使用 session.query 获取用户对象
        user = self.session.query(User).filter(User.contact_info == email).first()
        if user:
            # print(f"更新前密码：{user.passwd}")
            # print(f"新密码哈希：{new_password_hashed}")
            user.passwd = new_password_hashed
            self.session.commit()
            # 强制清空 session 缓存，确保下次查询获取的是最新数据
            self.session.expire_all()
            # 再次查询用户以验证更新结果
            updated_user = self.session.query(User).filter(User.contact_info == email).first()
            # print(f"更新后密码：{updated_user.passwd}")

            return True
        else:
            return False

