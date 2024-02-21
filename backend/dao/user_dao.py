# 导入SQLAlchemy模块
from sqlalchemy import SQLAlchemy

from basedao import BaseDao
from ..entity.user import User


# 定义UserDao类，继承自BaseDao
class UserDao(BaseDao):
    def __init__(self, session):
        super().__init__(session, User)

        db = SQLAlchemy(User)

        session = db.session()
        # 创建一个UserDao对象
        user_dao = UserDao(session)
        # 添加一个新用户
        user_dao.insert(
            User(
                id="id",
                contact_info="contact_info",
                name="name",
                passwd="passwd",
                verified="yes",
            )
        )
        # 根据邮箱查找一个用户
        user = BaseDao.select_all("contact_info")
        user = BaseDao.select_by_id("contact_info")
        # 修改用户的验证状态
        BaseDao.update(user, verified="no")
        # 删除用户
        BaseDao.delete(user)
