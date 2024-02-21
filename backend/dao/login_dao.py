# 导入SQLAlchemy模块
from sqlalchemy import SQLAlchemy

from basedao import BaseDao
from ..entity.user import Login


# 定义UserDao类，继承自BaseDao
class LoginDao(BaseDao):
    def __init__(self, session):
        super().__init__(session, Login)

        # 创建一个db对象，传入Flask应用对象
        db = SQLAlchemy(Login)

        session = db.session()
        # 创建一个LoginDao对象
        user_dao = LoginDao(session)
        # 添加一个新用户
        user_dao.insert(
            Login(
                id="id",
                contact_info="contact_info",
                login_time="login_time",
                verified="yes",
            )
        )
        # 根据邮箱查找一个用户
        login = BaseDao.select_all("contact_info")
        login = BaseDao.select_by_id("contact_info")
        # 修改用户的验证状态
        BaseDao.update(login, verified="no")
        # 删除用户
        BaseDao.delete(login)
