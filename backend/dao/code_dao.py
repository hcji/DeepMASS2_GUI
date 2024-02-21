# 导入SQLAlchemy模块
# 导入SQLAlchemy模块
from sqlalchemy import SQLAlchemy

from basedao import BaseDao
from ..entity.user import Code


# 定义UserDao类，继承自BaseDao
class CodeDao(BaseDao):
    def __init__(self, session):
        super().__init__(session, Code)

        # 创建一个db对象，传入Flask应用对象
        db = SQLAlchemy(Code)

        session = db.session()
        # 创建一个CodeDao对象
        code_dao = CodeDao(session)
        # 添加一个新用户
        code_dao.insert(
            Code(
                contact_info="contact_info",
                verify_code="verify_code",
                verify_time="verify_time",
                verified="yes",
            )
        )
        # 根据邮箱查找一个用户
        code = BaseDao.select_all("contact_infov")
        code = BaseDao.select_by_id("contact_info")
        # 修改用户的验证状态
        BaseDao.update(code, verified="no")
        # 删除用户
        BaseDao.delete(code)
