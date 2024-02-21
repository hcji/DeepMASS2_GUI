# 导入SQLAlchemy模块
from sqlalchemy import SQLAlchemy

from basedao import BaseDao
from ..entity.user import Work


# 定义UserDao类，继承自BaseDao
class WorkDao(BaseDao):
    def __init__(self, session):
        super().__init__(session, Work)

        # 创建一个db对象，传入Flask应用对象
        db = SQLAlchemy(Work)

        session = db.session()
        # 创建一个WorkDao对象
        work_dao = WorkDao(session)
        # 添加一个新用户
        work_dao.insert(
            Work(
                contact_info="contact_info",
                spectrum_count="spectrum_count",
                submit_time="submit_time",
                end_time="end_time",
                work_duration="work_duration",
                verified="yes",
            )
        )
        # 根据邮箱查找一个用户
        work = BaseDao.select_all("contact_info")
        work = BaseDao.select_by_id("contact_info")
        # 修改用户的验证状态
        BaseDao.update(work, verified="no")
        # 删除用户
        BaseDao.delete(work)
