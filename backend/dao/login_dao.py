# 导入SQLAlchemy模块

from basedao import BaseDao


# 定义UserDao类，继承自BaseDao
class LoginDao(BaseDao):
    def __init__(self):
        super().__init__()
