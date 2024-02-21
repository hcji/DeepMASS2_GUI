from sqlalchemy import create_engine
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import declarative_base, sessionmaker

Base = declarative_base()


# 定义一个 BaseDao 类，用于封装增删改查的方法
class BaseDao:
    # 初始化方法，创建一个引擎和一个会话
    def __init__(self):
        self.engine = create_engine("sqlite:///User_Information.db")
        self.session = sessionmaker(bind=self.engine)()

    # 插入数据的方法，接受一个对象作为参数，将其添加到会话并提交
    def insert(self, obj):
        try:
            self.session.add(obj)
            self.session.commit()
            print(f"插入成功：{obj}")
        except SQLAlchemyError as e:
            self.session.rollback()
            print(f"插入失败：{e}")

    # 删除数据的方法，接受一个对象作为参数，将其从会话中删除并提交
    def delete(self, obj):
        try:
            self.session.delete(obj)
            self.session.commit()
            print(f"删除成功：{obj}")
        except SQLAlchemyError as e:
            self.session.rollback()
            print(f"删除失败：{e}")

    # 更新数据的方法，接受一个对象作为参数，将其修改后添加到会话并提交
    def update(self, obj):
        try:
            self.session.add(obj)
            self.session.commit()
            print(f"更新成功：{obj}")
        except SQLAlchemyError as e:
            self.session.rollback()
            print(f"更新失败：{e}")

    # 查询所有数据的方法，接受一个类作为参数，返回一个列表
    def select_all(self, cls):
        try:
            result = self.session.query(cls).all()
            print(f"查询成功：{result}")
            return result
        except SQLAlchemyError as e:
            print(f"查询失败：{e}")
            return None

    # 查询单个数据的方法，接受一个类和一个主键值作为参数，返回一个对象或 None
    def select_by_id(self, cls, id):
        try:
            result = self.session.query(cls).get(id)
            print(f"查询成功：{result}")
            return result
        except SQLAlchemyError as e:
            print(f"查询失败：{e}")
            return None

    # 关闭会话的方法，释放资源
    def close(self):
        self.session.close()
