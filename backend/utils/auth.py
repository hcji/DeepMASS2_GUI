import hashlib
import time

from sqlalchemy import select
from sqlalchemy.orm import Session

from backend.entity.user import Login
from backend.entity.user import User
from backend.utils.connect import engine


# TODO
def auth_ps(contact_info, passwd):
    passwd = hashlib.new("md5", passwd.encode("utf-8")).hexdigest()
    with Session(engine) as session:
        # 在数据库中查询username与passwd是否匹配，username设置为unique
        results = session.execute(
            select(User.passwd).where(User.contact_info.in_([contact_info]))
        )
        res = results.scalars().one_or_none()
        if res == passwd:
            login_time = time.time()
            log = Login(contact_info=contact_info, login_time=login_time)
            session.add(log)
            session.commit()
            return True
        else:
            return False
        # raise gr.Error("Error in authentication process")


if __name__ == "__main__":
    print(auth_ps("aaa@qq.com", "aaa123"))
