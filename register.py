import datetime
import hashlib
import uuid

from fastapi import FastAPI
from fastapi import Form
from fastapi.responses import JSONResponse
from sqlalchemy import create_engine
from sqlalchemy import select
from sqlalchemy.orm import declarative_base
from sqlalchemy.orm import sessionmaker

from send_email import gen_random_code, Send_Email
from user import User, Code

engine = create_engine("sqlite:///./User_Information.db", echo=True)
Base = declarative_base()

app = FastAPI()


@app.post("/register")
async def register(
    name: str = Form(),
    passwd: str = Form(),
    confirmPassword: str = Form(),
    contact_info: str = Form(),
    vercode: str = Form(),
):
    Session = sessionmaker(bind=engine)
    session = Session()

    # TODO 对passwd取md5
    passwd = hashlib.new("md5", passwd.encode("utf-8")).hexdigest()

    # TODO 将用用户输入的验证码与数据库里的进行对比，正确则继续后面的注册流程，错误则提示用户
    # 查询验证码是否存在
    rescode = (
        session.execute(
            select(Code).where(Code.contact_info.in_([contact_info]))
            # .where(Code.verify_code == vercode)
        )
        .scalars()
        .one_or_none()
    )

    # 若验证码正确
    if rescode.verify_code == vercode:
        # 删除验证码
        session.query(Code).filter_by(contact_info=contact_info).delete()
        session.commit()
        # 计算当前时间戳
        now = datetime.datetime.now().timestamp()
        # 对比数据库时间戳是否超时，若超时则返回
        if now - rescode[0].verify_time > 60 * 10:
            return JSONResponse({"code": 202, "msg": "verification code expired"})

        #  查询contact_info是否重复，重复则返回错误信息并在前端进行处理
        results = (
            session.execute(select(User).where(User.contact_info.in_([contact_info])))
            .scalars()
            .all()
        )
        if len(results) > 0:
            return JSONResponse({"code": 420, "msg": "contact_info exists"})
        #  插入一条用户信息到数据库的用户表中，并返回200信息，由前端跳转到登陆界面
        ed_user = User(
            id=uuid.uuid1().hex, name=name, contact_info=contact_info, passwd=passwd
        )
        session.add(ed_user)
        session.commit()

        return JSONResponse({"code": 200, "msg": "注册成功"})
        # return RedirectResponse(url="http://218.245.102.112:12341/")

    else:
        return "验证码错误"


@app.get("/sendmail")
async def register(email: str):
    Session = sessionmaker(bind=engine)
    session = Session()

    code = gen_random_code()
    print(type(email))

    # 把验证码存入数据，附带当前时间，email、code、time，实际列名称以数据库为准

    verify_time = datetime.datetime.now().timestamp()

    ed_code = Code(contact_info=email, verify_code=code, verify_time=verify_time)
    # if len(session.query(Code).filter(Code.contact_info == email).all()) > 0:
    #     print("update")
    session.query(Code).filter(Code.contact_info == email).delete()
    session.add(ed_code)
    Send_Email(code, email)
    session.commit()

    return JSONResponse({"code": 200, "msg": "发送成功"})


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)
