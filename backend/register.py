from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from sqlalchemy import create_engine
from sqlalchemy.orm import declarative_base

from backend.service.email_service import EmailSenderService
from backend.service.user_service import UserService

engine = create_engine("sqlite:///./User_Information.db", echo=True)
Base = declarative_base()


app = FastAPI()

# 允许所有来源的CORS(解决跨域问题)
origins = ["*"]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


class Register(BaseModel):
    contact_info: str
    vercode: str
    passwd: str
    confirmPassword: str
    name: str
    agreement: str


@app.post("/register")
def register(reg: Register):

    if reg.passwd != reg.confirmPassword:
        return JSONResponse({"msg": "两次输入的密码不一致"})
    
    res = UserService().user_register(
        reg.contact_info, reg.passwd, reg.name, reg.vercode
    )
    return JSONResponse(res)


# 密码重置接口使用的数据模型
class ResetPassword(BaseModel):
    contact_info: str
    vercode: str
    passwd: str
    confirmPassword: str

@app.post("/reset_password")
def reset_password(reg: ResetPassword):
    print("reset_password 接口参数：", reg)
    if reg.passwd != reg.confirmPassword:
        return JSONResponse({"msg": "两次输入的密码不一致"})
    res = UserService().reset_password(
        reg.contact_info, reg.passwd, reg.vercode
    )
    print("reset_password 返回：", res)
    return JSONResponse(res)


@app.get("/sendmail")
def register(email: str):
    EmailSenderService().send_captcha(email)
    return JSONResponse({"msg": "发送成功"})


if __name__ == "__main__":
    import uvicorn

    # app = gr.mount_gradio_app(app, io, path="/gradio")
    uvicorn.run(app, host="0.0.0.0", port=8000)
