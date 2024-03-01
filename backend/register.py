import gradio as gr
from fastapi import FastAPI
from fastapi import Form
from fastapi.responses import JSONResponse
from sqlalchemy import create_engine
from sqlalchemy.orm import declarative_base

from backend.gradio_app import demo as io
from backend.service.email_service import EmailSenderService
from backend.service.user_service import UserService

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
    res = UserService().user_register(contact_info, passwd, name, vercode)
    return JSONResponse(res)


@app.get("/sendmail")
def register(email: str):
    EmailSenderService().send_captcha(email)
    return JSONResponse({"code": 200, "msg": "发送成功"})


if __name__ == "__main__":
    import uvicorn

    app = gr.mount_gradio_app(app, io, path="/gradio")
    uvicorn.run(app, host="0.0.0.0", port=8000)
