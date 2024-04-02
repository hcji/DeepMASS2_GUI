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
    res = UserService().user_register(
        reg.contact_info, reg.passwd, reg.name, reg.vercode
    )
    return JSONResponse(res)


@app.get("/sendmail")
def register(email: str):
    EmailSenderService().send_captcha(email)
    return JSONResponse({"msg": "发送成功"})


if __name__ == "__main__":
    import uvicorn

    # app = gr.mount_gradio_app(app, io, path="/gradio")
    uvicorn.run(app, host="0.0.0.0", port=8000)
