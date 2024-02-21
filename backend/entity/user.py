from sqlalchemy import Column, Integer, String
from sqlalchemy import REAL
from sqlalchemy.orm import declarative_base

Base = declarative_base()


class User(Base):
    __tablename__ = "deepmass_user"
    id = Column(Integer, primary_key=True, unique=True, doc="user id")
    name = Column(String, doc="user nick name")
    contact_info = Column(String, unique=True, doc="user email")
    passwd = Column(String, doc="user password with md5")


class Work(Base):
    __tablename__ = "deepmass_work"
    contact_info = Column(String, primary_key=True, doc="user email")
    spectrum_count = Column(Integer, doc="spectrum processed in a work")
    submit_time = Column(REAL, doc="user submit time")
    end_time = Column(REAL, doc="work end time")
    work_duration = Column(REAL, doc="work duration in")


class Login(Base):
    __tablename__ = "deepmass_login"
    contact_info = Column(String, primary_key=True, doc="user email")
    login_time = Column(REAL, doc="user login time")


class Code(Base):
    __tablename__ = "deepmass_code"
    contact_info = Column(String, primary_key=True, doc="user email")
    verify_code = Column(String, doc="user verification code")
    verify_time = Column(REAL, doc="user verification time")
