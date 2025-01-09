from sqlalchemy import Column, String
from sqlalchemy import REAL

from backend.entity.base import Base


class Login(Base):
    __tablename__ = "deepmass_login"
    contact_info = Column(String, primary_key=True, doc="user email")
    login_time = Column(REAL, doc="user login time")