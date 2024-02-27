from sqlalchemy import Column, Integer, String

from backend.entity.base import Base


class User(Base):
    __tablename__ = "deepmass_user"
    id = Column(Integer, primary_key=True, unique=True, doc="user id")
    name = Column(String, doc="user nick name")
    contact_info = Column(String, unique=True, doc="user email")
    passwd = Column(String, doc="user password with md5")



