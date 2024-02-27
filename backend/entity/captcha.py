from sqlalchemy import Column, String
from sqlalchemy import REAL

from backend.entity.base import Base


class Code(Base):
    __tablename__ = "deepmass_code"
    contact_info = Column(String, primary_key=True, doc="user email")
    verify_code = Column(String, doc="user verification code")
    verify_time = Column(REAL, doc="user verification time")