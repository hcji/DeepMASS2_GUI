from sqlalchemy import Column, Integer, String
from sqlalchemy import REAL

from backend.entity.base import Base


class Work(Base):
    __tablename__ = "deepmass_work"
    contact_info = Column(String, primary_key=True, doc="user email")
    spectrum_count = Column(Integer, doc="spectrum processed in a work")
    submit_time = Column(REAL, doc="user submit time")
    end_time = Column(REAL, doc="work end time")
    work_duration = Column(REAL, doc="work duration in")