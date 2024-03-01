from sqlalchemy import Column, Integer, String
from sqlalchemy import REAL

from backend.entity.base import Base


class Job(Base):
    __tablename__ = "deepmass_work"
    id = Column(Integer, primary_key=True)
    email = Column(String, doc="user email")
    spectrum_count = Column(Integer, doc="spectrum processed in a work")
    submit_time = Column(REAL, doc="user submit time")
    end_time = Column(REAL, doc="work end time")
    work_duration = Column(REAL, doc="work duration in")
