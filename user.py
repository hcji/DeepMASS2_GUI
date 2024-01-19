from sqlalchemy import REAL
from sqlalchemy.orm import declarative_base

# engine = create_engine('sqlite:///./User_Information.db', echo=True)
Base = declarative_base()

from sqlalchemy import Column, Integer, String


class User(Base):
    __tablename__ = "deepmass_user"
    id = Column(Integer, primary_key=True, unique=True)
    name = Column(String)
    contact_info = Column(String, unique=True)
    passwd = Column(String)

    def _repr_(self):
        return "<User(name='%s', contact_info='%s', passwd='%s')>" % (
            self.name,
            self.contact_info,
            self.passwd,
        )


class Work(Base):
    __tablename__ = "deepmass_work"
    contact_info = Column(String, primary_key=True)
    spectrum_count = Column(Integer)
    submit_time = Column(REAL)
    end_time = Column(REAL)
    work_duration = Column(REAL)


class Login(Base):
    __tablename__ = "deepmass_login"
    contact_info = Column(String, primary_key=True)
    login_time = Column(REAL)

    def _repr_(self):
        return "<login(contact_info='%s', login_time='%s')>" % (
            self.contact_info,
            self.login_time,
        )


class Code(Base):
    __tablename__ = "deepmass_code"
    contact_info = Column(String, primary_key=True)
    verify_code = Column(String)
    verify_time = Column(REAL)

    def _repr_(self):
        return "<code(contact_info='%s', verify_code='%s', verify_time='%s')>" % (
            self.contact_info,
            self.verify_code,
            self.verify_time,
        )
