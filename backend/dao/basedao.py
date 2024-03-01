from sqlalchemy import create_engine
from sqlalchemy.orm import declarative_base, sessionmaker

Base = declarative_base()


class BaseDao:
    def __init__(self):
        self.engine = create_engine(
            "sqlite:///./backend/sqlite/User_Information.db", echo=True
        )
        self.session = sessionmaker(bind=self.engine)()
