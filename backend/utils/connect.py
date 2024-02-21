import os

from sqlalchemy import create_engine

engine = create_engine("sqlite:///./backend/sqlite/User_Information.db", echo=True)
# Session = sessionmaker(bind=engine)
# session = Session()
if __name__ == "__main__":
    print(os.getcwd())
    print(engine)
