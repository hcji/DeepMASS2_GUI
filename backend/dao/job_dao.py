# 导入SQLAlchemy模块
from backend.dao.basedao import BaseDao
from backend.entity.job import Job


# 定义UserDao类，继承自BaseDao
class JobDAO(BaseDao):
    def __init__(self):
        super().__init__()

    def insert_job(self, email, start_time):
        job = Job(
            email=email,
            submit_time=start_time,
        )
        self.session.add(job)
        self.session.commit()
        return job.id

    def query_by_id(self, id: int) -> Job:
        job = self.session.query(Job).filter(Job.id == id).first()
        return job

    def update_job(self, id, end_time):
        job = self.query_by_id(id)
        job.end_time = end_time
        job.work_duration = end_time - job.submit_time
        self.session.commit()
