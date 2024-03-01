import time

from backend.dao.job_dao import JobDAO


class JobService:
    def __init__(self):
        super().__init__()
        self.dao = JobDAO()

    def start_job(self, email):
        start_time = time.time()

        return self.dao.insert_job(email, start_time)

    def end_job(self, id):
        end_time = time.time()
        return self.dao.update_job(id, end_time)


# if __name__ == '__main__':
#     job = JobService()
#     id = job.start_job('aaa@qq.com')
#     job.end_job(id)
