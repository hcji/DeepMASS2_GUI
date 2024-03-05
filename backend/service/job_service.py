import logging
import time

from backend.dao.job_dao import JobDAO


class JobService:
    def __init__(self):
        super().__init__()
        self.dao = JobDAO()

    def start_job(self, email):
        logging.info("鉴定任务开始记录开始")
        start_time = time.time()
        id = self.dao.insert_job(email, start_time)
        logging.info("鉴定任务开始记录结束")
        return id

    def end_job(self, id):
        logging.info("鉴定任务完成记录开始")
        end_time = time.time()
        self.dao.update_job(id, end_time)
        logging.info("鉴定任务完成记录结束")


# if __name__ == '__main__':
#     job = JobService()
#     id = job.start_job('aaa@qq.com')
#     job.end_job(id)
