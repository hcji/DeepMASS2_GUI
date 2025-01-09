from backend.dao.login_log_dao import LoginLogDAO


class LoginLogService:
    def __init__(self):
        super().__init__()
        self.dao = LoginLogDAO()

    def insert_login_log(self, email):
        self.dao.insert_log(email)
