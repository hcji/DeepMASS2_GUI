from backend.service.user_service import UserService


def auth_ps(contact_info, password):
    login_flag = UserService().auth_login(contact_info, password)
    return login_flag


if __name__ == "__main__":
    print(auth_ps("aaa@qq.com", "aaa123"))
