from ..dao import code_dao


code_dao_instance = code_dao.CodeDao()
result = code_dao_instance.logindao("a@qq.com", "a123")
