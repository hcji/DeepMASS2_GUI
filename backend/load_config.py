import yaml


with open("backend/config/config.yaml", "r", encoding="utf-8") as f:
    GLOBAL_CONFIG = yaml.load(f.read(), Loader=yaml.FullLoader)

# print(type(GLOBAL_CONFIG), GLOBAL_CONFIG)
mail_user = GLOBAL_CONFIG["email"]["mail_user"]
mail_pwd = GLOBAL_CONFIG["email"]["mail_pwd"]
mail_sender = GLOBAL_CONFIG["email"]["mail_sender"]
