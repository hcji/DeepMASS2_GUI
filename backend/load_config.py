import yaml


with open("backend/config/config.yaml", "r", encoding="utf-8") as f:
    result = yaml.load(f.read(), Loader=yaml.FullLoader)

print(type(result), result)
mail_user = result["email"]["mail_user"]
mail_pwd = result["email"]["mail_pwd"]
mail_sender = result["email"]["mail_sender"]
