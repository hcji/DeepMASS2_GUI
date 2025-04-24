#!/usr/bin/env python3
import json
import os
import pickle
import socket
import struct
import sys
sys.path.append("..")
import time
from multiprocessing import Manager, Process
import numpy as np
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from sklearn.metrics.pairwise import cosine_similarity
from multiprocessing import Pool
from cryptography.hazmat.primitives import serialization, hashes
from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes
from cryptography.hazmat.backends import default_backend
from cryptography.hazmat.primitives.asymmetric import padding, rsa
def check_inputs(s):
    if s.get('annotation') is None:
        return False
    elif s.get('parent_mass') is None:
        return False
    elif s.get('precursor_mz') is None:
        return False
    else:
        return True
# AES 加密函数
def encrypt_data(data, key):
    iv = os.urandom(16)
    cipher = Cipher(algorithms.AES(key), modes.CFB(iv), backend=default_backend())
    encryptor = cipher.encryptor()
    encrypted = encryptor.update(data) + encryptor.finalize()
    return iv, encrypted

def decrypt_data(encrypted, iv, key):
    cipher = Cipher(algorithms.AES(key), modes.CFB(iv), backend=default_backend())
    decryptor = cipher.decryptor()
    data = decryptor.update(encrypted) + decryptor.finalize()
    return data

# RSA 加解密对称密钥函数
def encrypt_sym_key(sym_key, public_key):
    encrypted = public_key.encrypt(
        sym_key,
        padding.OAEP(
            mgf=padding.MGF1(algorithm=hashes.SHA256()),
            algorithm=hashes.SHA256(),
            label=None
        )
    )
    return encrypted

def decrypt_sym_key(encrypted_sym_key, private_key):
    sym_key = private_key.decrypt(
        encrypted_sym_key,
        padding.OAEP(
            mgf=padding.MGF1(algorithm=hashes.SHA256()),
            algorithm=hashes.SHA256(),
            label=None
        )
    )
    return sym_key


def send_rsa_pub(conn, server_pub_pem):
    conn.sendall(struct.pack('!I', len(server_pub_pem)))
    conn.sendall(server_pub_pem)

def receive_sym_key(conn, server_private_key):
    key_len = struct.unpack('!I', recv_exact(conn, 4))[0]
    encrypted_sym_key = recv_exact(conn, key_len)
    sym_key = decrypt_sym_key(encrypted_sym_key, server_private_key)
    return sym_key


# 生成服务端 RSA 密钥对
server_private_key = rsa.generate_private_key(
    public_exponent=65537,
    key_size=2048,
    backend=default_backend()
)
server_pub_pem = server_private_key.public_key().public_bytes(
    encoding=serialization.Encoding.PEM,
    format=serialization.PublicFormat.SubjectPublicKeyInfo
)


def recv_exact(conn, size):
    data = b""
    while len(data) < size:
        packet = conn.recv(size - len(data))
        if not packet:
            raise ConnectionError("连接中断")
        data += packet
    return data

# 从用户端接收加密文件数据
def receive_encrypted_file(conn, sym_key):
    fileinfo_size = struct.calcsize('!128sq')
    buf = recv_exact(conn, fileinfo_size)
    if not buf:
        return None, None
    filename_bytes, filesize = struct.unpack('!128sq', buf)
    filename = filename_bytes.strip(b'\0').decode('utf-8')
    print(f"服务器接收到文件：{filename}，大小：{filesize}")
    payload = recv_exact(conn, filesize)
    iv = payload[:16]
    encrypted_file = payload[16:]
    # print("接收到加密文件:", encrypted_file.hex())
    file_content = decrypt_data(encrypted_file, iv, sym_key)

    return filename, file_content

# def send_encrypted_result(conn, result_text, sym_key):
#     data = result_text.encode('utf-8')
#     # print("发送加密结果数据:\n", data)
#     iv, encrypted_data = encrypt_data(data, sym_key)
#     # print("加密结果数据:\n", encrypted_data.hex())
#     payload = iv + encrypted_data
#     header = struct.pack('!128sq', b'result'.ljust(128, b'\0'), len(payload))
#     conn.sendall(header)
#     conn.sendall(payload)

def send_encrypted_result(conn, payload_bytes, sym_key):

    iv, encrypted_data = encrypt_data(payload_bytes, sym_key)
    full_payload = iv + encrypted_data
    header = struct.pack('!128sq', b'result'.ljust(128, b'\0'), len(full_payload))
    conn.sendall(header)
    conn.sendall(full_payload)

# def aggregate_results(results):
#     return "\n".join(results)
def check_inputs(s):
    if s.get('annotation') is None:
        return False
    elif s.get('parent_mass') is None:
        return False
    elif s.get('precursor_mz') is None:
        return False
    else:
        return True
def aggregate_client_deepmass_score(agg_data):
    # database = pd.read_csv('/data2/xiongziyao/data/huaweiyun/data/DeepMassStructureDB-v1.1.csv')

    # s_obj = search_from_database(s, database, ppm = 50)
    # if not check_inputs(s_obj):
    #     return None
    # print('s_obj——', s_obj)
    # 如果总数超过300，则根据 distances 排序取前300

    # if len(agg_data['reference_spectrum']) > 100:
    distances = np.array(agg_data['distances'], dtype=float)
    # 使用 np.argsort 对 distances 进行升序排序，返回排序后的索引序列
    idx = np.argsort(distances, kind='stable')[:300]
    agg_data['reference_spectrum'] = [agg_data['reference_spectrum'][i] for i in idx]
    agg_data['reference_smile'] = [agg_data['reference_smile'][i] for i in idx]
    agg_data['reference_vector'] = np.array(agg_data['reference_vector'])[idx]
    agg_data['distances'] = np.array(agg_data['distances'])[idx]
        # distances = np.array(agg_data['distances'], dtype=float)

    s_obj = agg_data['agg_identify_unknown_s']
    if not check_inputs(s_obj):
        return None, None
    # 定义辅助函数：计算分子指纹与相似度
    get_fp = lambda x: AllChem.GetMorganFingerprintAsBitVect(x, radius=2)
    get_sim = lambda x, y: DataStructs.FingerprintSimilarity(x, y)
    get_corr = lambda query, vec: cosine_similarity([query], [vec])[0][0]
    candidate_mol = [Chem.MolFromSmiles(smi) for smi in s_obj.get('annotation')['CanonicalSMILES']]

    # 过滤掉没有 'smiles' 的参考谱
    # ref_spec = [r for r in agg_data['reference_spectrum'] if r.get('smiles') is not None]
    # ref_smile = [r.metadata['smiles'] for r in ref_spec]
    ref_spec = agg_data['reference_spectrum']
    ref_smile = agg_data['reference_smile']
    # 对应的参考向量
    ref_vec = np.array(agg_data['reference_vector'])

    # print('ref_spec', ref_spec)
    # print('ref_smile', ref_smile)
    # print('ref_vec', ref_vec)

    # 计算每个参考谱的指纹
    k, reference_fp = [], []
    for i, smi in enumerate(ref_smile):
        try:
            m = Chem.MolFromSmiles(smi)
            reference_fp.append(get_fp(m))
            k.append(i)
        except Exception as e:
            pass
    if len(k) != len(ref_smile):
        k = np.array(k)
        ref_smile = np.array(ref_smile)[k]
        ref_vec = np.array(ref_vec)[k, :]
    
    # 从上传谱对象 s_obj 中获得候选分子列表
        
    if len(candidate_mol) == 0:
        return None, None
    
    deepmass_score = []
    # 使用聚合数据中的 query_vector（取 agg_data['query_vector']）
    for i in range(len(candidate_mol)):
        try:
            candidate_fp_i = get_fp(candidate_mol[i])
        except Exception as e:
            deepmass_score.append(0)
            continue
        # 计算每个参考谱与 query_vector 的相似度
        candidate_vecsim_i = np.array([get_corr(agg_data['query_vector'], vec) for vec in ref_vec])
        # 计算候选分子与每个参考谱的指纹相似度
        candidate_fpsim_i = np.array([get_sim(candidate_fp_i, ref_fp) for ref_fp in reference_fp])
        # 取指纹相似度最高的20个参考谱
        top20 = np.argsort(-candidate_fpsim_i)[:20]
        candidate_score_i = np.sqrt(np.sum(candidate_vecsim_i[top20] * candidate_fpsim_i[top20]))
        deepmass_score.append(candidate_score_i / 20)
    
    deepmass_score = np.array(deepmass_score)
    deepmass_score /= (np.max(deepmass_score) + 1e-10)
    # 构造结果字典，使用 s_obj.get('annotation')['InChIKey'] 作为键
    result = dict(zip(s_obj.get('annotation')['InChIKey'], deepmass_score))
    # print('rdeepmass_score\n',result)
    return result, ref_spec, s_obj

def send_file_to_client(client_addr, client_port, filename, file_content):
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect((client_addr, client_port))
        
        # 向客户端发送服务端 RSA 公钥
        send_rsa_pub(s, server_pub_pem)
        
        # 接收客户端发送的经过 RSA 加密的 AES 对称密钥
        sym_key = receive_sym_key(s, server_private_key)
        
        # 使用接收到的对称密钥对文件内容进行 AES 加密
        iv, encrypted_data = encrypt_data(file_content, sym_key)
        payload = iv + encrypted_data
        
        fname_bytes = filename.encode('utf-8').ljust(128, b'\0')
        header = struct.pack('!128sq', fname_bytes, len(payload))
        s.sendall(header)
        s.sendall(payload)
        
        # 先读取返回数据的 header，再读取 payload
        header_size = struct.calcsize('!128sq')
        header = recv_exact(s, header_size)
        _, payload_size = struct.unpack('!128sq', header)
        result_payload = recv_exact(s, payload_size)
        
        s.close()
        
        # 从 payload 中提取 IV 和加密数据，并解密
        iv_result = result_payload[:16]
        encrypted_result = result_payload[16:]
        decrypted_result = decrypt_data(encrypted_result, iv_result, sym_key)
        return decrypted_result
    except Exception as e:
        print(f"向客户端 {client_addr}:{client_port} 发送文件失败:", e)
        return b""

def process_client(client, filename, file_content):
    """
    处理单个客户端：发送文件、接收结果、反序列化后返回结果。
    如果发送失败或转换失败，则返回 None
    """
    client_addr, client_file_port = client
    result_bytes = send_file_to_client(client_addr, client_file_port, filename, file_content)
    if not result_bytes:
        print(f"客户端 {client_addr}:{client_file_port} 无法连接，已从注册列表进行删除")
        return None
    try:
        result_data = pickle.loads(result_bytes)
        return result_data
    except Exception as e:
        print("转换客户端结果出错:", e)
        return None

def user_connection_handler(conn, reg_clients):
    try:
        # 1. 发送 RSA 公钥给用户端
        send_rsa_pub(conn, server_pub_pem)
        # 2. 接收加密后的 AES 密钥
        sym_key = receive_sym_key(conn, server_private_key)
        # 3. 接收用户端发送的加密文件
        filename, file_content = receive_encrypted_file(conn, sym_key)
        if filename is None:
            conn.close()
            return
        print(f"已从用户接收到文件 {filename}")

        # 将接收的文件保存到本地
        temp_path = filename
        with open(temp_path, "wb") as f:
            f.write(file_content)

        # 解析上传文件为谱对象，假设 s_obj 为解析后的谱对象
        # spetrums = load_from_files([filename])

        if os.path.exists(temp_path):
            os.remove(temp_path)
        # 一个谱图
        # s_obj = spetrums[0]
        # print('s_obj ', s_obj)
    except Exception as e:
        print("处理用户上传文件异常：", e)
        conn.close()
        return

    print("等待子服务器注册完成...")
    time.sleep(5)
    print("当前注册的子服务器：", list(reg_clients))
    results = []

    # for client in list(reg_clients):
    #     client_addr, client_file_port = client
    #     result_bytes = send_file_to_client(client_addr, client_file_port, filename, file_content)
    #     if not result_bytes:
    #         try:
    #             reg_clients.remove(client)
    #             print(f"客户端 {client_addr}:{client_file_port} 无法连接，已从注册列表进行删除")
    
    clients = list(reg_clients)
    params = [(client, filename, file_content) for client in clients]

    with Pool(processes=len(clients)) as pool:
    # 并行地处理每个客户端
        client_results = pool.starmap(process_client, params)

    for client, result_data in zip(clients, client_results):
        if result_data is None:
            try:
                reg_clients.remove(client)
            except Exception:
                pass
            continue
        try:
            # 使用 pickle.loads() 
            # result_data = pickle.loads(result_bytes)

            # print(result_dict)
            # if isinstance(result_data, list):
            #     selected_result = result_data[0]
            #     # 测试是否成功获取selected_result
            #     # print(selected_result)
            # else:
            #     selected_result = result_data
            # results.append(selected_result)
            results.append(result_data)
            # print("调试：results =", results)

        except Exception as e:
            print("转换客户端结果出错:", e)

    # if not results:
    #     all_results = "无有效客户端结果"
    # else:
    all_results = []  # 用于存储所有结果

    for idx in range(len(results[0])):
        
        # 如果第一个客户端返回的 query_vector 为空，则直接给定默认聚合结果
        if results[0][idx]['query_vector'] is None:
            print(f"警告：谱对象 {idx} 的 query_vector 为空，返回默认结果。")
            deepmass_score, ref_spec, identify_unknown_s = {}, [], results[0][idx]['identify_unknown_s']

        elif results[0][idx]['identify_unknown_s'].get('annotation') is None:
            print(f"警告：谱对象 {idx} 的 annotation 为空，返回默认结果。")
            deepmass_score, ref_spec, identify_unknown_s = {}, [], results[0][idx]['identify_unknown_s']

        elif not check_inputs(results[0][idx]['identify_unknown_s']):
            deepmass_score, ref_spec, identify_unknown_s = {}, [], results[0][idx]['identify_unknown_s']

        else:
            
            agg_query_vector = results[0][idx]['query_vector']
            agg_reference_spectrum = []
            agg_reference_vector = []
            agg_distances = []
            agg_reference_smile = []
            agg_identify_unknown_s = results[0][idx]['identify_unknown_s']
            for r in results:
                
                # print(" r[idx]['reference_smile']\n",  r[idx]['reference_smile'])
                if not r[idx]['reference_smile']:
                    continue
                agg_reference_spectrum.extend(r[idx]['reference_spectrum'])
                agg_reference_vector.extend(r[idx]['reference_vector'])
                agg_distances.extend(r[idx]['distances'])
                agg_reference_smile.extend(r[idx]['reference_smile'])
            agg_data = {
                'query_vector': agg_query_vector,
                'reference_spectrum': agg_reference_spectrum,
                'reference_smile': agg_reference_smile,
                'reference_vector': np.array(agg_reference_vector),
                'distances': np.array(agg_distances),
                'agg_identify_unknown_s': agg_identify_unknown_s
            }
            result = aggregate_client_deepmass_score(agg_data)
            if result is None:
                current_text = f"输入谱对象 {idx} 缺少必要字段，聚合失败。"
                print(current_text)  # 记录日志
                deepmass_score, ref_spec, identify_unknown_s = {}, [], results[0][idx]['identify_unknown_s']
            else:
                deepmass_score, ref_spec, identify_unknown_s = result
        # agg_result_text = f"deepmass_score: {deepmass_score}\nreference_spectrum: {ref_spec}\nidentify_unknown_s: {identify_unknown_s}\n"
        # agg_result_text = (
        #     f"deepmass_score: {deepmass_score}\n"
        #     f"reference_spectrum: {ref_spec}\n"
        #     f"identify_unknown_s: {identify_unknown_s}\n"
        # )
        agg_result_text =({
                "deepmass_score": deepmass_score,
                "reference_spectrum": ref_spec,        # 这里是 List[Spectrum]
                "identify_unknown_s": identify_unknown_s,  # 单个 Spectrum
            })


        # print("agg_result_text\n", agg_result_text)
        # print(f"deepmass_score: {deepmass_score}\nreference_spectrum: {ref_spec}\nidentify_unknown_s: {results[0][idx]['identify_unknown_s']}")
        all_results.append(agg_result_text)



    # final_result_text = ";\n".join(all_results)
    payload = pickle.dumps(all_results)
    send_encrypted_result(conn, payload, sym_key)

    # send_encrypted_result(conn, final_result_text, sym_key)
    conn.shutdown(socket.SHUT_WR)
    conn.close()



def user_service(SERVER_IP, FILE_PORT, reg_clients):
    # 建立一个持久监听的用户服务socket，用于处理每个用户上传连接
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind((SERVER_IP, FILE_PORT))
    s.listen(5)
    print(f"用户上传服务 启动于 {SERVER_IP}:{FILE_PORT}，等待上传文件...")
    while True:
        conn, addr = s.accept()
        # 使用多进程替代线程处理用户连接
        p = Process(target=user_connection_handler, args=(conn, reg_clients))
        # p.daemon = True
        p.start()

# 客户端注册到服务端，并且进行检测
def registration_server(SERVER_IP, REG_PORT, reg_clients):
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind((SERVER_IP, REG_PORT))
    s.listen(10)
    print(f"注册服务器启动于 {SERVER_IP}:{REG_PORT}")
    while True:
        conn, addr = s.accept()
        try:
            data = conn.recv(1024).decode('utf-8')
            reg_info = json.loads(data)
            if reg_info.get("identity") == "client":
                client_tuple = (reg_info["file_addr"], reg_info["file_port"])
                if client_tuple not in reg_clients:
                    reg_clients.append(client_tuple)
                    print(f"注册成功：{reg_info['file_addr']}:{reg_info['file_port']}")
                else:
                    print(f"客户端 {reg_info['file_addr']}:{reg_info['file_port']} 已注册, 无需重复注册")
            else:
                print("收到非客户端注册消息")
        except Exception as e:
            print("注册处理异常：", e)
        finally:
            conn.close()



if __name__ == '__main__':
    # 命令行参数：服务器 IP、注册端口（例如5001）、用户服务端口（例如5002）
    if len(sys.argv) < 4:
        print("用法: python server.py <server_ip> <reg_port> <file_port>")
        sys.exit(1)
    SERVER_IP = sys.argv[1]
    REG_PORT = int(sys.argv[2])
    FILE_PORT = int(sys.argv[3])
    manager = Manager()
    reg_clients = manager.list()   # 存放 (client_ip, client_file_port)
    # 启动注册服务器进程
    reg_proc = Process(target=registration_server, args=(SERVER_IP, REG_PORT, reg_clients))
    reg_proc.daemon = True
    reg_proc.start()
    # 启动用户服务（持续监听用户上传）
    user_service(SERVER_IP, FILE_PORT, reg_clients)

