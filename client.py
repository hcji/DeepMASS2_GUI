
#!/usr/bin/env python3
import json
import os
import socket
import struct
import sys
from multiprocessing import Process
import pickle

import hnswlib
import numpy as np
import pandas as pd
from core.annotating.candidates import search_from_database
from core.importing.load_from_files import load_from_files
from gensim.models import Word2Vec

from spec2vec import SpectrumDocument
from spec2vec.vector_operations import calc_vector

from cryptography.hazmat.primitives import serialization, hashes
from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes
from cryptography.hazmat.backends import default_backend
from cryptography.hazmat.primitives.asymmetric import padding

from core.annotating.candidates import search_from_database, search_from_pubchem

# AES 加解密函数
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

# 使用 RSA 公钥加密 AES 对称密钥
def encrypt_sym_key(sym_key, public_key):
    return public_key.encrypt(
        sym_key,
        padding.OAEP(
            mgf=padding.MGF1(algorithm=hashes.SHA256()),
            algorithm=hashes.SHA256(),
            label=None
        )
    )

def recv_exact(conn, size):
    data = b""
    while len(data) < size:
        packet = conn.recv(size - len(data))
        if not packet:
            raise ConnectionError("连接中断")
        data += packet
    return data


# 注册到服务端
def register_with_server(SERVER_IP, REG_PORT, client_ip, CLIENT_FILE_PORT):
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect((SERVER_IP, REG_PORT))
    reg_info = {
        "identity": "client",
        "file_addr": client_ip,
        "file_port": CLIENT_FILE_PORT,
        "message_addr": client_ip,
        # "message_port": CLIENT_MESSAGE_PORT
    }
    s.sendall(json.dumps(reg_info).encode('utf-8'))
    s.close()
    print("已向服务器注册。")

# 接收服务端发送的加密文件
def receive_encrypted_file(conn, sym_key):
    fileinfo_size = struct.calcsize('!128sq')
    buf = recv_exact(conn, fileinfo_size)
    if not buf:
        return None, None
    filename_bytes, filesize = struct.unpack('!128sq', buf)
    filename = filename_bytes.strip(b'\0').decode('utf-8')
    print(f"接收到来自服务端的文件: {filename}，大小: {filesize} 字节")
    payload = recv_exact(conn, filesize)
    iv = payload[:16]
    encrypted_data = payload[16:]
    # print(f"接收到加密数据\n: encrypted_data {encrypted_data}")
    return filename, decrypt_data(encrypted_data, iv, sym_key)

def check_inputs(s):
    if s.get('annotation') is None:
        return False
    elif s.get('parent_mass') is None:
        return False
    elif s.get('precursor_mz') is None:
        return False
    else:
        return True
# 客户端数据处理
def client_calc_deepmass_score(s, p, model, references):
    if not check_inputs(s):
       
        return {
            'query_vector': [],
            'reference_spectrum': [],
            'reference_smile': [],
            'reference_vector': [],
            'distances': [],
            'identify_unknown_s': s
        }
    # 计算查询向量
    query_vector = calc_vector(model, SpectrumDocument(s, n_decimals=2), allowed_missing_percentage=100)
    # print("打印query_vector:", query_vector)
    xq = np.array(query_vector).astype('float32')
    I, D = p.knn_query(xq, 300)
    # print("打印I\n:", I)
    # print("打印D\n:", D)
    # 根据索引从本地 references 中提取参考谱数据
    # reference_spectrum = np.array(references)[I[0, :]]
    # print("打印np.array(references)[I[0,:]]:", reference_spectrum)
    # 根据同样的索引从 p 中获取对应的向量数据

    raw_reference_spectrum = np.array(references)[I[0, :]]
    raw_reference_vector = np.array(p.get_items(I[0, :]))
    raw_distances = D[0, :]
    # print("打印np.array(references)[:10]]:\n", raw_reference_spectrum[:10])
    # 找到有效的索引：只保留那些在 raw_reference_spectrum 中存在 'smiles' 的记录
    valid_indices = [i for i, rec in enumerate(raw_reference_spectrum) if rec.get('smiles') is not None]

    # 同时过滤参考谱、参考向量和距离
    reference_spectrum = raw_reference_spectrum[valid_indices]
    reference_vector = raw_reference_vector[valid_indices]
    distances = raw_distances[valid_indices]
    # 从过滤后的参考谱中提取 SMILES 信息

    reference_smile = [rec.metadata['smiles'] for rec in reference_spectrum]

    result = {
        'query_vector': query_vector,
        'reference_spectrum':reference_spectrum,
        'reference_smile': reference_smile,
        'reference_vector': reference_vector,
        'distances': distances,
        'identify_unknown_s': s
    }
    return result


def client_cal(path, reference_positive_path='data/references_spectrums_positive.pickle',
               reference_negative_path='data/references_spectrums_negative.pickle',
               index_positive_path='data/references_index_positive_spec2vec.bin',
               index_negative_path='data/references_index_negative_spec2vec.bin'):
    deepmass_positive = Word2Vec.load("model/Ms2Vec_allGNPSpositive.hdf5")
    deepmass_negative = Word2Vec.load("model/Ms2Vec_allGNPSnegative.hdf5")

    database = pd.read_csv('data/DeepMassStructureDB-v1.1.csv')
    
    
    spectrums = load_from_files([path])
    # print("打印 load_from_files([path])spectrums:", spectrums)
    if not isinstance(spectrums, list):
        spectrums = [spectrums]
    # print(spectrums)
    # s = search_from_database(s, database, ppm=500)
    # spectrums = [search_from_database(s, database, ppm = 50) for s in spectrums]
    # print("search之后：", spectrums)
    all_client_data = []

    for idx, s in enumerate(spectrums):
        print(f"正在处理索引为{idx}的质谱：", s)

        s = search_from_database(s, database, ppm=50)
        candidate = s.get('annotation')
        if candidate is None:
            try:
                s = search_from_pubchem(s, ppm=50)
            except Exception as e:
                print(f"search_from_pubchem 异常: {e}")
                client_data = {
                    'query_vector': [],
                    'reference_spectrum': [],
                    'reference_smile': [],
                    'reference_vector': [],
                    'distances': [],
                    'identify_unknown_s': s
                }
                all_client_data.append(client_data)
                continue
        candidate = s.get('annotation')
        if candidate is None:
            print("未能获取到 annotation, 返回空结果")
            client_data = {
                'query_vector': [],
                'reference_spectrum': [],
                'reference_smile': [],
                'reference_vector': [],
                'distances': [],
                'identify_unknown_s': s
            }
            # client_data = {
            #     'formula_score': None,
            #     'structure_score': empty_deep_data
            # }
            all_client_data.append(client_data)
            continue
    
        if s.metadata['ionmode'] == 'negative':
            # print("s.metadata['ionmode']为负离子模式")
            p = hnswlib.Index(space='l2', dim=300)
            p.load_index(index_negative_path)
            with open(reference_negative_path, 'rb') as file:
                references = pickle.load(file)
            model = deepmass_negative
            references = np.array(references)
            client_data = client_calc_deepmass_score(s, p, model, references)
            # print("client_data:", client_data)
            all_client_data.append(client_data)
        else:
            p = hnswlib.Index(space='l2', dim=300)
            p.load_index(index_positive_path)
            with open(reference_positive_path, 'rb') as file:
                references = pickle.load(file)
            model = deepmass_positive
            references = np.array(references)
            client_data = client_calc_deepmass_score(s, p, model, references)
            all_client_data.append(client_data)
            # print("client_data:", client_data)
    # print("客户端数据结果：", all_client_data)
    print("客户端处理完成{0}条数据".format(len(all_client_data)))
    # for s in spectrums:
    #     # print("s:",s)
    #     client_data = client_calc_deepmass_score(s, p, model, references)
    #     # print("client_data:", client_data)
    #     all_client_data.append(client_data)
    # print("客户端处理完成" , all_client_data)
    return all_client_data


def process_file(filename, file_content):
    # 保存文件到临时路径
    temp_path = filename
    with open(temp_path, "wb") as f:
        f.write(file_content)
    # 调用 client_cal 对临时文件进行处理
    result_data = client_cal(temp_path,
                             reference_positive_path='data/references_spectrums_positive.pickle',
                             reference_negative_path='data/references_spectrums_negative.pickle',
                             index_positive_path='data/references_index_positive_spec2vec.bin',
                             index_negative_path='data/references_index_negative_spec2vec.bin')
    # 删除临时文件
    if os.path.exists(temp_path):
        os.remove(temp_path)

    # 将结果转换为 JSON 字符串返回给服务器
    pickle_result = pickle.dumps(result_data)
    return pickle_result


# 处理来自服务端的连接
def handle_client_connection(conn):
    try:
        # 1. 动态接收服务端 RSA 公钥
        pubkey_len = struct.unpack('!I', recv_exact(conn, 4))[0]
        server_pub_pem = recv_exact(conn, pubkey_len)
        server_public_key = serialization.load_pem_public_key(server_pub_pem, backend=default_backend())
        # 2. 生成 AES 对称密钥，并用服务端公钥加密后发送给服务端
        sym_key = os.urandom(32)
        encrypted_sym_key = encrypt_sym_key(sym_key, server_public_key)
        conn.sendall(struct.pack('!I', len(encrypted_sym_key)))
        conn.sendall(encrypted_sym_key)
        # 3. 接收服务端发送的加密文件数据
        filename, file_content = receive_encrypted_file(conn, sym_key)
        if filename is None:
            conn.close()
            return
        # 4. 调用处理函数，得到 pickle 序列化后的结果
        pickle_result = process_file(filename, file_content)
        # 5. 对结果进行 AES 加密后发送给服务端
        iv, encrypted_result = encrypt_data(pickle_result, sym_key)
        # print("打印加密后的结果:encrypted_result\n", encrypted_result)
        payload = iv + encrypted_result
        header = struct.pack('!128sq', b'result'.ljust(128, b'\0'), len(payload))
        conn.sendall(header)
        conn.sendall(payload)
        conn.shutdown(socket.SHUT_WR)
    except Exception as e:
        print("子服务器处理文件异常:", e)
    finally:
        conn.close()


def client_file_service(client_ip, CLIENT_FILE_PORT):
    # 建立一个持久监听socket，只绑定一次
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind((client_ip, CLIENT_FILE_PORT))
    s.listen(5)
    print(f"子服务器正在监听文件接收端口 {client_ip}:{CLIENT_FILE_PORT} ...")
    while True:
        conn, addr = s.accept()
        
        p = Process(target=handle_client_connection, args=(conn,))

        p.start()

if __name__ == '__main__':
    if len(sys.argv) < 6:
        print("用法: python client.py <server_ip> <reg_port> <file_port> <client_ip> <client_file_port>")
        sys.exit(1)
    SERVER_IP = sys.argv[1]
    REG_PORT = int(sys.argv[2])
    _ = int(sys.argv[3])
    client_ip = sys.argv[4]
    CLIENT_FILE_PORT = int(sys.argv[5])
    register_with_server(SERVER_IP, REG_PORT, client_ip, CLIENT_FILE_PORT)
    client_file_service(client_ip, CLIENT_FILE_PORT)