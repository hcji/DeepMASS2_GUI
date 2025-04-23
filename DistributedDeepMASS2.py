# -*- coding: utf-8 -*-
"""
集成本地与分布式识别逻辑的 DeepMASS2 GUI 主程序
"""

# ------------- DistributedDeepMASS2.py ----------------

import shutil
import pickle
from cryptography.hazmat.primitives.asymmetric import padding
from cryptography.hazmat.primitives import serialization, hashes
import os, random, re, socket, struct, uuid, ast, string
from itertools import chain
import pandas as pd
import numpy as np
import matchms.filtering as msfilters
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from molmass import Formula
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt, QThread, QVariant
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import (QApplication, QMainWindow, QFileDialog, QInputDialog,
                             QTableWidgetItem)
from rdkit import Chem
from rdkit.Chem import Draw, rdFMCS
from PyQt5.QtWidgets import QApplication, QGridLayout, QLabel, QMainWindow
# 加密库
from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes
from cryptography.hazmat.backends import default_backend

# UI 类
from uic.main_2 import Ui_MainWindow

# 本地识别函数
from core.importing.load_from_files import load_from_files

# 分布式识别辅助函数
from core.annotating.candidates import search_from_pubchem
from core.annotating.formula import calc_formula_score

# --- 分布式辅助函数 ---
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

# 精准接收固定长度数据
def recv_exact(conn, size):
    data = b""
    while len(data) < size:
        packet = conn.recv(size - len(data))
        if not packet:
            raise ConnectionError("连接中断")
        data += packet
    return data

# 发送加密后的文件
def send_encrypted_file(conn, file_path, sym_key):
    original_filename = os.path.basename(file_path)
    unique_filename = f"{uuid.uuid4()}_{original_filename}"
    print(f"上传的文件编号为：{unique_filename}")
    with open(file_path, "rb") as f:
        file_data = f.read()
    iv, encrypted_data = encrypt_data(file_data, sym_key)
    # print("文件加密完成", encrypt_data)
    payload = iv + encrypted_data
    filename_bytes = unique_filename.encode('utf-8').ljust(128, b'\0')
    filesize = len(payload)
    header = struct.pack('!128sq', filename_bytes, filesize)
    conn.sendall(header)
    conn.sendall(payload)

# 接收服务端返回的加密结果
# def receive_encrypted_result(conn, sym_key):
#     header_size = struct.calcsize('!128sq')
#     header = recv_exact(conn, header_size)
#     # print("接收到 header:", header.hex())
#     _, filesize = struct.unpack('!128sq', header)
#     # print("解析出的 filesize 长度:", filesize)

#     payload = recv_exact(conn, filesize)
#     iv = payload[:16]
#     encrypted_result = payload[16:]
#     # print("接收到加密结果:", encrypted_result.hex())
#     result = decrypt_data(encrypted_result, iv, sym_key)
#     return result.decode('utf-8')

def receive_encrypted_result(conn, sym_key):
    header_size = struct.calcsize('!128sq')
    header = recv_exact(conn, header_size)
    _, payload_size = struct.unpack('!128sq', header)
    payload = recv_exact(conn, payload_size)
    iv, encrypted = payload[:16], payload[16:]
    data = decrypt_data(encrypted, iv, sym_key)
    return data

def extract_blocks(final_result_text):

    results = []  # 用于存储每个块的解析结果
    # 先按照块进行分割
    blocks = final_result_text.strip().split(";\n")
    for block in blocks:
        # 去除多余的空格，并按换行分割成多行
        lines = block.strip().split("\n")
        block_dict = {}
        for line in lines:
            if line.startswith("deepmass_score:"):
                value = line.split("deepmass_score:", 1)[1].strip()
                try:
                    # 尝试转换为字典
                    value = ast.literal_eval(value)
                except Exception as e:
                    print("转换 deepmass_score 为字典失败:", e)
                block_dict["deepmass_score"] = value
            elif line.startswith("reference_spectrum:"):
                value = line.split("reference_spectrum:", 1)[1].strip()
                block_dict["reference_spectrum"] = value
            elif line.startswith("reference_smile:"):
                value = line.split("reference_smile:", 1)[1].strip()
                block_dict["reference_smile"] = value
        if block_dict:
            results.append(block_dict)
    return results

# def identify_unknown_dist(s, database, calc_deepmass_score_result):
def identify_unknown_dist(identify_unknown_s, deepmass_score, reference_spectra):

    # s = search_from_database(s, database, ppm=50)
    s = identify_unknown_s
    candidate = s.get('annotation')
    if candidate is None:
        try:
            s = search_from_pubchem(s, ppm=50)
        except Exception as e:
            print("从 PubChem 查询失败：", e)
            return s
    candidate = s.get('annotation')
    if candidate is None:
        return s
    formula_score = calc_formula_score(s)
    # deepmass_score, reference_spectrum = calc_deepmass_score_result


    for i in candidate.index:
        k = candidate.loc[i, 'InChIKey']
        f = candidate.loc[i, 'MolecularFormula']
        candidate.loc[i, 'Formula Score'] = formula_score.get(f, 0)
        candidate.loc[i, 'Structure Score'] = deepmass_score.get(k, 0)
        candidate.loc[i, 'Consensus Score'] = 0.3 * formula_score.get(f, 0) + 0.7 * deepmass_score.get(k, 0)
    candidate = candidate.sort_values('Consensus Score', ignore_index=True, ascending=False)
    s.set('annotation', candidate)
    # s.set('reference', reference_spectrum)
    s.set('reference', reference_spectra)
    # print(reference_smile,reference_spectrum)
    return s

# def distributed_identification(server_ip, user_service_port, file_path, database):
def distributed_identification(server_ip, user_service_port, file_path):
    # 建立连接
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect((server_ip, user_service_port))

    # 1. 接收服务端动态发送的 RSA 公钥（明文发送）
    pubkey_len_bytes = recv_exact(s, 4)
    pubkey_len = struct.unpack('!I', pubkey_len_bytes)[0]
    server_pub_pem = recv_exact(s, pubkey_len)
    server_public_key = serialization.load_pem_public_key(server_pub_pem, backend=default_backend())
    # print("成功接收并加载服务端公钥。")

    # 2. 生成 AES 对称密钥，并用服务端公钥加密后发送
    sym_key = os.urandom(32)
    encrypted_sym_key = encrypt_sym_key(sym_key, server_public_key)
    s.sendall(struct.pack('!I', len(encrypted_sym_key)))
    s.sendall(encrypted_sym_key)
    # print("成功生成并加密 AES 密钥，已发送给服务端。")

    # 3. 发送加密后的文件
    send_encrypted_file(s, file_path, sym_key)
    print("文件发送完成，等待服务端返回聚合结果...")

    # 4. 接收服务端返回的加密结果并解密
    # result = receive_encrypted_result(s, sym_key)
    # 4. 接收并反序列化结果
    payload_bytes = receive_encrypted_result(s, sym_key)


    s.close()
    remote_results = pickle.loads(payload_bytes)
    print("接收到服务端返回的结果:")
    # print(result)

    # blocks_info = extract_blocks(result)


    # spectras = load_from_files([file_path])
    # result_lists = []
    # # ref_smile_lists = []
    # for i, s_obj in enumerate(spectras):
    #     # 按索引取得对应的块信息
    #     # block = blocks_info[i]
        

    #     # calc_deepmass_score_result = (block.get("deepmass_score"), block.get("reference_spectrum"))
    #     # ref_smile = blocks_info.get("reference_smile")
    #     # 按照目前 identify_unknown 的定义，它需要三个参数：
    #     # 1. 光谱对象 s_obj
    #     # 2. 数据库 database
    #     # 3. 计算深度评分结果 calc_deepmass_score_result
        
    #     # print("updated_spec\n", updated_spec)
    #     ref_smile_list = (block.get("deepmass_score"), block.get("reference_spectrum"), block.get("reference_smile"))
    #     # ref_smile_lists.append(ref_smile_list)
    #     updated_spec = identify_unknown_dist(s_obj, database, ref_smile_list)
    #     result_lists.append(updated_spec)

    # return result_lists
    # spectras = load_from_files([file_path])
    
    aggregated = []
    for i, block in enumerate(remote_results):
    # for i, s_obj in enumerate(spectras):
        # print("block\n", block)
        block = remote_results[i]
        deepmass_score       = block['deepmass_score']
        reference_spectra    = block['reference_spectrum']
        identify_unknown_s   = block['identify_unknown_s']
        # 交给你的 identify_unknown_dist
        updated = identify_unknown_dist(identify_unknown_s, deepmass_score, reference_spectra)
        aggregated.append(updated)
    return aggregated

        


class DeepMASS2(QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None):
        super(DeepMASS2, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("DeepMASS2")
        self.setWindowIcon(QtGui.QIcon("icon/favicon.ico"))
        
        try:
            shutil.rmtree('temp')  
            os.mkdir('temp') 
        except:
            pass
        
        # window
        self.label_logo.setPixmap(QPixmap("icon/logo_deepmass.png"))
        
        # plot      
        self.LabelAnno = QLabel()
        self.gridlayoutAnno = QGridLayout(self.groupBoxAnno)
        self.gridlayoutAnno.addWidget(self.LabelAnno)
        self.LabelRef = QLabel()
        self.gridlayoutRef = QGridLayout(self.groupBoxRef)
        self.gridlayoutRef.addWidget(self.LabelRef)
        
        self.figSpe = MakeFigure(3.6, 2.4, dpi = 300)
        self.figSpe_ntb = NavigationToolbar(self.figSpe, self)
        self.gridlayoutfigSpec = QGridLayout(self.box_spectrum)
        self.gridlayoutfigSpec.addWidget(self.figSpe)
        self.gridlayoutfigSpec.addWidget(self.figSpe_ntb)

        self.server_ip=self.server_port=None
        self.file_path=None
        self.identified_spectrums=[]
        self.current_spectrum=None
        self.current_reference=[]
        # try:
        #     self.default_database=pd.read_csv('data/DeepMassStructureDB-v1.1.csv')
        # except Exception as e:
        #     QtWidgets.QMessageBox.critical(self,'Error',f'Missing database: {e}')
        #     return

        # ——— 一开始就禁用这些按钮 ———
        for b in (self.butt_open, self.butt_run, self.butt_save, self.butt_match):
            b.setEnabled(False)

        self.butt_IP.setEnabled(True)


        self.butt_IP.clicked.connect(self.prompt_server)
        self.butt_open.clicked.connect(self.load_spectrums)
        self.butt_save.clicked.connect(self.save_identification)
        self.butt_run.clicked.connect(self.run_identification)
        self.butt_spectrum.clicked.connect(self.plot_spectrum)
        self.butt_loss.clicked.connect(self.plot_loss)
        self.butt_plotComm.clicked.connect(self.plot_mol_with_highlight)

        self.list_spectrum.itemClicked.connect(self.fill_formula_table)
        self.tab_formula.itemClicked.connect(self.fill_structural_table)
        self.tab_structure.itemClicked.connect(self.fill_reference_table)
        self.tab_reference.itemClicked.connect(self.plot_spectrum)

        self.tab_formula.setSelectionBehavior(QtWidgets.QTableView.SelectRows)
        self.tab_structure.setSelectionBehavior(QtWidgets.QTableView.SelectRows)
        self.tab_reference.setSelectionBehavior(QtWidgets.QTableView.SelectRows)
        
        self.progressBar.setValue(0)
        self.progressBar.setFormat('Ready')

    def WarnMsg(self, Text):
        msg = QtWidgets.QMessageBox()
        msg.resize(550, 200)
        msg.setIcon(QtWidgets.QMessageBox.Warning)
        msg.setText(Text)
        msg.setWindowTitle("Warning")
        msg.exec_()    
    
    
    def ErrorMsg(self, Text):
        msg = QtWidgets.QMessageBox()
        msg.resize(550, 200)
        msg.setIcon(QtWidgets.QMessageBox.Critical)
        msg.setText(Text)
        msg.setWindowTitle("Error")
        msg.exec_()


    def InforMsg(self, Text):
        msg = QtWidgets.QMessageBox()
        msg.resize(550, 200)
        msg.setIcon(QtWidgets.QMessageBox.Information)
        msg.setText(Text)
        msg.setWindowTitle("Information")
        msg.exec_()
    def _set_table_widget(self, widget, data):
        widget.setRowCount(0)
        widget.setRowCount(data.shape[0])
        widget.setColumnCount(data.shape[1])
        widget.setHorizontalHeaderLabels(data.columns)
        widget.setVerticalHeaderLabels(data.index.astype(str))
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                if type(data.iloc[i,j]) == np.float64:
                    item = QtWidgets.QTableWidgetItem()
                    item.setData(Qt.EditRole, QVariant(float(data.iloc[i,j])))
                else:
                    item = QtWidgets.QTableWidgetItem(str(data.iloc[i,j]))
                widget.setItem(i, j, item)

    
    def _set_busy(self):
        for b in (self.butt_open,self.butt_run,self.butt_save,self.butt_IP): b.setDisabled(True)
    def _set_finished(self):
        self.progressBar.setValue(100); self.progressBar.setFormat('Ready')
        for b in (self.butt_open,self.butt_run,self.butt_save,self.butt_IP): b.setDisabled(False)

    def get_formula_mass(self, formula):
        f = Formula(formula)
        return f.isotope.mass

    def prompt_server(self):
        txt,ok=QInputDialog.getText(self,'分布式设置','请输入 IP:Port:')
        if not ok: 
            return
        inp=txt.strip()
        if inp:
            
            try:
                ip,pt=inp.split(':');self.server_ip, self.server_port=ip, int(pt)
                self.butt_IP.setText(f"远程→{ip}:{pt}")
                QtWidgets.QMessageBox.information(self,'Information',f"分布式模式: {ip}:{pt}")
                
                for b in (self.butt_open, self.butt_run, self.butt_save):
                    b.setEnabled(True)
            except:
                QtWidgets.QMessageBox.warning(self,'Warning','格式错误，请输入 IP:Port')

        else:
            parts = inp.split(':')
            if len(parts) != 2 or not parts[1].isdigit():
                QtWidgets.QMessageBox.warning(self, 'Warning', '格式错误，请输入 IP:Port（例如 192.168.1.10:6678）')
            
    def load_spectrums(self):
        self._set_busy()
        files,_=QFileDialog.getOpenFileNames(self,'Load','','MS Files (*.mgf *.mat *.msp)')
        if not files: self._set_finished();return
        self.file_path=files[0]
        specs=load_from_files([self.file_path])
        titles=[s.metadata.get('compound_name','') for s in specs]
        self.spectrums=pd.DataFrame({'title':titles,'spectrum':specs})
        self.list_spectrum.clear();[self.list_spectrum.addItem(t) for t in titles]
        self.list_spectrum.setCurrentRow(0)
        self._set_finished()

    def run_identification(self):
        # ——— 首先清空上一次的结果 ———
        self.identified_spectrums = []
        self.current_spectrum = None
        self.current_reference = []
        # 清空三个表格
        for table in (self.tab_formula, self.tab_structure, self.tab_reference):
            table.setRowCount(0)

        self._set_busy()
        QtWidgets.QMessageBox.information(self,'Information','正在运行中，请稍候...')
        self.progressBar.setValue(0)
        if not self.file_path:
            QtWidgets.QMessageBox.critical(self,'Error','请先上传文件');self._set_finished();return
        if self.server_ip and self.server_port:
            try:
                specs=distributed_identification(self.server_ip,self.server_port,self.file_path)
                self.identified_spectrums=specs
                for i in range(len(specs)):
                    self.progressBar.setValue(int(100*(i+1)/len(specs)))
                    QtWidgets.QApplication.processEvents()
                # 1. Navigator
                if self.list_spectrum.count() > 0:
                    self.list_spectrum.setCurrentRow(0)
                    # 2. Formula Finder
                    self.fill_formula_table()
                    if self.tab_formula.rowCount() > 0:
                        self.tab_formula.setCurrentCell(0, 0)
                        # 3. Structure Finder
                        self.fill_structural_table()
                        if self.tab_structure.rowCount() > 0:
                            self.tab_structure.setCurrentCell(0, 0)
                            # 4. Reference Spectrum
                            self.fill_reference_table()
                            if self.tab_reference.rowCount() > 0:
                                self.tab_reference.setCurrentCell(0, 0)
                                # 5. Spectrum 绘图
                                self.plot_spectrum()

            except Exception as e:
                QtWidgets.QMessageBox.critical(self,'Error',f'分布式出错: {e}')
            finally:
                self._set_finished()
        else:
            self._set_finished()
    def fill_reference_table(self):
        idx = self.list_spectrum.currentRow()
        if idx < 0 or idx >= len(self.identified_spectrums):
            return
        self.current_spectrum = self.identified_spectrums[idx]

        if 'reference' not in self.current_spectrum.metadata.keys():
            self.WarnMsg('No identification result for the selected spectrum')
            return

        self.current_reference = self.current_spectrum.metadata['reference']
        if not self.current_reference:
            self.WarnMsg('No valid reference spectra')
            return
        
        i = self.tab_structure.currentRow()
        header = [self.tab_structure.horizontalHeaderItem(c).text()
                  for c in range(self.tab_structure.columnCount())]
        col_smi = header.index('CanonicalSMILES')
        cell = self.tab_structure.item(i, col_smi)
        if cell is None:
            return
        smi_anno = cell.text()
        annotation = self.current_spectrum.metadata['annotation']
        matches = np.where(annotation['CanonicalSMILES'].values == smi_anno)[0]
        if len(matches) == 0:
            QtWidgets.QMessageBox.warning(self, 'Warning',
                f'在 annotation 中未找到 SMILES: {smi_anno}')
            return
        i = matches[0]
        reference_table = []
        for s in self.current_reference:
            if 'smiles' in s.metadata.keys():
                smiles = s.metadata['smiles']
            else:
                smiles = ''
            if 'compound_name' in s.metadata.keys():
                name = s.metadata['compound_name']
            else:
                name = smiles
            if 'adduct' in s.metadata.keys():
                adduct = s.metadata['adduct']
            else:
                adduct = ''
            if 'parent_mass' in s.metadata.keys():
                parent_mass = s.metadata['parent_mass']
            else:
                parent_mass = ''
            if 'database' in s.metadata.keys():
                ref_database = s.metadata['database']
            else:
                ref_database = ''
            reference_table.append([name, adduct, smiles, parent_mass, ref_database])
        reference_table = pd.DataFrame(reference_table, columns = ['name', 'adduct', 'smiles', 'parent_mass', 'database'])
        self._set_table_widget(self.tab_reference, reference_table)
        self.tab_reference.setCurrentCell(0, 0)
        self.plot_spectrum()
        # self._set_finished()



    def fill_formula_table(self):
        idx=self.list_spectrum.currentRow()
        if idx<0 or idx>=len(self.identified_spectrums):
            return
        spec=self.identified_spectrums[idx]
        self.current_spectrum=spec
        ann=spec.metadata.get('annotation')
        if ann is None or ann.empty:
            QtWidgets.QMessageBox.warning(self,'Warning','No available structures');return
        fm=np.unique(ann['MolecularFormula'])
        ms=[self.get_formula_mass(f) for f in fm]
        diffs=(np.abs(np.array(ms)-float(spec.metadata.get('parent_mass',0)))*1000 if 'parent_mass' in spec.metadata else np.zeros(len(fm)))
        df=pd.DataFrame({'formula':fm,'mass':ms,'error (mDa)':diffs}).sort_values('error (mDa)').reset_index(drop=True)
        self._set_table_widget(self.tab_formula,df)
        self.tab_formula.setCurrentCell(0,0)
        self.fill_structural_table();self.fill_information_table()

    def fill_structural_table(self):
        if self.current_spectrum is None: return
        idx=self.tab_formula.currentRow()
        if idx<0: return
        form=self.tab_formula.item(idx,0).text()
        ann=self.current_spectrum.metadata['annotation']
        sub=ann[ann['MolecularFormula']==form].reset_index(drop=True)
        if sub.empty: QtWidgets.QMessageBox.warning(self,'Warning','No structures for selected formula');return
        self._set_table_widget(self.tab_structure,sub)
        self.tab_structure.setCurrentCell(0,0)
        self.fill_reference_table()


    def fill_information_table(self):
        information = self.current_spectrum.metadata
        keys = [k for k in information.keys() if k in ['compound_name', 'precursor_mz', 'precursor_intensity', 'retention_time', 'inchikey', 
                                                       'formula', 'smiles', 'adduct', 'charge', 'parent_mass', 'ionmode']]
        values = [information[k] for k in keys]
        info_table = pd.DataFrame({'keys':keys, 'values':values})
        self._set_table_widget(self.tab_information, info_table)

    def plot_spectrum(self):
        try:
            i = self.tab_reference.currentRow()
            self.figSpe.PlotSpectrum(self.current_spectrum, self.current_reference[i], loss = False)
        except:
            self.WarnMsg('Cannot plot Spectrum')
        self.plot_mol()

    def plot_loss(self):
        try:
            i = self.tab_reference.currentRow()
            reference = self.current_spectrum.metadata['reference'][i]
            self.figSpe.PlotSpectrum(self.current_spectrum, reference, loss = True)
        except:
            self.WarnMsg('Cannot plot Losses')

    def plot_mol(self, highlight=False):

        i = self.tab_structure.currentRow()
        if i < 0:
            return
        header = [self.tab_structure.horizontalHeaderItem(i).text() for i in range(self.tab_structure.columnCount())]
        try:
           j = list(header).index('CanonicalSMILES')
        except ValueError:
            QtWidgets.QMessageBox.warning(self, 'Warning', '找不到 CanonicalSMILES 列')
            return
        mol_anno = self.tab_structure.item(i, j).text()
        mol_anno = Chem.MolFromSmiles(mol_anno)

        i = self.tab_reference.currentRow()

        if i < 0:
            return
        header = [self.tab_reference.horizontalHeaderItem(i).text() for i in range(self.tab_reference.columnCount())]

        j = list(header).index('smiles')
        smir = self.tab_reference.item(i, j).text()
        mol_ref = Chem.MolFromSmiles(smir)


        if highlight:
            mcs = rdFMCS.FindMCS([mol_anno, mol_ref], bondCompare=rdFMCS.BondCompare.CompareOrderExact,
                                 matchValences = True, ringMatchesRingOnly = True)
            mcs_str = mcs.smartsString
            mcs_mol = Chem.MolFromSmarts(mcs_str)
            allsubs_anno = tuple(chain.from_iterable(mol_anno.GetSubstructMatches(mcs_mol)))
            allsubs_ref = tuple(chain.from_iterable(mol_ref.GetSubstructMatches(mcs_mol)))
        else:
            allsubs_anno = ()
            allsubs_ref = ()
        
        if mol_anno is not None:
            file_name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
            Draw.MolToFile(mol_anno, 'temp/{}.png'.format(file_name), wedgeBonds=False, highlightAtoms=allsubs_anno)
            self.LabelAnno.setPixmap(QPixmap('temp/{}.png'.format(file_name)))
        
        if mol_ref is not None:
            file_name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
            Draw.MolToFile(mol_ref, 'temp/{}.png'.format(file_name), wedgeBonds=False, highlightAtoms=allsubs_ref)
            self.LabelRef.setPixmap(QPixmap('temp/{}.png'.format(file_name)))

    def plot_mol_with_highlight(self):
        self.plot_mol(highlight = True)
        self.InforMsg('Finished')
    def save_identification(self):
        save_dir=QFileDialog.getExistingDirectory(self,'Select Save Directory')
        if not save_dir: 
            return
        for spec in self.identified_spectrums:
            name=re.sub(r'[^\w\- ]+','',spec.metadata.get('compound_name','unknown'))
            ann=spec.metadata.get('annotation')
            if ann is None:
                ann = pd.DataFrame(columns=['Title', 'MolecularFormula', 'CanonicalSMILES', 'InChIKey'])
            df=ann.copy()

            df.to_csv(os.path.join(save_dir,f"{name}.csv"),index=True)
        QtWidgets.QMessageBox.information(self,'Information',f'已保存 {len(self.identified_spectrums)} 个结果到 {save_dir}')

class MakeFigure(FigureCanvas):
    def __init__(self,width=5, height=5, dpi=300):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.fig.subplots_adjust(top=0.95,bottom=0.3,left=0.18,right=0.95)
        super(MakeFigure,self).__init__(self.fig) 
        self.axes = self.fig.add_subplot(111)
        self.axes.spines['bottom'].set_linewidth(0.5)
        self.axes.spines['left'].set_linewidth(0.5)
        self.axes.spines['right'].set_linewidth(0.5)
        self.axes.spines['top'].set_linewidth(0.5)
        self.axes.tick_params(width=0.8,labelsize=3)
        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        
        
    def PlotSpectrum(self, spectrum, reference, loss = False):
        self.axes.cla()
        mz, abunds = spectrum.peaks.mz, spectrum.peaks.intensities
        mz1, abunds1 = reference.peaks.mz, reference.peaks.intensities
        if loss:
            try:
                spectrum = msfilters.add_parent_mass(spectrum)
                spectrum = msfilters.add_losses(spectrum, loss_mz_from=10.0, loss_mz_to=2000.0)
                reference = msfilters.add_parent_mass(reference)
                reference = msfilters.add_losses(reference, loss_mz_from=10.0, loss_mz_to=2000.0)
                mz, abunds = spectrum.losses.mz, spectrum.losses.intensities
                mz1, abunds1 = reference.losses.mz, reference.losses.intensities
            except:
                print('Cannot Plot Losses')
                return
        abunds /= np.max(abunds)
        abunds1 /= np.max(abunds1)
        self.axes.vlines(mz, ymin=0, ymax=abunds, color='r', lw = 0.5)
        self.axes.vlines(mz1, ymin=0, ymax=-abunds1, color='b', lw = 0.5)
        self.axes.axhline(y=0,color='black', lw = 0.5)
        self.axes.set_xlabel('m/z', fontsize = 3.5)
        self.axes.set_ylabel('abundance', fontsize = 3.5)
        self.draw()


if __name__=='__main__':
    import sys
    app=QApplication(sys.argv)
    ui=DeepMASS2()
    ui.show()
    sys.exit(app.exec_())
