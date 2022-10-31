# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 10:13:32 2022

@author: DELL
"""

import os
import re
import hnswlib
import string
import random
import shutil 
import numpy as np
import pandas as pd

from itertools import chain
from PyQt5.Qt import QThread
from PyQt5.QtCore import Qt, QVariant
from PyQt5 import QtCore, QtWidgets, QtGui
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import QApplication, QMainWindow, QGridLayout, QLabel

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar

import matchms.filtering as msfilters
from hnswlib import Index
from rdkit import Chem
from rdkit.Chem import Draw, rdFMCS
from rdkit.Chem import DataStructs, AllChem
from molmass import Formula
from matchms.Spectrum import Spectrum
from matchms.importing import load_from_mgf
from ms2deepscore import MS2DeepScore
from ms2deepscore.models import load_model

from uic import main, parameter
from core.identification import identify_unknown
from core.similarity import calc_aligned_similarity


class Parameter(QtWidgets.QWidget, parameter.Ui_Form):
    
    def __init__(self, parent=None): 
        super(Parameter, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("Parameters")
        self.comboBox.addItems(['Jaccard'])


class DeepMASS2(QMainWindow, main.Ui_MainWindow):
    
    def __init__(self, parent=None):
        super(DeepMASS2, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("DeepMASS2")
        
        try:
            shutil.rmtree('temp')  
            os.mkdir('temp') 
        except:
            pass
        
        # window
        self.ParameterUI = Parameter()
        self.ParameterUI.pushButton_OK.clicked.connect(self.set_parameter)
        self.ParameterUI.pushButton_Cancel.clicked.connect(self.ParameterUI.close)
        
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
        
        # parameters
        self.n_ref = 300
        self.n_neb = 20
        self.sim_metric = 'Jaccard'
        self.priority = []
        self.in_silicon_only = True
        
        # data
        self.pn = None
        self.pp = None
        self.database = pd.DataFrame()
        self.spectrums = pd.DataFrame(columns=['title', 'spectrum'])
        self.identified_spectrums = []
        self.reference_positive = None
        self.reference_negative = None
        self.chemical_space = 'biodatabase'
        self.default_database = 'data/MsfinderStructureDB-VS15-plus-GNPS.csv'
        self.default_index_positive = 'data/references_index_positive.bin'
        self.default_index_negative = 'data/references_index_negative.bin'
        self.default_deepmass_positive = 'model/MS2DeepScore_allGNPSpositive.hdf5'
        self.default_deepmass_negative = 'model/MS2DeepScore_allGNPSnegative.hdf5'
        self.default_reference_positive = 'data/references_spectrums_positive.npy'
        self.default_reference_negative = 'data/references_spectrums_negative.npy'
        self.current_spectrum = None
        
        # action
        self.butt_open.setDisabled(True)
        self.butt_open.setIcon(QtGui.QIcon('icon/opened-folder.png'))
        self.butt_run.setDisabled(True)
        self.butt_run.setIcon(QtGui.QIcon('icon/start.png'))
        self.butt_save.setDisabled(True)
        self.butt_save.setIcon(QtGui.QIcon('icon/save.png'))
        self.butt_para.setIcon(QtGui.QIcon('icon/parameters.png'))
        self.butt_open.clicked.connect(self.load_spectrums)
        self.butt_save.clicked.connect(self.save_identification)
        self.butt_run.clicked.connect(self.run_identification)
        self.butt_para.clicked.connect(self.ParameterUI.show)
        self.butt_spectrum.clicked.connect(self.plot_spectrum)
        self.butt_loss.clicked.connect(self.plot_loss)
        self.butt_plotComm.clicked.connect(self.plot_mol_with_highlight)
        self.butt_match.clicked.connect(self.calc_spectrum_matching)
        self.list_spectrum.itemClicked.connect(self.fill_reference_table)
        self.tab_formula.itemClicked.connect(self.fill_structural_table)
        self.tab_reference.itemClicked.connect(self.plot_spectrum)
        self.tab_structure.itemClicked.connect(self.plot_spectrum)
        
        self.tab_formula.setSelectionBehavior(QtWidgets.QTableView.SelectRows)
        self.tab_structure.setSelectionBehavior(QtWidgets.QTableView.SelectRows)
        self.tab_reference.setSelectionBehavior(QtWidgets.QTableView.SelectRows)
        
        # initial
        self.Thread_LoadIndexPositive = None
        self.Thread_LoadIndexNegative = None
        self.Thread_LoadDatabase = None
        self.Thread_LoadReference = None
        self.Thread_Identification = None
        
        self.progressBar.setValue(0)
        self.progressBar.setFormat('Loading database')
        self.Thread_LoadDatabase = Thread_LoadDatabase(self.default_database)
        self.Thread_LoadDatabase._compounds.connect(self._set_database)
        self.Thread_LoadDatabase.start()
        self.Thread_LoadDatabase.finished.connect(self.load_references_positive)
        self.deepmass_positive = MS2DeepScore(load_model("model/MS2DeepScore_allGNPSpositive.hdf5"))
        self.deepmass_negative = MS2DeepScore(load_model("model/MS2DeepScore_allGNPSnegative.hdf5"))


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

      
    def _set_index_positive(self, msg):
        self.pp = msg


    def _set_index_negative(self, msg):
        self.pn = msg

    
    def _set_database(self, msg):
        self.database = msg
        
        
    def _set_reference_positive(self, msg):
        self.reference_positive = msg


    def _set_reference_negative(self, msg):
        self.reference_negative = msg
    
    
    def _set_succeed_annotation(self, msg):
        self.identified_spectrums.append(msg)
        
        
    def _set_process_bar(self, msg):
        self.progressBar.setValue(int(msg))
        

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
        self.butt_open.setDisabled(True)
        self.butt_run.setDisabled(True)
        self.butt_save.setDisabled(True)
        

    def _set_finished(self):
        self.progressBar.setValue(100)
        self.progressBar.setFormat('Ready') 
        self.butt_open.setDisabled(False)
        self.butt_run.setDisabled(False)
        self.butt_save.setDisabled(False)
    
    
    def set_custom_database(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select custom database", "","CSV Files (*.csv)", options=options)
        if fileName:
            self.progressBar.setValue(50)
            self.progressBar.setFormat('Loading database')
            self.Thread_LoadDatabase = Thread_LoadDatabase(fileName)
            self.Thread_LoadDatabase._compounds.connect(self._set_database)
            self.Thread_LoadDatabase.start()
            self.Thread_LoadDatabase.finished.connect(self._set_finished)
            
    
    def set_parameter(self):
        self.n_ref = int(self.ParameterUI.spinBox_nref.value())
        self.n_neb = int(self.ParameterUI.spinBox_nneb.value())
        self.sim_metric = str(self.ParameterUI.comboBox.currentText())
        if self.ParameterUI.radioButton_1.isChecked():
            self.chemical_space = 'biodatabase'
        if self.ParameterUI.radioButton_2.isChecked():
            self.chemical_space = 'pubchem'
        if self.ParameterUI.radioButton_4.isChecked():
            self.chemical_space = 'biodatabase plus'
        if self.ParameterUI.radioButton_3.isChecked():
            self.chemical_space = 'custom'
            self.set_custom_database()
        self.priority = []
        if self.ParameterUI.checkBox_1.isChecked():
            self.priority.append('HMDB')
        if self.ParameterUI.checkBox_2.isChecked():
            self.priority.append('KNApSAcK')
        if self.ParameterUI.checkBox_3.isChecked():
            self.priority.append('ChEBI') 
        if self.ParameterUI.checkBox_4.isChecked():
            self.priority.append('DrugBank')
        if self.ParameterUI.checkBox_5.isChecked():
            self.priority.append('SMPDB')
        if self.ParameterUI.checkBox_6.isChecked():
            self.priority.append('STOFF')
        if self.ParameterUI.checkBox_7.isChecked():
            self.priority.append('YMDB')
        if self.ParameterUI.checkBox_8.isChecked():
            self.priority.append('T3DB')
        if self.ParameterUI.checkBox_9.isChecked():
            self.priority.append('NANPDB')
        if self.ParameterUI.checkBox_10.isChecked():
            self.priority.append('FooDB')
        if self.ParameterUI.checkBox_11.isChecked():
            self.priority.append('Urine')
        if self.ParameterUI.checkBox_12.isChecked():
            self.priority.append('BMDB')
        if self.ParameterUI.checkBox_13.isChecked():
            self.priority.append('LipidMAPS')
        if self.ParameterUI.checkBox_14.isChecked():
            self.priority.append('Saliva')
        if self.ParameterUI.checkBox_15.isChecked():
            self.priority.append('Feces')
        if self.ParameterUI.checkBox_16.isChecked():
            self.priority.append('PlantCyc')
        if self.ParameterUI.checkBox_17.isChecked():
            self.priority.append('ECMDB')
        if self.ParameterUI.checkBox_18.isChecked():
            self.priority.append('Serum')
        if self.ParameterUI.checkBox_19.isChecked():
            self.priority.append('PubChem-bio')
        if self.ParameterUI.checkBox_20.isChecked():
            self.priority.append('CSF')
        if self.ParameterUI.checkBox_21.isChecked():
            self.priority.append('BLEXP')
        if self.ParameterUI.checkBox_22.isChecked():
            self.priority.append('NPA')                                                                                                                                                                                                       
        if self.ParameterUI.checkBox_23.isChecked():
            self.priority.append('COCONUT')
        if self.ParameterUI.checkBox_24.isChecked():
            self.priority.append('UNPD')
        if self.ParameterUI.checkBox_25.isChecked():
            self.priority.append('GNPS')            
        self.ParameterUI.close()
        # print(self.priority)
    
    
    def get_formula_mass(self, formula):
        f = Formula(formula)
        return f.isotope.mass
    
        
    def load_references_positive(self):
        self.progressBar.setValue(22)
        self.progressBar.setFormat('Loading positive references')   
        self.Thread_LoadIndexPositive = Thread_LoadIndex(self.default_index_positive)
        self.Thread_LoadIndexPositive._index.connect(self._set_index_positive)
        self.Thread_LoadIndexPositive.start()
        self.Thread_LoadIndexPositive.finished.connect(self.load_references_negative)


    def load_references_negative(self):
        self.progressBar.setValue(44)
        self.progressBar.setFormat('Loading negative references')   
        self.Thread_LoadIndexNegative = Thread_LoadIndex(self.default_index_negative)           
        self.Thread_LoadIndexNegative._index.connect(self._set_index_negative)
        self.Thread_LoadIndexNegative.start()
        self.Thread_LoadIndexNegative.finished.connect(self.load_reference_spectrums)
        
        
    def load_reference_spectrums(self):
        self.progressBar.setValue(66)
        self.progressBar.setFormat('Loading reference spectrums')
        self.Thread_LoadReference = Thread_LoadReference(self.default_reference_positive, self.default_reference_negative)           
        self.Thread_LoadReference._reference_positive.connect(self._set_reference_positive)
        self.Thread_LoadReference._reference_negative.connect(self._set_reference_negative)
        self.Thread_LoadReference.start()
        self.Thread_LoadReference.finished.connect(self._set_finished)        

        
    def load_spectrums(self):
        self._set_busy()
        self.progressBar.setValue(24)
        self.progressBar.setFormat('Loading data')
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileNames, _ = QtWidgets.QFileDialog.getOpenFileNames(self, "Load", "","MGF Files (*.mgf)", options=options)
        if len(fileNames) == 0:
            self._set_finished()
            return
        spectrums = []
        for fileName in fileNames:
            spectrums += [s for s in load_from_mgf(fileName) if 'compound_name' in list(s.metadata.keys())]
        titles = [s.metadata['compound_name'] for s in spectrums]
        self.spectrums = pd.DataFrame({'title': titles, 'spectrum': spectrums})
        self.set_list_spectrums()
        self._set_finished()


    def set_list_spectrums(self):
        data = self.spectrums
        if len(data) == 0:
            return
        self.list_spectrum.clear()
        for i in data.index:
            self.list_spectrum.addItem(data.loc[i, 'title'])
        self.list_spectrum.show()
        self.list_spectrum.setCurrentRow(0)


    def run_identification(self):
        self._set_busy()
        self.identified_spectrums = []
        if len(self.spectrums) == 0:
            self.ErrorMsg('Please load unknown spectrums first')
            self._set_finished()
            return
        self.progressBar.setValue(0)
        self.progressBar.setFormat('Identifying unknowns')
        spectrums = self.spectrums['spectrum']
        p_positive = self.pp
        p_negative = self.pn
        n_ref = self.n_ref
        n_neb = self.n_neb
        database = self.database
        priority = self.priority
        model_positive = self.deepmass_positive
        model_negative = self.deepmass_negative
        reference_positive = self.reference_positive
        reference_negative = self.reference_negative
        chemical_space = self.chemical_space
        in_silicon_only = self.in_silicon_only
        self.Thread_Identification = Thread_Identification(spectrums, p_positive, p_negative, n_ref, n_neb, database, 
                                                           priority, model_positive, model_negative, 
                                                           reference_positive, reference_negative, chemical_space, 
                                                           in_silicon_only)  
        self.Thread_Identification._result.connect(self._set_succeed_annotation)
        self.Thread_Identification._i.connect(self._set_process_bar)
        self.Thread_Identification.start()
        self.Thread_Identification.finished.connect(self.fill_formula_table())
        

    def fill_reference_table(self):
        data = self.spectrums
        selectItem = self.list_spectrum.currentItem().text()
        w = np.where(data.loc[:, 'title'] == selectItem)[0][0]
        self.current_spectrum = self.identified_spectrums[w]
        if 'reference' not in self.current_spectrum.metadata.keys():
            self.WarnMsg('No identification result for the selected spectrum')
            return
        current_reference = self.current_spectrum.metadata['reference']
        
        reference_table = []
        for s in current_reference:
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
            reference_table.append([name, adduct, smiles, parent_mass])
        reference_table = pd.DataFrame(reference_table, columns = ['name', 'adduct', 'smiles', 'parent_mass'])
        self._set_table_widget(self.tab_reference, reference_table)
        self.tab_reference.setCurrentCell(0, 0)
        self.plot_spectrum()
        self._set_finished()
    
    
    def fill_formula_table(self):
        if 'annotation' not in self.current_spectrum.metadata.keys():
            self.WarnMsg('No available structures')
            return               
        annotation = self.current_spectrum.metadata['annotation']
        if len(annotation) == 0:
            self.WarnMsg('No available structures')
            return            
        formula = np.unique(annotation['MolecularFormula'])
        mass = [self.get_formula_mass(f) for f in formula]
        if 'parent_mass' in self.current_spectrum.metadata.keys():
            diff = np.array([abs(m - self.current_spectrum.metadata['parent_mass']) for m in mass])
        else:
            diff = np.repeat(np.nan, len(mass))
        formula_table = pd.DataFrame({'formula': formula, 'mass': mass, 'error (mDa)': 1000*diff})
        self._set_table_widget(self.tab_formula, formula_table)
        self.tab_formula.setCurrentCell(0, 0)
        self.fill_structural_table()
        self.fill_information_table()
        self._set_finished()
        
    
    def fill_structural_table(self):
        annotation = self.current_spectrum.metadata['annotation']
        i = self.tab_formula.currentRow()
        header = [self.tab_formula.horizontalHeaderItem(i).text() for i in range(self.tab_formula.columnCount())]
        j = list(header).index('formula')
        formula = self.tab_formula.item(i, j).text()
        structural_table = annotation.loc[annotation['MolecularFormula'] == formula,:]
        structural_table = structural_table.reset_index(drop = True)
        if len(structural_table) == 0:
            self.WarnMsg('No available structures')
            return
        self._set_table_widget(self.tab_structure, structural_table)
        self.tab_structure.setCurrentCell(0, 0)
        self.fill_reference_table()
        self._set_finished()
        
        
    def fill_information_table(self):
        information = self.current_spectrum.metadata
        keys = [k for k in information.keys() if type(information[k]) in [str, float, int]]
        values = [information[k] for k in keys]
        info_table = pd.DataFrame({'keys':keys, 'values':values})
        self._set_table_widget(self.tab_information, info_table)
        
    
    def plot_spectrum(self):
        try:
            i = self.tab_reference.currentRow()
            reference = self.current_spectrum.metadata['reference'][i]
            self.figSpe.PlotSpectrum(self.current_spectrum, reference, loss = False)
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
    
    
    def plot_mol(self, highlight = False):
        i = self.tab_structure.currentRow()
        header = [self.tab_structure.horizontalHeaderItem(i).text() for i in range(self.tab_structure.columnCount())]
        j = list(header).index('CanonicalSMILES')
        mol_anno = self.tab_structure.item(i, j).text()
        mol_anno = Chem.MolFromSmiles(mol_anno)
        
        i = self.tab_reference.currentRow()
        header = [self.tab_reference.horizontalHeaderItem(i).text() for i in range(self.tab_reference.columnCount())]
        j = list(header).index('smiles')
        mol_ref = self.tab_reference.item(i, j).text()
        mol_ref = Chem.MolFromSmiles(mol_ref)
        
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
        
        
    def calc_spectrum_matching(self):
        i = self.tab_reference.currentRow()
        spectrum = self.current_spectrum
        reference = self.current_spectrum.metadata['reference'][i]
        
        i = self.tab_structure.currentRow()
        header = [self.tab_structure.horizontalHeaderItem(i).text() for i in range(self.tab_structure.columnCount())]
        j = list(header).index('CanonicalSMILES')
        smiles_anno = self.tab_structure.item(i, j).text()
        smiles_ref = reference.metadata['smiles']
        
        similarity, matching_tab = calc_aligned_similarity(smiles_anno, smiles_ref, spectrum, reference)
        self._set_table_widget(self.tab_matching, matching_tab)
        
        
    def save_identification(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        savePath = QtWidgets.QFileDialog.getExistingDirectory(self, "Save", options=options)
        if savePath:
            if savePath == '':
                self.WarnMsg('Invalid path')
                return
            for s in self.identified_spectrums:
                name = s.metadata['compound_name']
                name = re.sub(r'[^ \w+]', '', name)
                if 'annotation' in s.metadata.keys():
                    annotation = s.metadata['annotation']
                else:
                    annotation = pd.DataFrame(columns=['Title', 'MolecularFormula', 'CanonicalSMILES', 'InChIKey'])
                path = "{}/{}.csv".format(savePath, name)
                annotation.to_csv(path)        
        self.InforMsg('Finished')


class Thread_LoadReference(QThread): 
    _reference_positive = QtCore.pyqtSignal(list)
    _reference_negative = QtCore.pyqtSignal(list)
    
    def __init__(self, positive, negative):
        super().__init__()
        self.positive = positive
        self.negative = negative

    def run(self):       
        reference_positive = np.load(self.positive, allow_pickle=True)
        reference_negative = np.load(self.negative, allow_pickle=True)
        self._reference_positive.emit(list(reference_positive))
        self._reference_negative.emit(list(reference_negative))


class Thread_LoadDatabase(QThread): 
    _compounds = QtCore.pyqtSignal(pd.DataFrame)
    
    def __init__(self, db_path):
        super().__init__()
        self.db_path = db_path

    def run(self):       
        compounds = pd.read_csv(self.db_path)
        self._compounds.emit(compounds)


class Thread_LoadIndex(QThread): 
    _index = QtCore.pyqtSignal(hnswlib.Index)
    
    def __init__(self, spec_path):
        super().__init__()
        self.spec_path = spec_path
        
    def run(self):       
        spec_bin = Index(space = 'l2', dim = 200)
        spec_bin.load_index(self.spec_path)
        self._index.emit(spec_bin)


class Thread_Identification(QThread):
    _i = QtCore.pyqtSignal(int)
    _result = QtCore.pyqtSignal(Spectrum)

    def __init__(self, spectrums, p_positive, p_negative, n_ref, n_neb, database, priority, 
                 model_positive, model_negative, reference_positive, reference_negative,
                 chemical_space, in_silicon_only):
        super(Thread_Identification, self).__init__()
        self.p_positive = p_positive
        self.p_negative = p_negative
        self.n_ref = n_ref
        self.n_neb = n_neb
        self.spectrums = spectrums
        self.database = database
        self.priority = priority
        self.model_positive = model_positive
        self.model_negative = model_negative
        self.reference_positive = reference_positive
        self.reference_negative = reference_negative
        self.chemical_space = chemical_space
        self.in_silicon_only = in_silicon_only

    def __del__(self):
        self.wait()
        self.working = False       

    def run(self):
        for i, s in enumerate(self.spectrums):
            if 'ionmode' in s.metadata.keys():
                if s.metadata['ionmode'] == 'negative':
                    sn = identify_unknown(s, self.p_negative, self.n_ref, self.n_neb, self.database, self.priority, 
                                          self.model_negative, self.reference_negative, self.chemical_space,
                                          self.in_silicon_only)
                else:
                    sn = identify_unknown(s, self.p_positive, self.n_ref, self.n_neb, self.database, self.priority, 
                                          self.model_positive, self.reference_positive, self.chemical_space,
                                          self.in_silicon_only)
            else:
                sn = identify_unknown(s, self.p_positive, self.n_ref, self.n_neb, self.database, self.priority, 
                                          self.model_positive, self.reference_positive, self.chemical_space,
                                          self.in_silicon_only)
            self._i.emit(int(100 * i / len(self.spectrums)))
            self._result.emit(sn)


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
        self.axes.vlines(mz1, ymin = 0, ymax = -abunds1, color='b', lw = 0.5)
        self.axes.axhline(y=0,color='black', lw = 0.5)
        self.axes.set_xlabel('m/z', fontsize = 3)
        self.axes.set_ylabel('abundance', fontsize = 3)
        self.draw()



if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = DeepMASS2()
    ui.show()
    sys.exit(app.exec_())