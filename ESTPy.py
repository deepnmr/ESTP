# -*- coding: utf-8 -*-
# GUI developed by Kyungdoe Han 2021

# Library loading
import os 
import sys
import time
from sys import argv, stderr
from json import dumps, encoder
from PyQt5.QtWidgets import QFileDialog, QMessageBox
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import QObject, pyqtSignal, QEventLoop, QProcess, QThread, QRunnable, QThreadPool, pyqtSlot
from PyQt5.QtGui import QPixmap
from sys import argv, stdout
from estp.est_model_noex_gui import est_model
import json
from time import ctime


# a subwindow instance for image viewer
class subwindow(QtWidgets.QWidget):
    def createWindow(self,WindowWidth,WindowHeight):
       parent=None
       super(subwindow,self).__init__(parent)
       self.resize(WindowWidth,WindowHeight)

class Ui_ESTPmain(object):
    
    def handleStdOut(self):
        data = self.process.readAllStandardOutput().data()
        self.messages.append(data.decode('utf-8'))

    def handleStdErr(self):
        data = self.process.readAllStandardError().data()
        self.messages.append(data.decode('utf-8'))
        
    def Subwindow1(self): # for images
       #self.dirchg = 0
       font = QtGui.QFont()
       font.setPointSize(10)        
       self.mySubwindow = subwindow()
       self.mySubwindow.createWindow(1400,900)
       self.mySubwindow.move(570,50)
       self.gridLayouts = QtWidgets.QGridLayout(self.mySubwindow)
       self.gridLayouts.setObjectName("gridLayout_s")
       self.labels1 = QtWidgets.QLabel(self.mySubwindow)
       self.labels1.setFont(font)
       self.labels1.setObjectName("label_7")
       self.gridLayouts.addWidget(self.labels1, 0, 0, 1, 1)
       self.labels2 = QtWidgets.QLabel(self.mySubwindow)
       self.labels2.setFont(font)
       self.labels2.setObjectName("label_8")
       self.gridLayouts.addWidget(self.labels2, 2, 0, 1, 1)
       self.graphicsViews1 = QtWidgets.QGraphicsView(self.mySubwindow)
       self.graphicsViews1.setObjectName("graphicsView")
       self.gridLayouts.addWidget(self.graphicsViews1, 3, 0, 1, 1)
       self.optionss1 = ('A1', 'A2', 'A3', 'Data A1', 'Data A2', 'Data A3', 'MC A1', 'MC A2', 'MC A3', 'Results')
       self.comboBoxs1 = QtWidgets.QComboBox(self.mySubwindow)
       self.comboBoxs1.addItems(self.optionss1)
       self.comboBoxs1.setObjectName("comboBox_2")
       self.gridLayouts.addWidget(self.comboBoxs1, 1, 0, 1, 1)
       self.pushButtons1 = QtWidgets.QPushButton(self.mySubwindow)
       self.pushButtons1.setObjectName("pushButton_3")
       self.pushButtons1.clicked.connect(self.showGraphs)
       self.gridLayouts.addWidget(self.pushButtons1, 1, 1, 1, 2)
       self.mySubwindow.setWindowTitle("ESTPy v1.0 graph viewer")
       self.mySubwindow.setWindowIcon(QtGui.QIcon('logo.png'))    
       self.mySubwindow.show()

    def Subwindow2(self): # for images
       #self.dirchg = 0
       font = QtGui.QFont()
       font.setPointSize(10)        
       self.mySubwindow1 = subwindow()
       self.mySubwindow1.createWindow(800,600)
       self.mySubwindow1.move(570,50)
       self.resultbrowser = QtWidgets.QTextBrowser(self.mySubwindow1)
       self.resultbrowser.setGeometry(QtCore.QRect(5, 5, 750, 500))
       self.resultbrowser.setObjectName("results")
       self.mySubwindow1.setWindowTitle("ESTPy v1.0 result viewer")
       self.mySubwindow1.setWindowIcon(QtGui.QIcon('logo.png'))    
       self.mySubwindow1.show()

    def setupUi(self, ESTPmain):    
        self.dirchg = 0    
        self.loadFileCheck = 0
        self.dir_l = os.getcwd()
        self.dir_s = os.getcwd()
        self.dirLoad = 0
        self.dirSave = 0
        self.path = os.getcwd()
        self.executeActivate = 0
        font = QtGui.QFont()
        font.setPointSize(10)        
        ESTPmain.setObjectName("ESTPmain")
        ESTPmain.resize(432, 833)
        ESTPmain.move(50,50)
        ESTPmain.setWindowIcon(QtGui.QIcon('logo.png'))    
        self.centralwidget = QtWidgets.QWidget(ESTPmain)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setObjectName("verticalLayout")
        self.scrollArea = QtWidgets.QScrollArea(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.scrollArea.sizePolicy().hasHeightForWidth())
        self.scrollArea.setSizePolicy(sizePolicy)
        self.scrollArea.setMinimumSize(QtCore.QSize(0, 0))
        self.scrollArea.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.scrollArea.setWidgetResizable(False)
        self.scrollArea.setObjectName("scrollArea")
        self.scrollAreaWidgetContents = QtWidgets.QWidget()
        self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 387, 758))
        self.scrollAreaWidgetContents.setObjectName("scrollAreaWidgetContents")

        # Input files
        self.fileList = []
        self.options = ('Open single input file', 'Open multiple input files')
        self.loadInputLabel = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.loadInputLabel.setGeometry(QtCore.QRect(5, 0, 221, 21))
        self.loadInputLabel.setFont(font)
        self.loadInputLabel.setObjectName("loadInputLabel")
        self.loadInputCombo = QtWidgets.QComboBox(self.scrollAreaWidgetContents)
        self.loadInputCombo.setGeometry(QtCore.QRect(5, 20, 251, 27))
        self.loadInputCombo.setFont(font)
        self.loadInputCombo.setObjectName("loadInputCombo")
        self.loadInputCombo.addItems(self.options)
        self.loadInputButton = QtWidgets.QPushButton(self.scrollAreaWidgetContents)
        self.loadInputButton.setGeometry(QtCore.QRect(260, 18, 121, 30))
        self.loadInputButton.setFont(font)
        self.loadInputButton.setObjectName("loadInputButton")
        self.loadInputButton.clicked.connect(self.launchDialog)

        # Line 1
        self.line1 = QtWidgets.QFrame(self.scrollAreaWidgetContents)
        self.line1.setGeometry(QtCore.QRect(8, 50, 371, 5))
        self.line1.setFrameShape(QtWidgets.QFrame.HLine)
        self.line1.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line1.setObjectName("line1")
        
        # Files loaded
        self.filesLabel = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.filesLabel.setGeometry(QtCore.QRect(5, 56, 100, 15))
        self.filesLabel.setFont(font)
        self.filesLabel.setObjectName("filesLabel")
        self.listFilesView = QtWidgets.QListView(self.scrollAreaWidgetContents)
        self.listFilesView.setGeometry(QtCore.QRect(5, 76, 379, 72))
        self.listFilesView.setObjectName("listFilesView")

        # Line 2
        self.line_2 = QtWidgets.QFrame(self.scrollAreaWidgetContents)
        self.line_2.setGeometry(QtCore.QRect(8, 152, 371, 5))
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")

        # Output control
        ## Model output control
        self.outputcontrol = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.outputcontrol.setGeometry(QtCore.QRect(5, 154, 221, 21))
        self.outputcontrol.setFont(font)
        self.outputcontrol.setObjectName("outputcontrol")
        ## output name
        self.labelOutput = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.labelOutput.setGeometry(QtCore.QRect(7, 179, 100, 21))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(9)
        font.setUnderline(True)
        self.labelOutput.setFont(font)
        self.labelOutput.setObjectName("labelOutput")
        ## LineEdit for output name
        self.outputname = QtWidgets.QLineEdit(self.scrollAreaWidgetContents)
        self.outputname.setGeometry(QtCore.QRect(110, 176, 271, 27))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.outputname.setFont(font)
        self.outputname.setObjectName("outputname")
        self.showResults = QtWidgets.QCheckBox(self.scrollAreaWidgetContents)
        self.showResults.setGeometry(QtCore.QRect(10, 207, 121, 25))
        self.showResults.setFont(font)
        self.showResults.setChecked(True)
        self.showResults.setObjectName("showResults")
        self.saveResults = QtWidgets.QCheckBox(self.scrollAreaWidgetContents)
        self.saveResults.setGeometry(QtCore.QRect(140, 207, 281, 25))
        self.saveResults.setFont(font)
        self.saveResults.setChecked(True)
        self.saveResults.setObjectName("saveResults")
        self.MCcheck = QtWidgets.QCheckBox(self.scrollAreaWidgetContents)
        self.MCcheck.setGeometry(QtCore.QRect(263, 207, 50, 25))
        self.MCcheck.setFont(font)
        self.MCcheck.setObjectName("MCcheck")
        self.mcline = QtWidgets.QLineEdit(self.scrollAreaWidgetContents)
        self.mcline.setGeometry(QtCore.QRect(315, 206, 51, 27))
        self.mcline.setFont(font)
        self.mcline.setObjectName("mcline")

        # Line 3
        self.line_3 = QtWidgets.QFrame(self.scrollAreaWidgetContents)
        self.line_3.setGeometry(QtCore.QRect(8, 237, 371, 5))
        self.line_3.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_3.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_3.setObjectName("line_3")

        # Residue control
        self.labelResidue = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.labelResidue.setGeometry(QtCore.QRect(5, 241, 150, 21))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.labelResidue.setFont(font)
        self.labelResidue.setObjectName("labelResidue")
        self.A1check = QtWidgets.QCheckBox(self.scrollAreaWidgetContents)
        self.A1check.setGeometry(QtCore.QRect(10, 263, 41, 25))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.A1check.setFont(font)
        self.A1check.setChecked(True)
        self.A1check.setObjectName("A1check")
        self.A2check = QtWidgets.QCheckBox(self.scrollAreaWidgetContents)
        self.A2check.setGeometry(QtCore.QRect(67, 263, 45, 25))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.A2check.setFont(font)
        self.A2check.setChecked(True)
        self.A2check.setObjectName("A2check")
        self.A3check = QtWidgets.QCheckBox(self.scrollAreaWidgetContents)
        self.A3check.setGeometry(QtCore.QRect(123, 263, 47, 25))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.A3check.setFont(font)
        self.A3check.setChecked(True)
        self.A3check.setObjectName("A3check")
        self.labelMethod = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.labelMethod.setGeometry(QtCore.QRect(260, 240, 61, 21))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.labelMethod.setFont(font)
        self.labelMethod.setObjectName("labelMethod")
        self.methodBox = QtWidgets.QComboBox(self.scrollAreaWidgetContents)
        self.methodBox.setGeometry(QtCore.QRect(260, 260, 101, 27))
        font = QtGui.QFont()
        font.setPointSize(10)        
        self.methodBox.setFont(font)
        self.methodBox.setObjectName("methodBox")
        self.options3 = ('Baldwin', 'Matrix', 'NoEx')
        self.methodBox.addItems(self.options3)

        # Line 4
        self.line_4 = QtWidgets.QFrame(self.scrollAreaWidgetContents)
        self.line_4.setGeometry(QtCore.QRect(10, 290, 371, 5))
        self.line_4.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_4.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_4.setObjectName("line_4")

        # Control variables
        self.labelCV = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.labelCV.setGeometry(QtCore.QRect(5, 295, 221, 21))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.labelCV.setFont(font)
        self.labelCV.setObjectName("labelCV")
        self.kexmaxline = QtWidgets.QLineEdit(self.scrollAreaWidgetContents)
        self.kexmaxline.setGeometry(QtCore.QRect(180, 315, 59, 24))
        font = QtGui.QFont()
        font.setPointSize(9)
        self.kexmaxline.setFont(font)
        self.kexmaxline.setObjectName("kexmaxline")
        self.labelkex = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.labelkex.setGeometry(QtCore.QRect(5, 315, 22, 24))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        font.setBold(False)
        font.setUnderline(True)
        font.setWeight(50)
        self.labelkex.setFont(font)
        self.labelkex.setObjectName("labelkex")
        self.pBmaxline = QtWidgets.QLineEdit(self.scrollAreaWidgetContents)
        self.pBmaxline.setGeometry(QtCore.QRect(180, 346, 59, 24))
        font = QtGui.QFont()
        font.setPointSize(9)
        self.pBmaxline.setFont(font)
        self.pBmaxline.setObjectName("pBmaxline")
        self.pBnstepsline = QtWidgets.QLineEdit(self.scrollAreaWidgetContents)
        self.pBnstepsline.setGeometry(QtCore.QRect(307, 346, 58, 24))
        font = QtGui.QFont()
        font.setPointSize(9)
        self.pBnstepsline.setFont(font)
        self.pBnstepsline.setObjectName("pBnstepsline")
        self.pBmax = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.pBmax.setGeometry(QtCore.QRect(143, 346, 30, 24))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        font.setBold(True)
        font.setUnderline(False)
        font.setWeight(75)
        self.pBmax.setFont(font)
        self.pBmax.setObjectName("pBmax")
        self.pBminline = QtWidgets.QLineEdit(self.scrollAreaWidgetContents)
        self.pBminline.setGeometry(QtCore.QRect(70, 346, 61, 24))
        font = QtGui.QFont()
        font.setPointSize(9)
        self.pBminline.setFont(font)
        self.pBminline.setObjectName("pBminline")
        self.kexminline = QtWidgets.QLineEdit(self.scrollAreaWidgetContents)
        self.kexminline.setGeometry(QtCore.QRect(70, 315, 61, 24))
        font = QtGui.QFont()
        font.setPointSize(9)
        self.kexminline.setFont(font)
        self.kexminline.setObjectName("kexminline")
        self.pBnsteps = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.pBnsteps.setGeometry(QtCore.QRect(253, 346, 47, 24))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        font.setBold(True)
        font.setUnderline(False)
        font.setWeight(75)
        self.pBnsteps.setFont(font)
        self.pBnsteps.setObjectName("pBnsteps")
        self.kexnstepsline = QtWidgets.QLineEdit(self.scrollAreaWidgetContents)
        self.kexnstepsline.setGeometry(QtCore.QRect(307, 315, 58, 24))
        font = QtGui.QFont()
        font.setPointSize(9)
        self.kexnstepsline.setFont(font)
        self.kexnstepsline.setObjectName("kexnstepsline")
        self.kexmin = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.kexmin.setGeometry(QtCore.QRect(34, 315, 31, 24))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        font.setBold(True)
        font.setUnderline(False)
        font.setWeight(75)
        self.kexmin.setFont(font)
        self.kexmin.setObjectName("kexmin")
        self.kexMax = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.kexMax.setGeometry(QtCore.QRect(143, 315, 30, 24))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        font.setBold(True)
        font.setUnderline(False)
        font.setWeight(75)
        self.kexMax.setFont(font)
        self.kexMax.setObjectName("kexMax")
        self.labelpB = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.labelpB.setGeometry(QtCore.QRect(5, 346, 22, 24))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        font.setBold(False)
        font.setUnderline(True)
        font.setWeight(50)
        self.labelpB.setFont(font)
        self.labelpB.setObjectName("labelpB")
        self.kexNsteps = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.kexNsteps.setGeometry(QtCore.QRect(253, 315, 47, 24))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        font.setBold(True)
        font.setUnderline(False)
        font.setWeight(75)
        self.kexNsteps.setFont(font)
        self.kexNsteps.setObjectName("kexNsteps")
        self.pBmin = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.pBmin.setGeometry(QtCore.QRect(34, 346, 31, 24))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        font.setBold(True)
        font.setUnderline(False)
        font.setWeight(75)
        self.pBmin.setFont(font)
        self.pBmin.setObjectName("pBmin")
        self.line_5 = QtWidgets.QFrame(self.scrollAreaWidgetContents)
        self.line_5.setGeometry(QtCore.QRect(10, 374, 371, 5))
        self.line_5.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_5.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_5.setObjectName("line_5")
        self.labelCV_2 = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.labelCV_2.setGeometry(QtCore.QRect(5, 380, 221, 21))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.labelCV_2.setFont(font)
        self.labelCV_2.setObjectName("labelCV_2")

        # Messagebox inherits former quickview and messages
        self.messages = QtWidgets.QTextBrowser(self.scrollAreaWidgetContents)
        self.messages.setGeometry(QtCore.QRect(5, 400, 380, 291))
        self.messages.setObjectName("messages")

        # Line 6
        self.line_6 = QtWidgets.QFrame(self.scrollAreaWidgetContents)
        self.line_6.setGeometry(QtCore.QRect(10, 698, 371, 5))
        self.line_6.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_6.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_6.setObjectName("line_6")

        # Prepare / Execute / Reset
        font = QtGui.QFont()
        font.setPointSize(10)    
        self.pushButton = QtWidgets.QPushButton(self.scrollAreaWidgetContents)
        self.pushButton.setGeometry(QtCore.QRect(2, 707, 111, 30))
        self.pushButton.setFont(font)
        self.pushButton.setObjectName("pushButton")
        self.pushButton.clicked.connect(self.prepare) # prepare 
        self.pushButton_2 = QtWidgets.QPushButton(self.scrollAreaWidgetContents)
        self.pushButton_2.setGeometry(QtCore.QRect(122, 707, 131, 30))
        self.pushButton_2.setFont(font)
        self.pushButton_2.setObjectName("pushButton_2")
        self.pushButton_2.clicked.connect(self.execute) # execute
        self.pushButton_2.setEnabled(False)
        self.pushButton_3 = QtWidgets.QPushButton(self.scrollAreaWidgetContents)
        self.pushButton_3.setGeometry(QtCore.QRect(262, 707, 124, 30))
        self.pushButton_3.setFont(font)
        self.pushButton_3.setObjectName("pushButton_3")
        self.pushButton_3.clicked.connect(self.onClick4) # reset

        # Line 7 & scroll area
        self.line_7 = QtWidgets.QFrame(self.scrollAreaWidgetContents)
        self.line_7.setGeometry(QtCore.QRect(10, 743, 371, 5))
        self.line_7.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_7.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_7.setObjectName("line_7")
        self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        self.verticalLayout.addWidget(self.scrollArea)

        # Menubar
        ESTPmain.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(ESTPmain)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 432, 26))
        self.menubar.setObjectName("menubar")
        self.menuESTP = QtWidgets.QMenu(self.menubar)
        self.menuESTP.setObjectName("menuESTP")
        self.menuHelp = QtWidgets.QMenu(self.menubar)
        self.menuHelp.setObjectName("menuHelp")
        ESTPmain.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(ESTPmain)
        self.statusbar.setObjectName("statusbar")
        ESTPmain.setStatusBar(self.statusbar)

        # Actions for menu
        self.actionHelp = QtWidgets.QAction(ESTPmain)
        self.actionHelp.setObjectName("actionHelp")
        self.actionHelp.triggered.connect(self.help)
        self.actionAbout = QtWidgets.QAction(ESTPmain)
        self.actionAbout.setObjectName("actionAbout")
        self.actionAbout.triggered.connect(self.about)
        self.actionSave_project = QtWidgets.QAction(ESTPmain)
        self.actionSave_project.setObjectName("actionSave_project")
        self.actionSave_project.triggered.connect(self.saveProject)
        self.actionLoad_project = QtWidgets.QAction(ESTPmain)
        self.actionLoad_project.setObjectName("actionLoad_project")
        self.actionLoad_project.triggered.connect(self.loadProject)
        # self.actionSave_input_file = QtWidgets.QAction(ESTPmain)
        # self.actionSave_input_file.setObjectName("actionSave_input_file")
        # self.actionSave_input_file.triggered.connect(self.saveInput)
        # self.actionLoad_input_file = QtWidgets.QAction(ESTPmain)
        # self.actionLoad_input_file.setObjectName("actionLoad_input_file")
        # self.actionLoad_input_file.triggered.connect(self.loadInput)
        #self.actionLoad_graphs = QtWidgets.QAction(ESTPmain)
        #self.actionLoad_graphs.setObjectName("actionLoad_graphs")
        self.actionExit = QtWidgets.QAction(ESTPmain)
        self.actionExit.setObjectName("Exit program")
        self.actionExit.triggered.connect(self.shutDown)
        self.menuESTP.addAction(self.actionSave_project)
        self.menuESTP.addAction(self.actionLoad_project)
        #self.menuESTP.addAction(self.actionSave_input_file)
        #self.menuESTP.addAction(self.actionLoad_input_file)
        #self.menuESTP.addAction(self.actionLoad_graphs)
        self.menuESTP.addAction(self.actionExit)        
        self.menuHelp.addAction(self.actionHelp)
        self.menuHelp.addAction(self.actionAbout)
        self.menubar.addAction(self.menuESTP.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())

        # call subwindow
        self.Subwindow1()

        self.retranslateUi(ESTPmain)
        QtCore.QMetaObject.connectSlotsByName(ESTPmain)

    def saveProject(self):
        path1 = os.getcwd()
        self.dir_s = QFileDialog.getExistingDirectory(None, 'Select a directory or make a directory in this window', path1, QFileDialog.ShowDirsOnly)
        self.dirSave = 1
        print(self.dir_s)
        return self.dir_s
    
    def loadProject(self):
        self.response1 = self.getInput()            
        self.dir_l = os.path.dirname(os.path.abspath(self.response1))
        print(self.dir_l)
        self.dirLoad = 1
        entries = [self.response1]
        model = QtGui.QStandardItemModel()
        self.listFilesView.setModel(model)
        for i in entries:
            item = QtGui.QStandardItem(i)
            model.appendRow(item)                    
        #print(self.fileList)
        if self.response != (''):
            self.fileList.append(entries[0])
            text=open(entries[0]).read()
            self.messages.setPlainText(text)
            self.icount = 1
        self.pushButton.setEnabled(False)
        self.pushButton_2.setEnabled(True)
        QtWidgets.QApplication.processEvents()
        return self.dir_l

    # def saveInput(self):
    #     msgbox = QMessageBox(QMessageBox.Question, "Confirmation", "Save current input file (*.estp)?")
    #     msgbox.addButton(QMessageBox.Yes)
    #     msgbox.addButton(QMessageBox.No)

    # def loadInput(self):
    #     msgbox = QMessageBox(QMessageBox.Question, "Confirmation", "Load previous input file (*.estp)?")
    #     msgbox.addButton(QMessageBox.Yes)
    #     msgbox.addButton(QMessageBox.No)

    def about(self):
        msgbox = QMessageBox(QMessageBox.Warning, "About ESTPy v1.0", "Donghan Lee 2021", buttons = QMessageBox.Ok, parent=None)
        msgbox.setIconPixmap(QPixmap("logo.png"))
        msgbox.exec()

    def help(self):
        msgbox = QMessageBox(QMessageBox.Warning, "Help", "Help list", buttons = QMessageBox.Ok, parent=None)
        msgbox.setIconPixmap(QPixmap("logo.png"))
        msgbox.exec() 

    def shutDown(self, ESTPmain):
        msgbox = QMessageBox(QMessageBox.Question, "Confirmation", "Are you sure you want to close?")
        msgbox.addButton(QMessageBox.Yes)
        msgbox.addButton(QMessageBox.No)
        msgbox.setDefaultButton(QMessageBox.No)
        reply = msgbox.exec()
        if reply != QMessageBox.Yes:
            pass
        else:
            sys.exit()
            QApplication.quit()
        
    def launchDialog(self, ESTPmain):
        option = self.options.index(self.loadInputCombo.currentText())
        if option == 0:
            self.response = self.getSingleInput()
            entries = [self.response]
            model = QtGui.QStandardItemModel()
            self.listFilesView.setModel(model)
            for i in entries:
                item = QtGui.QStandardItem(i)
                model.appendRow(item)                    
            #print(self.fileList)
            if self.response != (''):
                self.fileList.append(entries[0])
                text=open(entries[0]).read()
                self.messages.setPlainText(text)
                self.icount = 1
        elif option == 1:
            self.response = self.getMultipleInput()
            if self.response != ([]):
                entries = [self.response[0], self.response[1]]
                model = QtGui.QStandardItemModel()
                self.listFilesView.setModel(model)
                for i in entries:
                    item = QtGui.QStandardItem(i)
                    model.appendRow(item)  
                self.fileList.append(self.response[0])
                self.fileList.append(self.response[1])
                text=open(entries[0]).read()
                self.messages.append("File loading successful. See Quick View for: ")
                self.messages.append(str(entries[0]))
                self.messages.setPlainText(text)
                self.icount = i            
        else:
            print('Error: please select input files')  
            self.messages.append("Error: please select input files")
            self.icount = 0

    def showGraphs(self):
        #dirchg = 0
        option2 = self.optionss1.index(self.comboBoxs1.currentText())
        path = os.getcwd()
        if self.dir_s != None:
            if self.dir_s != os.getcwd() and self.dirSave == 1:
                path = self.dir_s
        if self.dir_l != None:
            if self.dir_l != os.getcwd() and self.dirLoad == 1:
                path = self.dir_l
        if option2 == 0:  
            pix = QPixmap(path + '/' + self.outputname.text() +'A1_noex.png')
            item = QtWidgets.QGraphicsPixmapItem(pix)
            scene = QtWidgets.QGraphicsScene()
            scene.addItem(item)
            self.graphicsViews1.setScene(scene)
            QtWidgets.QApplication.processEvents()
        elif option2 == 1:
            pix = QPixmap(path + '/' + self.outputname.text() +'A2_noex.png')
            item = QtWidgets.QGraphicsPixmapItem(pix)
            scene = QtWidgets.QGraphicsScene()
            scene.addItem(item)
            self.graphicsViews1.setScene(scene)  
            QtWidgets.QApplication.processEvents()
        elif option2 == 2:
            pix = QPixmap(path + '/' + self.outputname.text() +'A3_noex.png')
            item = QtWidgets.QGraphicsPixmapItem(pix)
            scene = QtWidgets.QGraphicsScene()
            scene.addItem(item)
            self.graphicsViews1.setScene(scene) 
            QtWidgets.QApplication.processEvents()
        elif option2 == 3:
            pix = QPixmap(path + '/' + self.outputname.text() +'_dataA1_noex.png')
            item = QtWidgets.QGraphicsPixmapItem(pix)
            scene = QtWidgets.QGraphicsScene()
            scene.addItem(item)
            self.graphicsViews1.setScene(scene) 
        elif option2 == 4:
            pix = QPixmap(path + '/' + self.outputname.text() +'_dataA2_noex.png')
            item = QtWidgets.QGraphicsPixmapItem(pix)
            scene = QtWidgets.QGraphicsScene()
            scene.addItem(item)
            self.graphicsViews1.setScene(scene) 
        elif option2 == 5:
            pix = QPixmap(path + '/' + self.outputname.text() +'_dataA3_noex.png')
            item = QtWidgets.QGraphicsPixmapItem(pix)
            scene = QtWidgets.QGraphicsScene()
            scene.addItem(item)
            self.graphicsViews1.setScene(scene)             
        elif option2 == 6:
            pix = QPixmap(path + '/' + self.outputname.text() +'_MC_A1_noex.png')
            item = QtWidgets.QGraphicsPixmapItem(pix)
            scene = QtWidgets.QGraphicsScene()
            scene.addItem(item)
            self.graphicsViews1.setScene(scene)  
        elif option2 == 7:
            pix = QPixmap(path + '/' + self.outputname.text() +'_MC_A2_noex.png')
            item = QtWidgets.QGraphicsPixmapItem(pix)
            scene = QtWidgets.QGraphicsScene()
            scene.addItem(item)
            self.graphicsViews1.setScene(scene)  
        elif option2 == 8:
            pix = QPixmap(path + '/' + self.outputname.text() +'_MC_A3_noex.png')
            item = QtWidgets.QGraphicsPixmapItem(pix)
            scene = QtWidgets.QGraphicsScene()
            scene.addItem(item)
            self.graphicsViews1.setScene(scene)  
        elif option2 == 9:
            self.Subwindow2()
            file1 = open(path + "/"+self.outputname.text()+'_result.txt', 'r')
            Lines = file1.readlines() 
            file1.close()
            for i in range(len(Lines)):
                self.resultbrowser.append(Lines[i])
            QtWidgets.QApplication.processEvents()
        else:
            print('Error: please select input files')    
            self.messages.append("Error: please select input files") 

    def defImage1(self):
        path = os.getcwd()
        # if dirchg == 1:
        #     path = dir_l
        # elif dirchg == 2:
        #     path = dir_s
        pix = QPixmap(path+'defaultA1_noex.png')
        item = QtWidgets.QGraphicsPixmapItem(pix)
        scene = QtWidgets.QGraphicsScene()
        scene.addItem(item)
        self.graphicsViews1.setScene(scene)  

    def getSingleInput(self):
        file_filter = 'Data File (*.txt)'
        self.response = QFileDialog.getOpenFileName(
            parent= None,
            caption='Select a data file',
            directory=os.getcwd(),
            filter=file_filter,
            initialFilter='Data File (*.txt)'
        )
        print(self.response[0])
        self.loadFileCheck = 1
        return self.response[0]
    
    def getMultipleInput(self):
        file_filter ='Data File (*.txt)'
        self.response = QFileDialog.getOpenFileNames(
            parent=None,
            caption='Select data files',
            directory=os.getcwd(),
            filter=file_filter,
            initialFilter='Data File (*.txt)'
        )
        print(self.response)
        self.loadFileCheck = 2
        return self.response[0]    

    def getInput(self):
        file_filter = 'Data File (*.estp)'
        self.response = QFileDialog.getOpenFileName(
            parent= None,
            caption='Select an input data file',
            directory=os.getcwd(),
            filter=file_filter,
            initialFilter='Data File (*.estp)'
        )
        self.loadFileCheck = 1
        return self.response[0]

    # Execute button
    def execute(self, dirchg):
        
        self.methodBox.setEnabled(False)
        self.A1check.setEnabled(False)
        self.A2check.setEnabled(False)
        self.A3check.setEnabled(False)
        self.kexmaxline.setEnabled(False)
        self.pBmaxline.setEnabled(False)
        self.pBnstepsline.setEnabled(False)
        self.pBminline.setEnabled(False)
        self.kexminline.setEnabled(False)
        self.kexnstepsline.setEnabled(False)
        self.messages.clear() 
        QtWidgets.QApplication.processEvents()
        path = os.getcwd()
 
        if self.dir_s != None:
            if self.dir_s != os.getcwd() and self.dirSave == 1:
                path = self.dir_s
        if self.dir_l != None:
            if self.dir_l != os.getcwd() and self.dirLoad == 1:
                path = self.dir_l
        if self.dirLoad == 1:
           path = self.dir_l
        # if os.name == 'posix':
        #     pathstr = '\'
        # elif os.name == 'nt'
        #     pathstr = '/' 
        file1 = open(path+'/temp.estp', 'r')
        Lines = file1.readlines() 
        file1.close()
        file1 = open(path+'/'+self.outputname.text()+'.estp', 'w')
        for i in range(len(Lines)):
            file1.write(Lines[i])
        file1.close()
        cfgFile = path + '/' + self.outputname.text()+'.estp'
        self.messages.append('Temporary data input created')
        QtWidgets.QApplication.processEvents()
        self.threadpool = QThreadPool()
        worker = executeThread(cfgFile, self.messages, path)
        self.messages.append('#####################')
        self.messages.append('Model finished')
        self.messages.append('#####################')
        QtWidgets.QApplication.processEvents()
        self.defImage1()

    def messages1(self):
        self.messages.clear()    
        self.messages.append("Running normal mode")
        self.messages.append("Please wait for model run...")
            
    # reset
    def onClick4(self):
        self.pushButton.setEnabled(True)
        self.pushButton_2.setEnabled(False)
        self.messages.clear() 
        self.messages.append("### The program has been reset ###")
        self.icount = 0 
        self.fileList = []

    # Prepare input            
    def prepare(self):
        if self.loadFileCheck == 0:
            msgbox1 = QMessageBox()
            msgbox1.setText("No input files loaded")
            msgbox1.exec_()
        else:
            encoder.FLOAT_REPR = lambda o: format(o, '.4f')
            m2 = est_model()
            m2.verbose = True
            newDict = {'Header' : [],
            'Project Name' : self.outputname.text(),
            'init': {
                    'kex':{'min' : float(self.kexminline.text()), 'max' : float(self.kexmaxline.text()), 'nsteps' : int(self.kexnstepsline.text())},
                     'pB':{'min' : float(self.pBminline.text()), 'max' : float(self.pBmaxline.text()), 'nsteps' : int(self.pBnstepsline.text())},
                     'Method' : self.methodBox.currentText()
                    }
            }
            newDict['Header'].append(m2.programName)
            newDict['Header'].append('Time: ' + ctime())
            expList = []
        # if len(argv) < 2:
        #     stderr.write('Please specify the input file(s):\nprepare.py inputfile1 inputfile2 ...\n')
        #     exit()
        # define datasets
            if self.loadFileCheck == 1:
                m2.dataset.addData(self.response)
                expList.append(self.response)
            elif self.loadFileCheck == 2:
                for i in range(0, len(self.response)):
                    m2.dataset.addData(self.response[i])
                    expList.append(self.response[i])
            newDict['datasets'] = expList
            rdlist = m2.dataset.getResidues()
            newDict['residues'] = rdlist
            buf = dumps(newDict, sort_keys=True, indent=4)
            #print(buf)

            path = os.getcwd()
            if self.dir_s != None:
                if self.dir_s != os.getcwd() and self.dirSave == 1:
                    path = self.dir_s
            if self.dir_l != None:
                if self.dir_l != os.getcwd() and self.dirLoad == 1:
                    path = self.dir_l
            file1 = open(path + "/"+self.outputname.text()+'.estp', 'w')
            for i in range(len(buf)):
                file1.write(buf[i])
            file1.close()
            #stderr.write('Done.\n')
            file1 = open(path + "/"+self.outputname.text()+'.estp', 'r')
            Lines = file1.readlines() 
            file1.close()
            file1 = open(path + '/temp.estp', 'w')
            for i in range(len(Lines)):
                file1.write(Lines[i])
            file1.close()
            if self.A1check.isChecked() == False:
                file1 = open(path + '/temp.estp', 'r')
                Lines = file1.readlines() 
                file1.close()
                for i in range(len(Lines)):
                    if "A1" in Lines[i]:
                        Lines[i-1] = '            "flag": "off",\n'
                        file1 = open(path + '/temp.estp', 'w')
                for i in range(len(Lines)):
                    file1.write(Lines[i])
                file1.close()
            #print(Lines)
            self.messages.append("Residue A1 disabled")     
            if self.A2check.isChecked() == False:
                file1 = open(path + '/temp.estp', 'r')
                Lines = file1.readlines() 
                file1.close()
                for i in range(len(Lines)):
                    if "A2" in Lines[i]:
                        Lines[i-1] = '            "flag": "off",\n'
                file1 = open(path + '/temp.estp', 'w')
                for i in range(len(Lines)):
                    file1.write(Lines[i])
                file1.close()
                #print(Lines)
                self.messages.append("Residue A2 disabled")             
            if self.A3check.isChecked() == False:
                file1 = open(path + '/temp.estp', 'r')
                Lines = file1.readlines() 
                file1.close()
                for i in range(len(Lines)):
                    if "A3" in Lines[i]:
                        Lines[i-1] = '            "flag": "off",\n'
                file1 = open(path + '/temp.estp', 'w')
                for i in range(len(Lines)):
                    file1.write(Lines[i])
                file1.close()
                #print(Lines)
                self.messages.append("Residue A3 disabled")         
            option = self.options3.index(self.methodBox.currentText())
            if option == 1:
                file1 = open(path + '/temp.estp', 'r')
                Lines = file1.readlines() 
                file1.close()
                for i in range(len(Lines)):
                    if "Baldwin" in Lines[i] or "NoEx" in Lines[i]:
                        Lines[i] = '        "Method": "Matrix",\n'
                file1 = open(path + '/temp.estp', 'w')
                for i in range(len(Lines)):
                    file1.write(Lines[i])
                file1.close()
                    #print(Lines)
                self.messages.append("Matrix method applied")     
            elif option == 2:
                file1 = open(path + '/temp.estp', 'r')
                Lines = file1.readlines() 
                file1.close()
                for i in range(len(Lines)):
                    if "Baldwin" in Lines[i] or "Matrix" in Lines[i]:
                        Lines[i] = '        "Method": "NoEx",\n'
                file1 = open(path + '/temp.estp', 'w')
                for i in range(len(Lines)):
                    file1.write(Lines[i])
                file1.close()
                #print(Lines)
                self.messages.append("NoEx method applied")                 
                        
            text=open(path +  "/"+ "temp.estp").read()
            self.messages.setPlainText(text) #quickview
            file1 = open(path + '/temp.estp', 'r')
            Lines = file1.readlines() 
            file1.close()
            file1 = open(path + "/"+ self.outputname.text()+'.estp', 'w')
            for i in range(len(Lines)):
                file1.write(Lines[i])
            file1.close()
            self.messages.append("### Input file generated ###")   
            self.pushButton_2.setEnabled(True)

    def retranslateUi(self, ESTPmain):
        _translate = QtCore.QCoreApplication.translate
        ESTPmain.setWindowTitle(_translate("ESTPmain", "ESTPy v1.0"))
        self.pushButton.setText(_translate("ESTPmain", "Prepare input"))
        self.pushButton_2.setText(_translate("ESTPmain", "Execute"))
        self.pushButton_3.setText(_translate("ESTPmain", "Reset"))
        self.loadInputLabel.setText(_translate("ESTPmain", "Input files"))
        self.loadInputButton.setText(_translate("ESTPmain", "Load"))
        self.filesLabel.setText(_translate("ESTPmain", "Files loaded"))
        self.outputcontrol.setText(_translate("ESTPmain", "Model output control"))
        self.labelOutput.setText(_translate("ESTPmain", "Output name"))
        self.outputname.setText(_translate("ESTPmain", "Default"))
        self.showResults.setText(_translate("ESTPmain", "Show results"))
        self.saveResults.setText(_translate("ESTPmain", "Save results"))
        self.MCcheck.setText(_translate("ESTPmain", "MC"))
        self.mcline.setText(_translate("ESTPmain", "3"))
        self.labelResidue.setText(_translate("ESTPmain", "Residue control"))
        self.A1check.setText(_translate("ESTPmain", "A1"))
        self.A2check.setText(_translate("ESTPmain", "A2"))
        self.A3check.setText(_translate("ESTPmain", "A3"))
        self.labelMethod.setText(_translate("ESTPmain", "Method"))
        self.methodBox.setItemText(0, _translate("ESTPmain", "Baldwin"))
        self.methodBox.setItemText(1, _translate("ESTPmain", "Matrix"))
        self.methodBox.setItemText(2, _translate("ESTPmain", "NoEx"))
        self.labelCV.setText(_translate("ESTPmain", "Control variables"))
        self.kexmaxline.setText(_translate("ESTPmain", "400.0"))
        self.labelkex.setText(_translate("ESTPmain", "kex"))
        self.pBmaxline.setText(_translate("ESTPmain", "0.1"))
        self.pBnstepsline.setText(_translate("ESTPmain", "6"))
        self.pBmax.setText(_translate("ESTPmain", "max"))
        self.pBminline.setText(_translate("ESTPmain", "0.01"))
        self.kexminline.setText(_translate("ESTPmain", "10.0"))
        self.pBnsteps.setText(_translate("ESTPmain", "nsteps"))
        self.kexnstepsline.setText(_translate("ESTPmain", "6"))
        self.kexmin.setText(_translate("ESTPmain", "min"))
        self.kexMax.setText(_translate("ESTPmain", "max"))
        self.labelpB.setText(_translate("ESTPmain", "pB"))
        self.kexNsteps.setText(_translate("ESTPmain", "nsteps"))
        self.pBmin.setText(_translate("ESTPmain", "min"))
        self.labelCV_2.setText(_translate("ESTPmain", "Messages"))
        self.menuESTP.setTitle(_translate("ESTPmain", "ESTP"))
        self.menuHelp.setTitle(_translate("ESTPmain", "Help"))
        self.actionHelp.setText(_translate("ESTPmain", "Help"))
        self.actionAbout.setText(_translate("ESTPmain", "About"))
        self.actionSave_project.setText(_translate("ESTPmain", "Save project"))
        self.actionLoad_project.setText(_translate("ESTPmain", "Load project"))
        #self.actionSave_input_file.setText(_translate("ESTPmain", "Save input file"))
        #self.actionLoad_input_file.setText(_translate("ESTPmain", "Load input file"))
        #self.actionLoad_graphs.setText(_translate("ESTPmain", "Load graphs"))
        self.actionExit.setText(_translate("ESTPmain", "Exit program"))
        # subwindow items
        self.labels1.setText(_translate("ESTPmain", "Select a graph"))
        self.labels2.setText(_translate("ESTPmain", "Results"))
        self.pushButtons1.setText(_translate("ESTPmain", "Load"))


class executeThread(QThread):
    finished = pyqtSignal()    # finish signal
    progress = pyqtSignal(int) # progress signal
    def __init__(self, *args, **kw):
        #super(executeThread, self).__init__()
        QThread.__init__(self, parent=None)
        #self.signal = WorkerSignals()
        
        self.run(*args, **kw)
    
    def run(self, fname, messages, path):     
        QtWidgets.QApplication.processEvents()
        self.progress.emit(1) # consider it as 1% increase in the whole progress bar
        m2=est_model()
        m2.verbose = False
        # print( '**************')
        # stdout.write(m2.programName + '\n')
        # print( '**************')
        # define datasets 
        configFile = open(fname)
        conf = json.load(configFile)
        configFile.close()
        projectName = conf['Project Name']
        datasetsNames = conf['datasets']
        for dataset in datasetsNames:
            m2.dataset.addData(dataset)
        residues = conf['residues']
        for r in residues:
            resid = r['name']
            active = r['flag']
            #print( 'Residue: ' + resid + ' ' + active)
            messages.append('Residue: ' + resid + ' ' + active)
            QtWidgets.QApplication.processEvents()
            for r in m2.dataset.res:
                if r.label == resid:
                    if active == 'on':
                        r.active = True
                    elif active == 'off':
                        r.active = False
                    else:
                        print ('Error: wrong flag for the residue' + resid + ' ' + active)
                        exit()
        m2.datapdf(path + '/' + projectName + '_data.pdf')
        p0=m2.initGuessAll(conf['init'], messages)
        out=m2.fit(p0, messages)
        logBuf = m2.getLogBuffer(out, messages)
        file1 = open(path + '/' + projectName + '_result.txt', 'w')
        print (logBuf)
        file1.write(logBuf)
        file1.close()
        m2.pdf(out, path + '/' + projectName + '.pdf', messages)
        print ('########\n')
        self.finished.emit()
        self.progress.emit(100)

