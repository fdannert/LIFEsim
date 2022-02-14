from PyQt5.QtCore import (Qt, QRegExp)
from PyQt5.QtWidgets import (QLabel, QLineEdit, QWidget, QHBoxLayout, QDoubleSpinBox, QFileDialog,
                             QPushButton, QButtonGroup, QRadioButton, QVBoxLayout)
from PyQt5.QtGui import (QDoubleValidator,QRegExpValidator)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure


class RadioButtonWidget(QWidget):
    def __init__(self, *args, **kwargs):
        super(RadioButtonWidget, self).__init__(*args, **kwargs)
        hlayout = QHBoxLayout(self)
        vlayout = QVBoxLayout()
        bg = QButtonGroup(self)

        self.label = QLabel('Array Baseline')

        self.bl_HZ = QRadioButton("Optimize for eHZ", self)
        self.bl_HZ.setChecked(True)

        self.bl_pl = QRadioButton("Optimize for Planet", self)

        bg.addButton(self.bl_HZ)
        bg.addButton(self.bl_pl)

        vlayout.addWidget(self.bl_HZ)
        vlayout.addWidget(self.bl_pl)
        vlayout.setAlignment(Qt.AlignRight)

        hlayout.addWidget(self.label)
        hlayout.addLayout(vlayout)

class DoubleBoxLabel(QWidget):
    def __init__(self, label, mini, maxi, step, value, suffix, *args, **kwargs):
        super(DoubleBoxLabel, self).__init__(*args, **kwargs)
        self.box = QLineEdit()
        self.label = QLabel(text=label)

        self.val = QDoubleValidator(bottom=mini,
                                    top=maxi)

        rx = QRegExp()
        rx.setPattern('\d*\.?\d*')
        self.val = QRegExpValidator(rx)

        self.box.setValidator(self.val)
        self.box.setAlignment(Qt.AlignRight)
        self.suf = QLabel(text=suffix)
        self.box.setText(str(float(value)))

        layout = QHBoxLayout(self)
        # layout.addWidget(self.label)
        layout.addWidget(self.box)
        layout.addWidget(self.suf)

    def value(self):
        return float(self.box.text())

    def set(self, value):
        self.box.setText(str(float(value)))


class BoxLabel(QWidget):
    def __init__(self, label, mini, maxi, value, suffix, *args, **kwargs):
        super(BoxLabel, self).__init__(*args, **kwargs)
        self.box = QLineEdit()
        self.label = QLabel(text=label)

        rx = QRegExp()
        rx.setPattern('\d*')
        self.val = QRegExpValidator(rx)

        self.box.setValidator(self.val)
        self.box.setAlignment(Qt.AlignRight)
        self.suf = QLabel(text=suffix)
        self.box.setText(str(int(value)))

        layout = QHBoxLayout(self)
        layout.addWidget(self.box)
        layout.addWidget(self.suf)

    def value(self):
        return int(self.box.text())

    def set(self, value):
        self.box.setText(str(int(value)))


class StringBoxLabel(QWidget):
    def __init__(self, label, *args, **kwargs):
        super(StringBoxLabel, self).__init__(*args, **kwargs)
        self.box = QLineEdit()
        self.label = QLabel(text=label)

        layout = QHBoxLayout(self)
        layout.addWidget(self.label)
        layout.addWidget(self.box)


class DoubleBoxRange(QWidget):
    def __init__(self, label, *args, **kwargs):
        super(DoubleBoxRange, self).__init__(*args, **kwargs)
        self.label = QLabel(text=label)
        self.lower = QDoubleSpinBox()
        self.upper = QDoubleSpinBox()
        self.dash = QLabel('-')

        self.lower.setMinimum(1.)
        self.upper.setMaximum(20.)
        self.upper.valueChanged.connect(self.lower.setMaximum)
        self.lower.valueChanged.connect(self.upper.setMinimum)
        self.lower.setValue(4.)
        self.upper.setValue(18.5)
        self.lower.setSingleStep(0.5)
        self.upper.setSingleStep(0.5)
        self.lower.setSuffix('μm')
        self.upper.setSuffix('μm')

        layout = QHBoxLayout(self)
        # layout.addWidget(self.label)
        layout.addWidget(self.lower)
        layout.addWidget(self.dash)
        layout.addWidget(self.upper)


class FileBrowser(QWidget):
    def __init__(self, label, *args, **kwargs):
        super(FileBrowser, self).__init__(*args, **kwargs)
        self.label = QLabel(label)
        self.filepath = QLineEdit()
        button = QPushButton('Browse...')

        button.clicked.connect(self.open_browse)

        layout = QHBoxLayout(self)
        layout.addWidget(self.label)
        layout.addWidget(self.filepath)
        layout.addWidget(button)

    def open_browse(self):
        data_path, _ = QFileDialog.getOpenFileName(self, 'Open File')
        self.filepath.setText(data_path)


class FileSaver(QWidget):
    def __init__(self, *args, **kwargs):
        super(FileSaver, self).__init__(*args, **kwargs)
        self.filepath = QLineEdit()
        button = QPushButton('Browse...')
        button.clicked.connect(self.open_browse)

        layout = QHBoxLayout(self)
        layout.addWidget(self.filepath)
        layout.addWidget(button)

    def open_browse(self):
        dialog = QFileDialog()
        dialog.setDefaultSuffix('txt')
        dialog.setAcceptMode(1)
        dialog.exec_()
        data_path = dialog.selectedFiles()
        self.filepath.setText(data_path[0])


class PlotCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=10, height=10, dpi=200):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(PlotCanvas, self).__init__(fig)
