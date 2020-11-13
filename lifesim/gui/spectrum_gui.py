import math

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QApplication, QDialog, QGroupBox, QSlider, QGridLayout, QLabel,
                             QVBoxLayout, QLineEdit, QWidget, QSpinBox, QHBoxLayout,
                             QDoubleSpinBox, QFileDialog, QPushButton, QTabWidget, QCheckBox,
                             QProgressBar)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure


class CustomSlider(QWidget):
    def __init__(self, *args, **kwargs):
        super(CustomSlider, self).__init__(*args, **kwargs)
        self.slider = QSlider(Qt.Horizontal)
        self.numbox = QDoubleSpinBox()
        self.numbox.setRange(self.slider.minimum(), self.slider.maximum())
        self.slider.valueChanged.connect(self.numbox.setValue)
        self.slider.rangeChanged.connect(self.numbox.setRange)
        self.numbox.valueChanged.connect(self.slider.setValue)
        layout = QHBoxLayout(self)
        layout.addWidget(self.numbox)
        layout.addWidget(self.slider)


class DoubleBoxLabel(QWidget):
    def __init__(self, label, mini, maxi, step, value, suffix, *args, **kwargs):
        super(DoubleBoxLabel, self).__init__(*args, **kwargs)
        self.box = QDoubleSpinBox()
        self.label = QLabel(text=label)

        self.box.setMinimum(mini)
        self.box.setMaximum(maxi)
        self.box.setSingleStep(step)
        self.box.setValue(value)
        self.box.setSuffix(suffix)

        layout = QHBoxLayout(self)
        layout.addWidget(self.label)
        layout.addWidget(self.box)


class BoxLabel(QWidget):
    def __init__(self, label, mini, maxi, value, suffix, *args, **kwargs):
        super(BoxLabel, self).__init__(*args, **kwargs)
        self.box = QSpinBox()
        self.label = QLabel(text=label)

        self.box.setMinimum(mini)
        self.box.setMaximum(maxi)
        self.box.setValue(value)
        self.box.setSuffix(suffix)

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
        layout.addWidget(self.label)
        layout.addWidget(self.lower)
        layout.addWidget(self.dash)
        layout.addWidget(self.upper)


class FileBrowser(QWidget):
    def __init__(self, label, *args, **kwargs):
        super(FileBrowser, self).__init__(*args, **kwargs)
        label = QLabel(label)
        self.filepath = QLineEdit()
        button = QPushButton('Browse...')

        button.clicked.connect(self.open_browse)

        layout = QHBoxLayout(self)
        layout.addWidget(label)
        layout.addWidget(self.filepath)
        layout.addWidget(button)

    def open_browse(self):
        data_path, _ = QFileDialog.getOpenFileName(self, 'Open File')
        self.filepath.setText(data_path)


class PlotCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(PlotCanvas, self).__init__(fig)




class Frame(QDialog):
    def __init__(self, parent=None):
        super(Frame, self).__init__(parent)

        self.create_instrument()
        self.create_star()
        self.create_planet()
        self.preview_spectrum()
        self.create_simulation()
        self.create_progress_bar()

        settings = QWidget()

        s_layout = QGridLayout(settings)
        s_layout.addWidget(self.instrument, 0, 0)
        s_layout.addWidget(self.star, 1, 0)
        s_layout.addWidget(self.planet, 0, 1)
        s_layout.addWidget(self.s_preview, 1, 1)

        results = QWidget()

        tabs = QTabWidget()
        tabs.addTab(settings, 'Settings')
        tabs.addTab(results, 'Results')

        sidebar = QWidget()
        sb_layout = QVBoxLayout(sidebar)
        sb_layout.addWidget(self.simulation, alignment=Qt.AlignTop)
        sb_layout.addWidget(self.pb)


        body_layout = QGridLayout()
        body_layout.addWidget(tabs, 0, 1)
        body_layout.addWidget(sidebar, 0, 0)
        self.setLayout(body_layout)

        self.setWindowTitle('LIFEsim lite: Spectrum Simulator')

    def create_instrument(self):
        self.instrument = QGroupBox('Instrument')

        diameter = DoubleBoxLabel(label='Aperture Diameter',
                                  maxi=10.,
                                  mini=0.25,
                                  step=0.5,
                                  value=2.5,
                                  suffix='m')

        wl_range = DoubleBoxRange(label='Wavelength Range')

        spec_res = BoxLabel(label='Spectral Resolution',
                            mini=1,
                            maxi=100,
                            value=20,
                            suffix='')


        layout = QGridLayout()
        layout.addWidget(diameter, 0, 0)
        layout.addWidget(wl_range, 1, 0)
        layout.addWidget(spec_res, 2, 0)
        # layout.addWidget(b_diameter, 0, 1)
        # layout.addWidget(b_diameter, 0, 2)
        self.instrument.setLayout(layout)

    def create_star(self):
        self.star = QGroupBox('Star')
        temp = QWidget()

        temp_s = DoubleBoxLabel(label='Temperature',
                                mini=0.,
                                maxi=100000.,
                                step=100.,
                                value=5778.,
                                suffix='K')
        temp_s.box.setDecimals(0)

        radius_s = DoubleBoxLabel(label='Radius',
                                  mini=0.,
                                  maxi=2158.,
                                  step=0.5,
                                  value=1.,
                                  suffix='R☉')

        distance_s = DoubleBoxLabel(label='Distance',
                                    mini=0.,
                                    maxi=50.,
                                    step=1.,
                                    value=10.,
                                    suffix='pc')
        distance_s.box.setDecimals(0)

        lat = DoubleBoxLabel(label='Galactic Latitude',
                             mini=0.,
                             maxi=2*math.pi,
                             step=0.1,
                             value=0.,
                             suffix='rad')

        z = DoubleBoxLabel(label='Exozodis',
                           mini=0.,
                           maxi=20.,
                           step=0.5,
                           value=1.,
                           suffix='z')

        layout = QVBoxLayout(temp)
        layout.addWidget(temp_s)
        layout.addWidget(radius_s)
        layout.addWidget(distance_s)
        layout.addWidget(lat)
        layout.addWidget(z)

        s_layout = QVBoxLayout(self.star)
        s_layout.addWidget(temp, alignment=Qt.AlignTop)

    def create_planet(self):
        self.planet = QGroupBox('Planet')

        angsep = DoubleBoxLabel(label='Angular Separation',
                                mini=0.,
                                maxi=1.,
                                value=0.1,
                                step=0.1,
                                suffix='arcsec')

        radius_p = DoubleBoxLabel(label='Radius',
                                  mini=0.,
                                  maxi=10.,
                                  value=1.,
                                  step=0.5,
                                  suffix='R⊕')

        browse = FileBrowser(label='Spectrum')

        # browse = QPushButton('Browse...')
        # browse.clicked.connect(self.open_browse)
        #
        # self.filepath1 = QLineEdit()

        layout = QGridLayout()
        layout.addWidget(angsep, 0, 0)
        layout.addWidget(radius_p, 1, 0)
        layout.addWidget(browse, 2, 0)
        self.planet.setLayout(layout)

    def preview_spectrum(self):
        self.s_preview = QGroupBox('Spectrum Preview')
        plot = PlotCanvas(self, width=5, height=4, dpi=100)
        plot.axes.plot([0, 1, 2, 3, 4], [10, 1, 20, 3, 40])

        layout = QGridLayout()
        layout.addWidget(plot, 0, 0)
        self.s_preview.setLayout(layout)

    def create_simulation(self):
        self.simulation = QGroupBox('Simulation')

        time_l = QLabel('Integration Time')

        time_b = QDoubleSpinBox()
        time_b.setMinimum(0.)
        time_b.setMaximum(17520.)
        time_b.setSuffix('h')
        time_b.setSingleStep(1.)
        time_b.setValue(10.)

        t = QWidget()
        t_layout = QVBoxLayout(t)
        t_layout.addWidget(time_l)
        t_layout.addWidget(time_b)

        noise_l = QLabel('Noise Sources')

        sl = QCheckBox('Stellar Leakage')
        lz = QCheckBox('Localzodi')
        ez = QCheckBox('Exozodi')
        sl.setChecked(True)
        lz.setChecked(True)
        ez.setChecked(True)

        noise = QWidget()
        n_layout = QVBoxLayout(noise)
        n_layout.addWidget(noise_l)
        n_layout.addWidget(sl)
        n_layout.addWidget(lz)
        n_layout.addWidget(ez)

        simulate = QPushButton('Run Simulation')

        layout = QVBoxLayout(self.simulation)
        layout.addWidget(t, alignment=Qt.AlignBottom)
        layout.addWidget(noise, alignment=Qt.AlignBottom)
        layout.addWidget(simulate, alignment=Qt.AlignBottom)

    def create_progress_bar(self):
        self.pb = QGroupBox('Progress')

        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.progress.setValue(50.)
        self.progress.setOrientation(Qt.Vertical)

        layout = QVBoxLayout(self.pb)
        layout.addWidget(self.progress, alignment=Qt.AlignHCenter)

if __name__ == '__main__':
    app = QApplication([])
    gallery = Frame()
    gallery.show()
    app.exec_()
