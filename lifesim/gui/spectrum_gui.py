import math

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QApplication, QDialog, QGroupBox, QSlider, QGridLayout, QLabel,
                             QVBoxLayout, QLineEdit, QWidget, QSpinBox, QHBoxLayout,
                             QDoubleSpinBox, QFileDialog, QPushButton, QTabWidget, QCheckBox,
                             QProgressBar)
from PyQt5.QtGui import QPixmap
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

import lifesim as life
from lifesim.modules.util import import_spectrum


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

        # TODO: Remove
        self.filepath.setText('/home/felix/Documents/MA/life_sim_share/input_data/planet_spectra/Earth_Clear_R1000_10pc_Björn_Konrad.txt')

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

        # LIFEsim simulator setup
        self.bus = life.Bus()
        self.options = life.Options()
        self.options.set_scenario('baseline')

        inst = life.Instrument(name='LIFE',
                               options=self.options)
        self.bus.add_module(module=inst)

        tm = life.TransmissionMap(name='tm')
        self.bus.add_module(module=tm)

        self.bus.connect(module_names=('LIFE', 'tm'))

        sl = life.PhotonNoiseStar(name='sl')
        self.bus.add_module(module=sl)
        lz = life.PhotonNoiseLocalzodi(name='lz')
        self.bus.add_module(module=lz)
        ez = life.PhotonNoiseExozodi(name='ez')
        self.bus.add_module(module=ez)

        # ----- GUI Setup ------

        # create all sub-widgets (sorted after functionality)
        self.create_instrument()
        self.create_star()
        self.create_planet()
        self.preview_spectrum()
        self.create_simulation()
        self.create_scenario()
        self.create_logo()
        self.view_spectrum()

        # structure all setting into one widget
        planet_star = QWidget()

        ps_layout = QHBoxLayout(planet_star)
        ps_layout.addWidget(self.star)
        ps_layout.addWidget(self.planet)

        settings = QWidget()

        s_layout = QVBoxLayout(settings)
        s_layout.addWidget(self.instrument, alignment=Qt.AlignTop)
        s_layout.addWidget(planet_star, alignment=Qt.AlignTop)

        # structure the results in one widget
        results = QWidget()

        r_layout = QVBoxLayout(results)
        r_layout.addWidget(self.s_result)

        # create the tab widget
        tabs = QTabWidget()
        tabs.addTab(settings, 'Settings')
        tabs.addTab(results, 'Results')

        #structure the sidebar in one widget
        sidebar = QWidget()
        sb_layout = QVBoxLayout(sidebar)
        sb_layout.addWidget(self.simulation, alignment=Qt.AlignTop)
        sb_layout.addWidget(self.scenario, alignment=Qt.AlignTop)
        sb_layout.addWidget(self.logo, alignment=Qt.AlignBottom)

        # layout the main body of the frame
        body_layout = QGridLayout()
        body_layout.addWidget(tabs, 0, 1)
        body_layout.addWidget(sidebar, 0, 0)
        self.setLayout(body_layout)

        self.setWindowTitle('LIFEsim lite: Spectrum Simulator')

        a=1

    def create_instrument(self):
        self.instrument = QGroupBox('Instrument')

        self.diameter = DoubleBoxLabel(label='Aperture Diameter',
                                       maxi=10.,
                                       mini=0.25,
                                       step=0.5,
                                       value=2,
                                       suffix='m')

        self.wl_range = DoubleBoxRange(label='Wavelength Range')

        self.spec_res = BoxLabel(label='Spectral Resolution',
                                 mini=1,
                                 maxi=100,
                                 value=20,
                                 suffix='')

        line = QWidget()
        l_layout = QHBoxLayout(line)
        l_layout.addWidget(self.diameter, alignment=Qt.AlignLeft)
        l_layout.addWidget(self.spec_res, alignment=Qt.AlignLeft)

        layout = QVBoxLayout(self.instrument)
        layout.addWidget(line)
        layout.addWidget(self.wl_range, alignment=Qt.AlignHCenter)

    def create_star(self):
        self.star = QGroupBox('Star')
        temp = QWidget()

        self.temp_s = DoubleBoxLabel(label='Temperature',
                                     mini=0.,
                                     maxi=100000.,
                                     step=100.,
                                     value=5778.,
                                     suffix='K')
        self.temp_s.box.setDecimals(0)
        self.temp_s.label.setToolTip('Test Tip')

        self.radius_s = DoubleBoxLabel(label='Radius',
                                       mini=0.,
                                       maxi=2158.,
                                       step=0.5,
                                       value=1.,
                                       suffix='R☉')

        self.distance_s = DoubleBoxLabel(label='Distance',
                                         mini=0.,
                                         maxi=50.,
                                         step=1.,
                                         value=10.,
                                         suffix='pc')
        self.distance_s.box.setDecimals(0)

        self.lat = DoubleBoxLabel(label='Galactic Latitude',
                                  mini=0.,
                                  maxi=2*math.pi,
                                  step=0.1,
                                  value=0.79,
                                  suffix='rad')

        self.z = DoubleBoxLabel(label='Exozodis',
                                mini=0.,
                                maxi=20.,
                                step=0.5,
                                value=1.,
                                suffix='z')

        layout = QVBoxLayout(temp)
        layout.addWidget(self.temp_s)
        layout.addWidget(self.radius_s)
        layout.addWidget(self.distance_s)
        layout.addWidget(self.lat)
        layout.addWidget(self.z)

        s_layout = QVBoxLayout(self.star)
        s_layout.addWidget(temp, alignment=Qt.AlignTop)

    def create_planet(self):
        self.planet = QGroupBox('Planet')

        self.angsep = DoubleBoxLabel(label='Angular Separation',
                                mini=0.,
                                maxi=1.,
                                value=0.1,
                                step=0.1,
                                suffix='arcsec')

        self.radius_p = DoubleBoxLabel(label='Radius',
                                  mini=0.,
                                  maxi=10.,
                                  value=1.,
                                  step=0.5,
                                  suffix='R⊕')

        self.radius_spec = DoubleBoxLabel(label='Reference Radius',
                                       mini=0.,
                                       maxi=10.,
                                       value=1.,
                                       step=0.5,
                                       suffix='R⊕')

        self.d_spec = DoubleBoxLabel(label='Reference Distance',
                                         mini=0.,
                                         maxi=50.,
                                         step=1.,
                                         value=10.,
                                         suffix='pc')
        self.d_spec.box.setDecimals(0)

        self.browse = FileBrowser(label='Spectrum')

        layout = QVBoxLayout()
        layout.addWidget(self.angsep)
        layout.addWidget(self.radius_p)
        layout.addWidget(self.radius_spec)
        layout.addWidget(self.d_spec)
        layout.addWidget(self.browse)
        self.planet.setLayout(layout)

    def preview_spectrum(self):
        self.s_preview = QGroupBox('Spectrum Preview')
        self.p_plot = PlotCanvas(self, width=5, height=4, dpi=100)
        # self.p_plot.axes.plot([0, 1, 2, 3, 4], [10, 1, 20, 3, 40])

        button = QPushButton('Preview Spectrum')
        button.clicked.connect(self.show_preview)

        layout = QGridLayout()
        layout.addWidget(button, 0, 0)
        layout.addWidget(self.p_plot, 1, 0)
        self.s_preview.setLayout(layout)

    def view_spectrum(self):
        self.s_result = QWidget()

        self.r_plot = PlotCanvas(self, width=5, height=4, dpi=100)
        # self.p_plot.axes.plot([0, 1, 2, 3, 4], [10, 1, 20, 3, 40])

        layout = QVBoxLayout()
        layout.addWidget(self.r_plot)
        self.s_result.setLayout(layout)

    def create_simulation(self):
        self.simulation = QGroupBox('Simulation')

        time_l = QLabel('Integration Time')

        self.time_b = QDoubleSpinBox()
        self.time_b.setMinimum(0.)
        self.time_b.setMaximum(17520.)
        self.time_b.setSuffix('h')
        self.time_b.setSingleStep(1.)
        self.time_b.setValue(10.)

        t = QWidget()
        t_layout = QVBoxLayout(t)
        t_layout.addWidget(time_l)
        t_layout.addWidget(self.time_b)

        noise_l = QLabel('Noise Sources')

        self.box_sl = QCheckBox('Stellar Leakage')
        self.box_lz = QCheckBox('Localzodi')
        self.box_ez = QCheckBox('Exozodi')
        self.box_sl.setChecked(True)
        self.box_lz.setChecked(True)
        self.box_ez.setChecked(True)

        noise = QWidget()
        n_layout = QVBoxLayout(noise)
        n_layout.addWidget(noise_l)
        n_layout.addWidget(self.box_sl)
        n_layout.addWidget(self.box_lz)
        n_layout.addWidget(self.box_ez)

        simulate = QPushButton('Run Simulation')

        simulate.clicked.connect(self.run_simulation)

        layout = QVBoxLayout(self.simulation)
        layout.addWidget(t, alignment=Qt.AlignBottom)
        layout.addWidget(noise, alignment=Qt.AlignBottom)
        layout.addWidget(simulate, alignment=Qt.AlignBottom)

    def create_scenario(self):
        self.scenario = QGroupBox('Set Scenario')

        optimisic = QPushButton('Optimistic')
        optimisic.clicked.connect(self.set_scenario_opt)
        baseline = QPushButton('Baseline')
        baseline.clicked.connect(self.set_scenario_bas)
        pessimistic = QPushButton('Pessimistic')
        pessimistic.clicked.connect(self.set_scenario_pes)

        layout = QVBoxLayout(self.scenario)
        layout.addWidget(optimisic, alignment=Qt.AlignBottom)
        layout.addWidget(baseline, alignment=Qt.AlignBottom)
        layout.addWidget(pessimistic, alignment=Qt.AlignBottom)

    def create_progress(self):
        self.pb = QGroupBox('Progress')

        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.progress.setValue(50.)
        self.progress.setOrientation(Qt.Vertical)

        layout = QVBoxLayout(self.pb)
        layout.addWidget(self.progress, alignment=Qt.AlignHCenter)

    def create_logo(self):
        self.logo = QWidget()

        image = QLabel()
        pixmap = QPixmap('logo.png').scaledToWidth(200)
        image.setPixmap(pixmap)

        website = QLabel('life-space-mission.com')

        layout = QVBoxLayout(self.logo)
        layout.addWidget(image)
        layout.addWidget(website)


    def update_options(self):
        self.options.set_manual(diameter=self.diameter.box.value(),
                                wl_min=self.wl_range.lower.value(),
                                wl_max=self.wl_range.upper.value(),
                                spec_res=self.spec_res.box.value())
        self.bus.modules['LIFE'].apply_options(self.options)

    def show_preview(self):
        self.update_options()
        p_spec = import_spectrum(pathtofile=self.browse.filepath.text(),
                                 wl_bin_edges=self.bus.modules['LIFE'].wl_bin_edges,
                                 radius_p=self.radius_p.box.value(),
                                 distance_s=self.distance_s.box.value(),
                                 radius_spec=self.radius_spec.box.value(),
                                 distance_spec=self.d_spec.box.value(),
                                 clean=True)
        self.p_plot.axes.cla()
        self.p_plot.axes.plot(p_spec[0], p_spec[1])
        self.p_plot.axes.set_xlabel('Wavelength [m]')
        self.p_plot.axes.set_ylabel('Photon Count [s$^{-1}$]')
        self.p_plot.draw()

    def show_spectrum(self,
                      spectrum):
        self.r_plot.axes.cla()
        self.r_plot.axes.plot(spectrum[0], spectrum[1])
        self.r_plot.axes.set_xlabel('Wavelength [m]')
        self.r_plot.axes.set_ylabel('Signal to Noise Ratio')
        self.r_plot.draw()

    def run_simulation(self):
        self.update_options()

        self.bus.disconnect(module_names=('LIFE', 'sl'))
        self.bus.disconnect(module_names=('LIFE', 'lz'))
        self.bus.disconnect(module_names=('LIFE', 'ez'))
        if self.box_sl.isChecked():
            self.bus.connect(module_names=('LIFE', 'sl'))
        if self.box_lz.isChecked():
            self.bus.connect(module_names=('LIFE', 'lz'))
        if self.box_ez.isChecked():
            self.bus.connect(module_names=('LIFE', 'ez'))

        r_spec = self.bus.modules['LIFE'].get_spectrum(pathtofile=self.browse.filepath.text(),
                                                       temp_s=self.temp_s.box.value(),
                                                       radius_s=self.radius_s.box.value(),
                                                       distance_s=self.distance_s.box.value(),
                                                       radius_spec=self.radius_spec.box.value(),
                                                       distance_spec=self.d_spec.box.value(),
                                                       lat_s=self.lat.box.value(),
                                                       z=self.z.box.value(),
                                                       angsep=self.angsep.box.value(),
                                                       radius_p=self.radius_p.box.value(),
                                                       integration_time=self.time_b.value()*60*60)

        self.show_spectrum(r_spec)

    def set_values(self):
        self.diameter.box.setValue(self.options.array['diameter'])
        self.spec_res.box.setValue(self.options.array['spec_res'])
        self.wl_range.lower.setValue(self.options.array['wl_min'])
        self.wl_range.upper.setValue(self.options.array['wl_max'])

    def set_scenario_opt(self):
        self.options.set_scenario(case='optimistic')
        self.set_values()

    def set_scenario_bas(self):
        self.options.set_scenario(case='baseline')
        self.set_values()

    def set_scenario_pes(self):
        self.options.set_scenario(case='pessimistic')
        self.set_values()


if __name__ == '__main__':
    app = QApplication([])
    gallery = Frame()
    gallery.show()
    app.exec_()

    # <div style="color:darkred;"> <\div>
