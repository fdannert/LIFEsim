import math
import os

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QApplication, QDialog, QGroupBox, QSlider, QGridLayout, QLabel,
                             QVBoxLayout, QLineEdit, QWidget, QSpinBox, QHBoxLayout,
                             QDoubleSpinBox, QFileDialog, QPushButton, QTabWidget, QCheckBox,
                             QProgressBar, QComboBox, QSizePolicy)
from PyQt5.QtGui import QPixmap
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import numpy as np
from astropy import units as u
from astropy.visualization import quantity_support

import lifesim as life
from lifesim.modules.util import import_spectrum

# change the directory to the path of the spectrum_gui.py file to ensure that the logo is imported
# correctly
import os
os.chdir(os.path.dirname(__file__))

quantity_support()


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


class Frame(QDialog):
    def __init__(self, parent=None):
        super(Frame, self).__init__(parent)

        # LIFEsim simulator setup
        self.bus = life.Bus()
        self.options = life.Options()
        self.options.set_scenario('baseline')

        self.spec_import = life.SpectrumImporter()

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
        self.create_input()

        # structure all setting into one widget
        planet_inst = QWidget()

        ps_layout = QVBoxLayout(planet_inst)
        ps_layout.addWidget(self.instrument)
        ps_layout.addWidget(self.planet)

        settings = QWidget()

        s_layout = QHBoxLayout(settings)
        s_layout.addWidget(self.star)
        s_layout.addWidget(planet_inst)

        # Inputs
        preview = QWidget()

        p_layout = QHBoxLayout(preview)
        p_layout.addWidget(self.input)
        p_layout.addWidget(self.s_preview)

        # create the tab widget
        tabs = QTabWidget()
        tabs.addTab(settings, 'Settings')
        tabs.addTab(preview, 'Preview')
        tabs.addTab(self.s_result, 'Results')

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

        self.r_spec = None

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

        # line = QWidget()
        # l_layout = QHBoxLayout(line)
        # l_layout.addWidget(self.diameter, alignment=Qt.AlignLeft)
        # l_layout.addWidget(self.spec_res, alignment=Qt.AlignLeft)

        layout = QVBoxLayout(self.instrument)
        layout.addWidget(self.diameter)
        layout.addWidget(self.spec_res)
        layout.addWidget(self.wl_range)

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
        s_layout.addWidget(temp)

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

        layout = QVBoxLayout()
        layout.addWidget(self.angsep)
        layout.addWidget(self.radius_p)
        self.planet.setLayout(layout)

    def preview_spectrum(self):
        self.s_preview = QWidget()
        self.p_plot = PlotCanvas(self, width=5, height=4, dpi=100)
        # self.p_plot.axes.plot([0, 1, 2, 3, 4], [10, 1, 20, 3, 40])

        top = QWidget()
        self.prev_kind = QComboBox()
        self.prev_kind.addItems(['given units', 'converted units'])

        button = QPushButton('Preview Spectrum')
        button.clicked.connect(self.show_preview)

        layout_t = QHBoxLayout(top)
        layout_t.addWidget(self.prev_kind)
        layout_t.addWidget(button)

        self.error_field = QLabel()
        self.error_field.setStyleSheet("QLabel { color : red; }")

        self.p_plot.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        # self.p_plot.setMinimumWidth(self.p_plot.height())

        layout = QVBoxLayout(self.s_preview)
        layout.addWidget(top, alignment=Qt.AlignTop)
        layout.addWidget(self.p_plot)
        layout.addWidget(self.error_field, alignment=Qt.AlignBottom)

    def view_spectrum(self):
        self.s_result = QWidget()

        self.r_plot = PlotCanvas(self, width=5, height=4, dpi=100)
        # self.p_plot.axes.plot([0, 1, 2, 3, 4], [10, 1, 20, 3, 40])

        self.save_dialog = QGroupBox('Save Spectrum to File')

        self.save = FileSaver()
        button = QPushButton('Save')
        button.clicked.connect(self.save_spectrum)

        layout = QHBoxLayout(self.save_dialog)
        layout.addWidget(self.save)
        layout.addWidget(button)
        self.r_plot.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        layout = QVBoxLayout()
        layout.addWidget(self.r_plot)
        layout.addWidget(self.save_dialog, alignment=Qt.AlignBottom)
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
        pixmap = QPixmap('logo_blue.png').scaledToWidth(300)
        image.setPixmap(pixmap)

        website = QLabel('life-space-mission.com')

        layout = QVBoxLayout(self.logo)
        layout.addWidget(image)
        layout.addWidget(website)

    # def create_save(self):
    #     self.save_dialog = QGroupBox('Save Spectrum to File')
    #
    #     self.save = FileSaver()
    #     button = QPushButton('Save')
    #     button.clicked.connect(self.save_spectrum)
    #
    #     layout = QHBoxLayout(self.save_dialog)
    #     layout.addWidget(self.save)
    #     layout.addWidget(button)

    def create_input(self):
        self.input = QWidget()

        load = QGroupBox('Import')

        self.spec_kind = QComboBox()
        self.spec_kind.addItems(['absolute', 'additive'])

        self.browse = FileBrowser(label='Spectrum')

        self.x_units = StringBoxLabel('x-axis units')
        self.y_units = StringBoxLabel('y-axis units')
        self.x_units.box.setText('micron')
        self.y_units.box.setText('photon micron-1 h-1 m-2')


        layout_l = QVBoxLayout(load)
        layout_l.addWidget(self.spec_kind)
        layout_l.addWidget(self.browse)
        layout_l.addWidget(self.x_units)
        layout_l.addWidget(self.y_units)

        param = QGroupBox('Spectrum Parameters')

        self.distance_spec = DoubleBoxLabel(label='Distance',
                                       mini=0.,
                                       maxi=50.,
                                       step=1.,
                                       value=0.,
                                       suffix='pc')
        self.distance_spec.box.setDecimals(0)
        self.distance_spec.box.setValue(10.)

        self.radius_spec = DoubleBoxLabel(label='Radius',
                                     mini=0.,
                                     maxi=10.,
                                     value=0.,
                                     step=0.5,
                                     suffix='R⊕')
        self.radius_spec.box.setValue(1.)

        self.time_spec = DoubleBoxLabel(label='Integration Time',
                                   mini=0.,
                                   maxi=60.*60.*24.*365.*10.,
                                   step=60.,
                                   value=0.,
                                   suffix='s')

        layout_p = QVBoxLayout(param)
        layout_p.addWidget(self.distance_spec)
        layout_p.addWidget(self.radius_spec)
        layout_p.addWidget(self.time_spec)

        layout = QVBoxLayout(self.input)
        layout.addWidget(load)
        layout.addWidget(param)

    def update_options(self):
        self.options.set_manual(diameter=self.diameter.box.value(),
                                wl_min=self.wl_range.lower.value(),
                                wl_max=self.wl_range.upper.value(),
                                spec_res=self.spec_res.box.value())
        self.bus.modules['LIFE'].apply_options(self.options)

    def show_preview(self):
        self.error_field.setText('')
        self.update_options()
        try:
            self.spec_import.do_import(pathtotext=self.browse.filepath.text(),
                                       x_string=self.x_units.box.text(),
                                       y_string=self.y_units.box.text(),
                                       distance_s_spectrum=self.distance_spec.box.value(),
                                       distance_s_target=self.distance_s.box.value(),
                                       radius_p_spectrum=self.radius_spec.box.value(),
                                       radius_p_target=self.radius_p.box.value(),
                                       integration_time=self.time_spec.box.value())
            if self.prev_kind.currentIndex() == 0:
                p_spec = [self.spec_import.x_raw, self.spec_import.y_raw]
            else:
                p_spec = [self.spec_import.x_data.to(u.micron), self.spec_import.y_data]

            self.flux_planet_spectrum = [self.spec_import.x_data, self.spec_import.y_data]
            self.p_plot.axes.cla()
            self.p_plot.axes.plot(p_spec[0], p_spec[1], color="darkblue", linestyle="-")
            if self.prev_kind.currentIndex() == 1:
                self.p_plot.axes.set_xlim(self.options.array['wl_min'], self.options.array['wl_max'])
            self.p_plot.axes.set_ylim(
                (5e-1 * np.min(
                    p_spec[1]
                    [int(np.argwhere(p_spec[0] < (self.options.array['wl_min']*u.micron))[-1]):
                     int(np.argwhere(p_spec[0] > (self.options.array['wl_max']*u.micron))[0])]).value),
                (0.5e1 * np.max(p_spec[1]
                             [int(np.argwhere(p_spec[0] < (self.options.array['wl_min']*u.micron))[-1])
                              :int(np.argwhere(p_spec[0] > (self.options.array['wl_max']*u.micron))[0])])).value)
            # self.p_plot.axes.set_xlabel(r"$\lambda$ [$\mu$m]")
            # self.p_plot.axes.set_ylabel(r"Input signal [ph s$^{-1}$m$^{-2}$m$^{-1}$]")
            self.p_plot.axes.set_yscale('log')
            self.p_plot.axes.grid()
            self.p_plot.axes.ticklabel_format(axis='x', style='sci')
            self.p_plot.draw()
        except:
            self.error_field.setText('An error occured during import')

    def show_spectrum(self,
                      spectrum,
                      planet,
                      noise):
        self.r_plot.axes.cla()
        self.r_plot.axes.step(spectrum[0] * 1e6, planet / (self.time_b.value()*60*60),
                              where="mid", color="r", label=f"Planet")
        self.r_plot.axes.step(spectrum[0] * 1e6, noise / (self.time_b.value()*60*60),
                              where="mid", color="black", label=f"Noise sources")
        self.r_plot.axes.set_xlabel(r"$\lambda$ [$\mu$m]")
        self.r_plot.axes.set_ylabel(r"Detected signal per bin [e$^-$ s$^{-1}$ bin$^{-1}$]")

        self.r_plot.axes.set_xlim(self.options.array['wl_min']-0.5,
                                  self.options.array['wl_max']+0.5)
        self.r_plot.axes.set_yscale('log')
        self.r_plot.axes.grid()

        ax2a = self.r_plot.axes.twinx()
        ax2a.set_yscale('log')
        ax2a.step(spectrum[0] * 1e6, spectrum[1],
                  where="mid", color="darkblue", linestyle="--",
                  label=f"SNR per bin \nTotal: {np.sqrt((spectrum[1] ** 2).sum()):.2f}")
        ax2a.set_ylabel('SNR per bin')
        lines, labels = self.r_plot.axes.get_legend_handles_labels()
        lines2, labels2 = ax2a.get_legend_handles_labels()
        ax2a.legend(lines + lines2, labels + labels2, framealpha=1)
        ax2a.grid(False)

        d_min = np.min(np.array((spectrum[1],
                                 planet / (self.time_b.value()*60*60),
                                 noise / (self.time_b.value()*60*60)))) * 1e-1
        d_max = np.max(np.array((spectrum[1],
                                 planet / (self.time_b.value()*60*60),
                                 noise / (self.time_b.value()*60*60)))) * 1e1
        self.r_plot.axes.set_ylim(d_min, d_max)
        ax2a.set_ylim(d_min, d_max)

        self.r_plot.draw()

        ax2a.remove()

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

        self.r_spec, self.flux_p, self.flux_n = self.bus.modules['LIFE'].get_spectrum(
            temp_s=self.temp_s.box.value(),
            radius_s=self.radius_s.box.value(),
            distance_s=self.distance_s.box.value(),
            flux_planet_spectrum=self.flux_planet_spectrum,
            lat_s=self.lat.box.value(),
            z=self.z.box.value(),
            angsep=self.angsep.box.value(),
            integration_time=self.time_b.value()*60*60)

        print(self.r_spec[1][np.argmin(np.abs(self.r_spec[0]-11e-6))])

        self.show_spectrum(self.r_spec,
                           self.flux_p,
                           self.flux_n)

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

    def save_spectrum(self):
        if self.r_spec is not None:
            # file = open(self.save.filepath.text(), 'w')
            # file.write('Wavelength  SNR per bin\n')
            # file.writelines([self.r_spec[0], self.r_spec[1]])
            # file.close()
            np.savetxt(fname=self.save.filepath.text(),
                       X=np.array([self.r_spec[0], self.r_spec[1], self.flux_p]).T,
                       header='Wavelength [m]   SNR per bin for 1h  Input flux')


if __name__ == '__main__':
    app = QApplication([])
    gallery = Frame()
    gallery.show()
    app.exec_()

    # <div style="color:darkred;"> <\div>
