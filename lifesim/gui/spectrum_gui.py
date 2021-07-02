import math
import warnings
import urllib.request
from urllib.error import HTTPError, URLError
import os

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QApplication, QDialog, QGroupBox, QGridLayout, QLabel,
                             QVBoxLayout, QWidget, QHBoxLayout, QProgressBar,
                             QDoubleSpinBox, QPushButton, QTabWidget, QCheckBox,
                             QProgressBar, QComboBox, QSizePolicy)
from PyQt5.QtGui import QPixmap, QGuiApplication
import numpy as np
from astropy import units as u
from astropy.visualization import quantity_support

import lifesim as ls
from lifesim.util.radiation import black_body
from lifesim.gui.custom_widgets import (DoubleBoxLabel, BoxLabel, StringBoxLabel, DoubleBoxRange,
                                        FileBrowser, FileSaver, PlotCanvas, RadioButtonWidget)

# change the directory to the path of the spectrum_gui.py file to ensure that the logo is imported
# correctly
os.chdir(os.path.dirname(__file__))

quantity_support()

try:
    urllib.request.urlretrieve("https://github.com/fdannert/LIFEsim/raw/master/"
                               "lifesim/gui/logo_blue.png", "logo_blue.png")
except (HTTPError, URLError) as e:
    pass


class Frame(QDialog):
    def __init__(self, parent=None):
        super(Frame, self).__init__(parent)

        self.bus = ls.Bus()
        self.bus.data.options.set_scenario('baseline')

        self.spec_import = ls.SpectrumImporter()

        inst = ls.Instrument(name='life')
        self.bus.add_module(module=inst)

        transm = ls.TransmissionMap(name='transm')
        self.bus.add_module(module=transm)

        sl = ls.PhotonNoiseStar(name='sl')
        self.bus.add_module(module=sl)

        lz = ls.PhotonNoiseLocalzodi(name='lz')
        self.bus.add_module(module=lz)

        ez = ls.PhotonNoiseExozodi(name='ez')
        self.bus.add_module(module=ez)

        self.bus.connect(('life', 'transm'))
        self.bus.connect(('sl', 'transm'))

        # LIFEsim simulator setup

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

        self.change_visibility()

    def create_instrument(self):
        self.instrument = QGroupBox('Instrument')

        self.diameter = DoubleBoxLabel(label='Aperture Diameter',
                                       maxi=10.,
                                       mini=0.25,
                                       step=0.5,
                                       value=2,
                                       suffix='m')
        self.diameter.label.setToolTip('Diameter of a single aperture of one of the collector spacecraft. The overall'
                                       'collecting area is spanned by four collector spacecraft')

        self.wl_range = DoubleBoxRange(label='Wavelength Range')
        self.wl_range.label.setToolTip('Lower and upper wavelength bound of the spectrograph')

        self.spec_res = BoxLabel(label='Spectral Resolution',
                                 mini=1,
                                 maxi=100,
                                 value=20,
                                 suffix='')
        self.spec_res.label.setToolTip('Spectral resolution of the spectrograph')

        self.baseline_mode = RadioButtonWidget()
        self.baseline_mode.label.setToolTip('The distance between the collector spacecraft (i.e. the baseline) <br> '
                                            'has a major influence on the detectability of exoplanets. It has to be '
                                            'adjusted '
                                            'to match the separation between the target exoplanet and its host star. '
                                            'If the separation of the exoplanet is assumed to be unkown prior to the '
                                            'observation, the baseline can be set to optimal for planets in the HZ. '
                                            'If the separation of the target planet is known, the baseline can be set '
                                            'optimally for the planet position')

        self.baseline_used = QLabel(text='Baseline Used: ')

        layout = QGridLayout(self.instrument)
        layout.addWidget(self.diameter.label, 0, 0)
        layout.addWidget(self.diameter, 0, 1)
        layout.addWidget(self.spec_res. label, 1, 0)
        layout.addWidget(self.spec_res, 1, 1)
        layout.addWidget(self.wl_range.label, 2, 0)
        layout.addWidget(self.wl_range, 2, 1)
        layout.addWidget(self.baseline_mode.label, 3, 0)
        layout.addWidget(self.baseline_mode.bl_HZ, 3, 1)
        layout.addWidget(self.baseline_mode.bl_pl, 4, 1)
        layout.addWidget(self.baseline_used, 5, 0)

    def create_star(self):
        self.star = QGroupBox('Star')
        temp = QWidget()

        self.temp_s = DoubleBoxLabel(label='Temperature',
                                     mini=0.,
                                     maxi=100000.,
                                     step=100.,
                                     value=5778.,
                                     suffix='K')
        self.temp_s.label.setToolTip('Temperature of the host star')

        self.radius_s = DoubleBoxLabel(label='Radius',
                                       mini=0.,
                                       maxi=2158.,
                                       step=0.5,
                                       value=1.,
                                       suffix='R☉')
        self.radius_s.label.setToolTip('Radius of the host star')

        self.distance_s = DoubleBoxLabel(label='Distance',
                                         mini=0.,
                                         maxi=50.,
                                         step=1.,
                                         value=10.,
                                         suffix='pc')
        self.distance_s.label.setToolTip('Distance between the solar system and the host star of the target')

        self.lat = DoubleBoxLabel(label='Ecliptic Latitude',
                                  mini=0.,
                                  maxi=2*math.pi,
                                  step=0.1,
                                  value=0.79,
                                  suffix='rad')
        self.lat.label.setToolTip('Ecliptic latitude of the target. Since the localzodiacal dust lies in the <br>'
                                  'ecplitic of the solar system, this will influence how much localzodi light will '
                                  'leak into the measurement')

        self.z = DoubleBoxLabel(label='Exozodis',
                                mini=0.,
                                maxi=20.,
                                step=0.5,
                                value=1.,
                                suffix='z')
        self.z.label.setToolTip('Exozodi level in the target system. The zodi-number indicates how much more <br>'
                                'zodiacal dust is present in the target system when compared to the solar system in '
                                'terms of surface density')

        layout = QGridLayout(temp)
        layout.addWidget(self.temp_s.label, 0, 0)
        layout.addWidget(self.temp_s, 0, 1)

        layout.addWidget(self.radius_s.label, 1, 0)
        layout.addWidget(self.radius_s, 1, 1)

        layout.addWidget(self.distance_s.label, 2, 0)
        layout.addWidget(self.distance_s, 2, 1)

        layout.addWidget(self.lat.label, 3, 0)
        layout.addWidget(self.lat, 3, 1)

        layout.addWidget(self.z.label, 4, 0)
        layout.addWidget(self.z, 4, 1)

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
        self.angsep.label.setToolTip('Angular separation between target exoplanet and its host star')


        self.radius_p = DoubleBoxLabel(label='Radius',
                                       mini=0.,
                                       maxi=10.,
                                       value=1.,
                                       step=0.5,
                                       suffix='R⊕')
        self.radius_p.label.setToolTip('Radius of the target exoplanets')

        layout = QGridLayout(self.planet)
        layout.addWidget(self.angsep.label, 0, 0)
        layout.addWidget(self.angsep, 0, 1)
        layout.addWidget(self.radius_p.label, 1, 0)
        layout.addWidget(self.radius_p, 1, 1)

    def preview_spectrum(self):
        self.s_preview = QWidget()
        self.p_plot = PlotCanvas(self, width=5, height=4, dpi=100)

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
        self.error_field.setWordWrap(True)


        self.p_plot.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        layout = QVBoxLayout(self.s_preview)
        layout.addWidget(top, alignment=Qt.AlignTop)
        layout.addWidget(self.p_plot)
        layout.addWidget(self.error_field, alignment=Qt.AlignBottom)

    def view_spectrum(self):
        self.s_result = QWidget()

        self.r_plot = PlotCanvas(self, width=5, height=4, dpi=100)

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

        self.simulate = QPushButton('Run Simulation')

        self.simulate.clicked.connect(self.run_simulation)

        self.progress = QProgressBar()

        layout = QVBoxLayout(self.simulation)
        layout.addWidget(t, alignment=Qt.AlignBottom)
        layout.addWidget(noise, alignment=Qt.AlignBottom)
        layout.addWidget(self.simulate, alignment=Qt.AlignBottom)
        layout.addWidget(self.progress, alignment=Qt.AlignHCenter)

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
        pixmap = QPixmap('logo_blue.png').scaledToWidth(150)
        image.setPixmap(pixmap)

        website = QLabel('life-space-mission.com')

        layout = QVBoxLayout(self.logo)
        layout.addWidget(image)
        layout.addWidget(website)

    def create_input(self):
        self.input = QWidget()

        load = QGroupBox('Import')

        self.spec_kind = QComboBox()
        self.spec_kind.addItems(['absolute', 'additive', 'contrast'])
        self.spec_kind.setToolTip('For absolute, the absolute values of the given spectrum are used. For additive, <br>'
                                  'the values of the given spectrum are added onto a blackbody spectrum based on the '
                                  'temperature specified for the target planet. For constrast, the values of the given '
                                  'spectrum are interpreted as contrast values to a blackbody spectrum of the host '
                                  'star')

        self.browse = FileBrowser(label='Spectrum')
        self.browse.label.setToolTip('Location of the file containing the spectrum. The file should contain two <br>'
                                     'space-separated columns. The first column should contain the position of the '
                                     'wavelength bins. The second should contain the flux-like values corresponding to '
                                     'these wavelength bins')

        self.x_units = StringBoxLabel('x-axis units')
        self.x_units.label.setToolTip(r'Units of the wavelength bin location in the given spectrum. Units can be typed '
                                      r'in as free text, e.g. "micron"')
        self.y_units = StringBoxLabel('y-axis units')
        self.y_units.label.setToolTip(r'Units of the flux-like values in the given spectrum. Units can be typed '
                                      r'in as free text, e.g. "photon micron-1 s-1 m-2"')
        self.temp_p = DoubleBoxLabel(label='Temperature Planet',
                                     mini=0.,
                                     maxi=10000.,
                                     step=1.,
                                     value=0.,
                                     suffix='K')
        self.temp_p.label.setToolTip('Temperature of the target planet')

        temp = QWidget()
        ly = QHBoxLayout(temp)
        ly.addWidget(self.temp_p.label)
        ly.addWidget(self.temp_p)

        layout_l = QVBoxLayout(load)
        layout_l.addWidget(self.spec_kind)
        layout_l.addWidget(self.browse)
        layout_l.addWidget(self.x_units)
        layout_l.addWidget(self.y_units)
        layout_l.addWidget(temp)

        self.temp_p.hide()
        self.spec_kind.currentTextChanged.connect(self.change_visibility)

        param = QGroupBox('Spectrum Parameters')

        self.distance_spec = DoubleBoxLabel(label='Distance',
                                            mini=0.,
                                            maxi=50.,
                                            step=1.,
                                            value=0.,
                                            suffix='pc')
        self.distance_spec.label.setToolTip('If a distance to the target system was used during the creation of <br>'
                                            ' the spectrum, specify it here. This will allow LIFEsim to simulate the '
                                            'planet at different distances to the solar system specified in the '
                                            'settings tab')

        self.radius_spec = DoubleBoxLabel(label='Radius Planet',
                                          mini=0.,
                                          maxi=10.,
                                          value=0.,
                                          step=0.5,
                                          suffix='R⊕')
        self.radius_spec.label.setToolTip('If a radius of the target planet was used during the creation of <br>'
                                          ' the spectrum, specify it here. This will allow LIFEsim to simulate the '
                                          'planet in idfferent sizes specified in the settings tab')

        self.time_spec = DoubleBoxLabel(label='Integration Time',
                                        mini=0.,
                                        maxi=60.*60.*24.*365.*10.,
                                        step=60.,
                                        value=0.,
                                        suffix='s')
        self.time_spec.label.setToolTip('If an integration time was used during the creation of the spectrum, specify '
                                        'it here')

        layout_p = QGridLayout(param)
        layout_p.addWidget(self.distance_spec.label, 0, 0)
        layout_p.addWidget(self.distance_spec, 0, 1)
        layout_p.addWidget(self.radius_spec.label, 1, 0)
        layout_p.addWidget(self.radius_spec, 1, 1)
        layout_p.addWidget(self.time_spec.label, 2, 0)
        layout_p.addWidget(self.time_spec, 2, 1)

        layout = QVBoxLayout(self.input)
        layout.addWidget(load)
        layout.addWidget(param)

    def update_options(self):
        self.bus.data.options.set_manual(diameter=self.diameter.value(),
                                         wl_min=self.wl_range.lower.value(),
                                         wl_max=self.wl_range.upper.value(),
                                         spec_res=self.spec_res.value())
        self.bus.modules['life'].apply_options()

    def show_preview(self):
        self.error_field.setText('')
        self.update_options()
        try:
            if ((self.browse.filepath.text() != '')
                    and (not self.spec_kind.currentText() == 'contrast')):
                self.spec_import.do_import(pathtotext=self.browse.filepath.text(),
                                           x_string=self.x_units.box.text(),
                                           y_string=self.y_units.box.text(),
                                           distance_s_spectrum=self.distance_spec.value(),
                                           distance_s_target=self.distance_s.value(),
                                           radius_p_spectrum=self.radius_spec.value(),
                                           radius_p_target=self.radius_p.value(),
                                           integration_time=self.time_spec.value())

                self.flux_planet_spectrum = [self.spec_import.x_data, self.spec_import.y_data]
                if self.spec_kind.currentText() == 'additive':
                    widths = (self.spec_import.x_data.value[1:]
                              - self.spec_import.x_data.value[:-1])
                    widths = np.append(widths, widths[-1])
                    bins = self.spec_import.x_data.value + widths/2
                    fgamma = black_body(mode='planet',
                                        bins=bins,
                                        width=widths,
                                        temp=self.temp_p.value(),
                                        radius=self.radius_p.value(),
                                        distance=self.distance_s.value()) \
                             / widths \
                             * u.photon/u.second/(u.meter**3)
                    self.flux_planet_spectrum[1] += fgamma
                if self.prev_kind.currentIndex() == 0:
                    p_spec = [self.spec_import.x_raw, self.spec_import.y_raw]
                else:
                    p_spec = [self.flux_planet_spectrum[0].to(u.micron),
                              self.flux_planet_spectrum[1]]

            elif self.spec_kind.currentText() == 'additive':
                wl_edge = 1.
                wl_bins = []
                wl_bin_widths = []
                wl_bin_edges = [wl_edge]
                R = 200
                wl_max = 30

                while wl_edge < wl_max:

                    # set the wavelength bin width according to the spectral resolution
                    wl_bin_width = wl_edge / R / \
                                   (1 - 1 / R / 2)

                    # make the last bin shorter when it hits the wavelength limit
                    if wl_edge + wl_bin_width > wl_max:
                        wl_bin_width = wl_max - wl_edge

                    # calculate the center and edges of the bins
                    wl_center = wl_edge + wl_bin_width / 2
                    wl_edge += wl_bin_width

                    wl_bins.append(wl_center)
                    wl_bin_widths.append(wl_bin_width)
                    wl_bin_edges.append(wl_edge)

                # convert everything to [m]
                wl_bins = np.array(wl_bins) * 1e-6  # in m
                wl_bin_widths = np.array(wl_bin_widths) * 1e-6  # in m

                fgamma = black_body(mode='planet',
                                    bins=wl_bins,
                                    width=wl_bin_widths,
                                    temp=self.temp_p.value(),
                                    radius=self.radius_p.value(),
                                    distance=self.distance_s.value()) \
                         / wl_bin_widths \
                         * u.photon / u.second / (u.meter ** 3)

                self.flux_planet_spectrum = [wl_bins * u.meter, fgamma]
                p_spec = [self.flux_planet_spectrum[0].to(u.micron),
                          self.flux_planet_spectrum[1]]

            elif self.spec_kind.currentText() == 'contrast':
                self.spec_import.do_import(pathtotext=self.browse.filepath.text(),
                                           x_string=self.x_units.box.text(),
                                           y_string='photon s-1 m-3',
                                           distance_s_spectrum=None,
                                           distance_s_target=0,
                                           radius_p_spectrum=None,
                                           radius_p_target=0,
                                           integration_time=0)

                self.flux_planet_spectrum = [self.spec_import.x_data, self.spec_import.y_data]

                widths = (self.spec_import.x_data.value[1:]
                          - self.spec_import.x_data.value[:-1])
                widths = np.append(widths, widths[-1])
                bins = self.spec_import.x_data.value + widths / 2
                fgamma = black_body(mode='star',
                                    bins=bins,
                                    width=widths,
                                    temp=self.temp_s.value(),
                                    radius=self.radius_s.value(),
                                    distance=self.distance_s.value()) \
                         / widths \
                         * u.photon / u.second / (u.meter ** 3)
                self.flux_planet_spectrum[1] *= fgamma
                if self.prev_kind.currentIndex() == 0:
                    p_spec = [self.spec_import.x_raw, self.spec_import.y_raw]
                else:
                    p_spec = [self.flux_planet_spectrum[0].to(u.micron),
                              self.flux_planet_spectrum[1]]

            else:
                raise ValueError('Given file cannot be imported')


            self.p_plot.axes.cla()
            self.p_plot.axes.plot(p_spec[0], p_spec[1], color="darkblue", linestyle="-")
            if self.prev_kind.currentIndex() == 1:
                self.p_plot.axes.set_xlim(self.bus.data.options.array['wl_min'],
                                          self.bus.data.options.array['wl_max'])
            self.p_plot.axes.set_ylim(
                (5e-1 * np.min(
                    p_spec[1]
                    [int(np.argwhere(p_spec[0] <
                                     (self.bus.data.options.array['wl_min']*u.micron))[-1]):
                     int(np.argwhere(p_spec[0] >
                                     (self.bus.data.options.array['wl_max']*u.micron))[0])]).value),
                (0.5e1 * np.max(p_spec[1]
                                [int(np.argwhere(p_spec[0] <
                                                 (self.bus.data.options.array['wl_min']
                                                  *u.micron))[-1])
                                 :int(np.argwhere(p_spec[0] >
                                                  (self.bus.data.options.array['wl_max']
                                                   *u.micron))[0])])).value)

            self.p_plot.axes.set_yscale('log')
            self.p_plot.axes.grid()
            self.p_plot.axes.ticklabel_format(axis='x', style='sci')
            self.p_plot.draw()
        except Exception as e:
            self.error_field.setText('Import Error: ' + str(e))

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

        self.r_plot.axes.set_xlim(self.bus.data.options.array['wl_min']-0.5,
                                  self.bus.data.options.array['wl_max']+0.5)
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
        self.progress.setValue(10)
        QGuiApplication.processEvents()

        self.show_preview()
        self.update_options()

        self.progress.setValue(20)
        QGuiApplication.processEvents()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.bus.disconnect(module_names=('life', 'sl'))
            self.bus.disconnect(module_names=('life', 'lz'))
            self.bus.disconnect(module_names=('life', 'ez'))
        if self.box_sl.isChecked():
            self.bus.connect(module_names=('life', 'sl'))
        if self.box_lz.isChecked():
            self.bus.connect(module_names=('life', 'lz'))
        if self.box_ez.isChecked():
            self.bus.connect(module_names=('life', 'ez'))

        self.progress.setValue(30)
        QGuiApplication.processEvents()

        self.r_spec, self.flux_p, self.flux_n = self.bus.modules['life'].get_spectrum(
            temp_s=self.temp_s.value(),
            radius_s=self.radius_s.value(),
            distance_s=self.distance_s.value(),
            flux_planet_spectrum=self.flux_planet_spectrum,
            lat_s=self.lat.value(),
            z=self.z.value(),
            angsep=self.angsep.value(),
            integration_time=self.time_b.value()*60*60,
            baseline_to_planet=self.baseline_mode.bl_pl.isChecked(),
            pbar=self.progress)

        self.progress.setValue(90)
        QGuiApplication.processEvents()

        self.show_spectrum(self.r_spec,
                           self.flux_p,
                           self.flux_n)

        self.baseline_used.setText('Nulling Baseline Used: '
                                   + str(np.round(self.bus.data.inst['bl'],
                                                  decimals=1))
                                   + ' m')

        self.progress.setValue(100)
        QGuiApplication.processEvents()

    def set_values(self):
        self.diameter.set(self.bus.data.options.array['diameter'])
        self.spec_res.set(self.bus.data.options.array['spec_res'])
        self.wl_range.lower.setValue(self.bus.data.options.array['wl_min'])
        self.wl_range.upper.setValue(self.bus.data.options.array['wl_max'])

    def set_scenario_opt(self):
        self.bus.data.options.set_scenario(case='optimistic')
        self.set_values()

    def set_scenario_bas(self):
        self.bus.data.options.set_scenario(case='baseline')
        self.set_values()

    def set_scenario_pes(self):
        self.bus.data.options.set_scenario(case='pessimistic')
        self.set_values()

    def save_spectrum(self):
        if self.r_spec is not None:
            np.savetxt(fname=self.save.filepath.text(),
                       X=np.array([self.r_spec[0], self.r_spec[1], self.flux_p]).T,
                       header='Wavelength [m]   SNR per bin for 1h  Input flux')

    def change_visibility(self):
        if self.spec_kind.currentText() == 'additive':
            self.temp_p.show()
            self.temp_p.label.setText('Temperature Planet')
            self.temp_p.label.show()
            self.x_units.show()
            self.x_units.label.show()
            self.y_units.show()
            self.y_units.label.show()
            self.distance_spec.show()
            self.distance_spec.label.show()
            self.radius_spec.show()
            self.radius_spec.label.show()
            self.time_spec.show()
            self.time_spec.label.show()
        elif self.spec_kind.currentText() == 'absolute':
            self.temp_p.hide()
            self.temp_p.label.hide()
            self.x_units.show()
            self.x_units.label.show()
            self.y_units.show()
            self.y_units.label.show()
            self.distance_spec.show()
            self.distance_spec.label.show()
            self.radius_spec.show()
            self.radius_spec.label.show()
            self.time_spec.show()
            self.time_spec.label.show()
        elif self.spec_kind.currentText() == 'contrast':
            self.temp_p.show()
            self.temp_p.label.setText('Temperature Star')
            self.temp_p.label.show()
            self.x_units.show()
            self.x_units.label.show()
            self.y_units.hide()
            self.y_units.label.hide()
            self.distance_spec.hide()
            self.distance_spec.label.hide()
            self.radius_spec.hide()
            self.radius_spec.label.hide()
            self.time_spec.hide()
            self.time_spec.label.hide()


class Gui(object):
    def __init__(self):
        os.chdir(os.path.dirname(__file__))
        quantity_support()
        app = QApplication([])
        gallery = Frame()
        gallery.show()
        app.exec_()


if __name__ == '__main__':
    app = QApplication([])
    gallery = Frame()
    gallery.show()
    app.exec_()

    # <div style="color:darkred;"> <\div>
