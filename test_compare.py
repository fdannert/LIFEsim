# %%
import time
from numpy import cos, sin, pi
import numpy as np
from numpy.lib.function_base import diff
import pandas as pd
from tqdm import tqdm
from matplotlib.colors import Colormap, Normalize
import matplotlib.pyplot as plt
import lifesim as ls
# ipython magic
%matplotlib

"""
Test and compare the result of timedependent simulation vs original time independent simulation
"""
# %%


def setup_bus():
    bus = ls.Bus()
    bus.data.options.set_scenario('baseline')
    bus.data.import_catalog(
        "./output/suitable.hdf5")
    return bus

# %%


bus = setup_bus()  # new sim

# %%
# ---------- MEASUREMENT SIMULATION ---------- copied from lifesim_demo.py by fdannert


# create modules and add to bus
inst = ls.Instrument(name='inst')
bus.add_module(inst)

transm = ls.TransmissionMap(name='transm')
bus.add_module(transm)

exo = ls.PhotonNoiseExozodi(name='exo')
bus.add_module(exo)
local = ls.PhotonNoiseLocalzodi(name='local')
bus.add_module(local)
star = ls.PhotonNoiseStar(name='star')
bus.add_module(star)

# connect all modules
bus.connect(('inst', 'transm'))
bus.connect(('inst', 'exo'))
bus.connect(('inst', 'local'))
bus.connect(('inst', 'star'))

bus.connect(('star', 'transm'))


# %%


def test_anomalies(n_angles):
    snrs = []
    true_anomalies = np.arange(0, 2 * np.pi, np.pi * 2 / n_angles)
    for anomaly in true_anomalies:
        bus.data.catalog["theta_p"] = anomaly
        inst.get_snr_t()
        snrs.append(bus.data.catalog[["theta_p", "snr_1h", "inc_p"]].copy())
    return snrs
# %%


def plot_transmission_curve(index, inst, rotations=1, rotation_steps=360, label="", time_dependent=True):
    plt.figure()
    tr_chop_one_wv = inst.get_transmission_curve(
        index, time_dependent=time_dependent)[0][0, 0]
    plt.plot(np.linspace(0, rotations * 2, rotations * rotation_steps),
             tr_chop_one_wv, label=label, alpha=0.5)
    # tr_chop_one_wv/np.max(np.abs(tr_chop_one_wv)), label=label, alpha=0.5)
    plt.xlabel(r'$\pi$')
    plt.ylabel("normalised interference")


# %%

def test_one_anomalies(index, instrument, n_angles=10, rotation_period=12, rotations=1, rotation_steps=360):
    thetas = np.arange(0, 2 * np.pi, np.pi * 2 / n_angles)
    inst.set_rotation(rotation_period, rotations, rotation_steps)
    for theta in thetas:
        bus.data.catalog["theta_p"].iat[index] = theta
        plot_transmission_curve(
            instrument, index, rotations, rotation_steps, label=str(theta))
    plt.legend()
    plt.show()
# %%


def plot_all(instrument, rotation_period, rotations, rotation_steps=360):
    instrument.set_rotation(rotation_period, rotations, rotation_steps)
    for i in range(len(bus.data.catalog)):
        # for i in [1, 2, 3]:
        plot_transmission_curve(i, instrument, rotations,
                                rotation_steps, label=str(i))
    plt.legend()
    plt.show()
# %%


def compare_verticalplot(rotation_period=12, rotations=1, rotation_steps=360):
    inst.set_rotation(rotation_period, rotations, rotation_steps)
    plt.figure()
    indeces = [0, 1, 3, 4]
    labels = ["Edge on", "Face On Clockwise ",
              "Face On Anticlockwise", "Face On Long Period"]
    inst.set_rotation(rotation_period, rotations, rotation_steps)
    for i, index in enumerate(indeces):
        tr_chop_one_wv = inst.get_transmission_curve(
            index, time_dependent=False)[0][0, 0]
        plt.subplot(2, len(indeces), i+1)
        if i == 0:
            plt.ylabel("Transmission w/o Time")
        plt.plot(np.linspace(0, rotations * 2 * pi, rotations * rotation_steps),
                 tr_chop_one_wv, color="blue")
        plt.xticks([0, pi/2, pi, 3/2*pi, 2*pi],
                   ["0", r'$\pi/2$', r'$\pi$', r'$3/2\pi$', r'$2\pi$'])
        plt.title(labels[i])
    for i, index in enumerate(indeces):
        tr_chop_one_wv = inst.get_transmission_curve(
            index, time_dependent=True)[0][0, 0]
        plt.subplot(2, len(indeces), i+5)
        if i == 0:
            plt.ylabel("Transmission With Time")
        plt.plot(np.linspace(0, rotations * 2 * pi, rotations * rotation_steps),
                 tr_chop_one_wv, color="orange")
        plt.xticks([0, pi/2, pi, 3/2*pi, 2*pi],
                   ["0", r'$\pi/2$', r'$\pi$', r'$3/2\pi$', r'$2\pi$'])
        plt.xlabel("Radians")
    plt.tight_layout()
    plt.show()


def compare_t_o(index, rotation_period, rotations, rotation_steps=360):
    inst.set_rotation(rotation_period, rotations, rotation_steps)
    plot_transmission_curve(index, inst, rotations, rotation_steps,
                            label="t independent", time_dependent=False)
    plot_transmission_curve(
        index, inst, rotations, rotation_steps, label="t dependent", time_dependent=True)
    plt.legend()
    plt.show()
# %%


def transmission_and_path(index, rotation_period=12, rotations=1, rotation_steps=360):
    image_size = 512
    inst.set_rotation(rotation_period, rotations, rotation_steps)
    inst.apply_options()
    inst.adjust_bl_to_hz(hz_center=float(bus.data.catalog.hz_center.iloc[index]),
                         distance_s=float(bus.data.catalog.distance_s.iloc[index]))
    _, _, _, _, tm_chop = inst.run_socket(s_name='transmission',
                                          method='transmission_map',
                                          map_selection=[
                                              "tm_3", "tm_4", "tm_chop"],
                                          direct_mode=False,
                                          hfov=None,
                                          image_size=image_size)
    angsep = bus.data.catalog.iloc[index]["maxangsep"]
    theta_p = bus.data.catalog.iloc[index]["theta_p"]
    inc = bus.data.catalog.iloc[index]["inc_p"]
    p_orb = bus.data.catalog.iloc[index]["p_orb"] * 24
    thetas = np.linspace(theta_p, theta_p + 2*pi * rotation_period / p_orb * rotations,
                         rotations * rotation_steps, endpoint=False, dtype="float")
    v = inst.run_socket(s_name="transmission",
                        method="projected_vector",
                        theta=thetas,
                        inc=inc,
                        rad=angsep)
    v = np.array(v)*2000
    phi_lin = np.linspace(0, 2 * (rotations) * np.pi,
                          rotations * rotation_steps, endpoint=False)
    xs = v[0] * np.cos(phi_lin) + v[1] * (-np.sin(phi_lin))
    ys = v[0] * np.sin(phi_lin) + v[1] * np.cos(phi_lin)
    mas_pix = bus.data.inst["mas_pix"][0]
    xs = xs/mas_pix + image_size//2
    ys = ys/mas_pix + image_size//2
    plt.figure()
    plt.imshow(tm_chop[0], vmin=-1, vmax=1, cmap="twilight")
    print(mas_pix)
    plt.plot(xs[0],ys[0],"x")
    plt.plot(image_size//2,image_size//2,"*", color = "yellow")
    plt.plot(xs, ys, color = "white")
    # plt.ylim(-1, 1)
    # plt.xlim(-1, 1)
    plt.show()

#%%

def artificial_t_a_p(angsep, p_orb, inc_p):
    index = 0
    rotation_period = 12
    rotations = 1
    rotation_steps = 360
    image_size = 512
    inst.set_rotation(rotation_period, rotations, rotation_steps)
    inst.apply_options()
    inst.adjust_bl_to_hz(hz_center=float(bus.data.catalog.hz_center.iloc[index]),
                         distance_s=float(bus.data.catalog.distance_s.iloc[index]))
    _, _, _, _, tm_chop = inst.run_socket(s_name='transmission',
                                          method='transmission_map',
                                          map_selection=[
                                              "tm_3", "tm_4", "tm_chop"],
                                          direct_mode=False,
                                          hfov=None,
                                          image_size=image_size)
    angsep = angsep
    theta_p = bus.data.catalog.iloc[index]["theta_p"]
    inc = inc_p
    p_orb = p_orb
    thetas = np.linspace(theta_p, theta_p + 2*pi * rotation_period / p_orb * rotations,
                         rotations * rotation_steps, endpoint=False, dtype="float")
    v = inst.run_socket(s_name="transmission",
                        method="projected_vector",
                        theta=thetas,
                        inc=inc,
                        rad=angsep)

    phi_lin = np.linspace(0, 2 * (rotations) * np.pi,
                          rotations * rotation_steps, endpoint=False)
    xs = v[0] * np.cos(phi_lin) + v[1] * (-np.sin(phi_lin))
    ys = v[0] * np.sin(phi_lin) + v[1] * np.cos(phi_lin)
    mas_pix = bus.data.inst["mas_pix"][0]
    xs = xs/mas_pix + image_size//2
    ys = ys/mas_pix + image_size//2
    plt.figure()
    plt.imshow(tm_chop[0], vmin=-1, vmax=1, cmap="twilight")
    print(mas_pix)

    plt.plot(xs, ys, "--", color="white")
    plt.plot(256, 256, color="yellow", marker="*")
    size = 100
    plt.ylim(256-size, 256+size)
    plt.xlim(256-size, 256+size)
    plt.show()
# %%


def periods(p):
    trc = []
    inst.set_rotation(12, 1, 360)
    inst.apply_options()
    for period in p:
        bus.data.catalog["p_orb"] = period
        plot_transmission_curve(3, inst, rotations=1, rotation_steps=360, label=str(
            period), time_dependent=True)
    plot_transmission_curve(3, inst, rotations=1, rotation_steps=360,
                            label="stationary", time_dependent=False)
    plt.legend()
    plt.show()

# %%


def ffts(periods):
    trc = []
    inst.set_rotation(12, 1, 360)
    inst.apply_options()
    for period in periods:
        bus.data.catalog["p_orb"] = period
        tr_chop_one_wv = inst.get_transmission_curve(3)[0][0, 0]
        ft = np.fft.rfft(tr_chop_one_wv)
        plt.plot(np.arange(0, len(ft)), np.abs(ft), label=str(period))
    plt.xlim(-1, 25)
    plt.legend()
    plt.show()
# %%


def simulate(rotation_period=12, ratio=1):
    bus = ls.Bus()
    bus.data.options.set_scenario('baseline')
    bus.data.import_catalog(
        input_path='./output/TestSet.hdf5')
    inst = ls.Instrument(name='inst')
    bus.add_module(inst)
    if ratio > 0 and ratio < 1:
        number = int(len(bus.data.catalog)*ratio)
        if number != 0:
            bus.data.catalog = bus.data.catalog[:number]

    transm = ls.TransmissionMap(name='transm')
    bus.add_module(transm)

    exo = ls.PhotonNoiseExozodi(name='exo')
    bus.add_module(exo)
    local = ls.PhotonNoiseLocalzodi(name='local')
    bus.add_module(local)
    star = ls.PhotonNoiseStar(name='star')
    bus.add_module(star)

    # connect all modules
    bus.connect(('inst', 'transm'))
    bus.connect(('inst', 'exo'))
    bus.connect(('inst', 'local'))
    bus.connect(('inst', 'star'))

    bus.connect(('star', 'transm'))

    inst.set_rotation(rotation_period, 1, 360)

    detected_o = []
    inst.get_snr_t(time_dependent=False)
    data_o = bus.data.catalog.copy()
    m = data_o[data_o["stype"] == 4]
    k = data_o[data_o["stype"] == 3]
    g = data_o[data_o["stype"] == 2]
    f = data_o[data_o["stype"] == 1]
    # TODO inst.get_snr_t() observes the whole rotation period of the instrument therefore its snr_12h not snr_1h
    detected_o.extend([len(m[m["snr_1h"] > 7]), len(k[k["snr_1h"] > 7]), len(
        g[g["snr_1h"] > 7]), len(f[f["snr_1h"] > 7])])

    detected_t = []
    inst.get_snr_t()
    data_t = bus.data.catalog.copy()
    m = data_t[data_t["stype"] == 4]
    k = data_t[data_t["stype"] == 3]
    g = data_t[data_t["stype"] == 2]
    f = data_t[data_t["stype"] == 1]

    detected_t.extend([len(m[m["snr_1h"] > 7]), len(k[k["snr_1h"] > 7]), len(
        g[g["snr_1h"] > 7]), len(f[f["snr_1h"] > 7])])
    return detected_o, detected_t, data_o, data_t


def compare_detected(o, t):
    plt.figure()
    plt.bar([0, 1, 2, 3], o, 0.25, tick_label=[
            "M", "K", "G", "F"], label="Original")

    plt.bar([0.25, 1.25, 2.25, 3.25], t, 0.25, tick_label=[
            "M", "K", "G", "F"], label="Time_dependent")
    plt.ylabel("#Detected")
    plt.xlabel("Startype")
    plt.legend()
    plt.show()
# %%


def trunc(n, o=2):
    return int(10**o*n)/10**o


def inper(n):
    return int(n*100)


def analysis(do, dt):
    # Periods histogram
    if False:

        plt.figure()
        plt.title("Simulated Orbit Period Distribution")
        min_period = min(do["p_orb"])
        max_period = max(do["p_orb"])
        p_bins = np.arange(min_period, max_period, 1)
        _, _, patches = plt.hist(do["p_orb"], p_bins, alpha=0.5,
                                 color="blue", density=True)
        for i in range(len(patches)):
            if i < 100:
                patches[i].set_fc("yellow")
            if i < 50:
                patches[i].set_fc("orange")
            if i < 4:
                patches[i].set_fc("r")
        plt.xlabel("Orbit Period in Days")
        plt.ylabel("Percent of planets")
        plt.axvline(x=4.17, label="1h Instrument: " +
                    str(inper(len(dt[dt.p_orb < 4.17])/len(dt))) + "%", color="red")
        plt.axvline(x=50, label="12h Instrument:" +
                    str(inper(len(dt[dt.p_orb < 50])/len(dt))) + "%", color="orange")
        plt.axvline(x=100, label="24h Instrument:" +
                    str(inper(len(dt[dt.p_orb < 100])/len(dt))) + "%", color="yellow")
        # plt.ylim(0, 5000)
        plt.xlim(-1, 150)
        plt.legend()
        plt.show()
    # Periods histogram of only detected
    if False:
        plt.figure()
        plt.title("Simulated Orbit Period Distribution of only detected")
        detected = dt[dt["detected"]]
        min_period = min(detected["p_orb"])
        max_period = max(detected["p_orb"])
        p_bins = np.arange(min_period, max_period, 1)
        _, _, patches = plt.hist(detected["p_orb"], p_bins, alpha=0.5,
                                 color="blue", density=True)
        for i in range(len(patches)):
            if i < 100:
                patches[i].set_fc("yellow")
            if i < 50:
                patches[i].set_fc("orange")
            if i < 4:
                patches[i].set_fc("r")
        plt.xlabel("Orbit Period in Days")
        plt.ylabel("Percent of planets")
        plt.axvline(x=4.17, label="1h Instrument: " +
                    str(inper(len(detected[detected.p_orb < 4.17])/len(detected))) + "%", color="red")
        plt.axvline(x=50, label="12h Instrument:" +
                    str(inper(len(detected[detected.p_orb < 50])/len(detected))) + "%", color="orange")
        plt.axvline(x=100, label="24h Instrument:" +
                    str(inper(len(detected[detected.p_orb < 100])/len(detected))) + "%", color="yellow")
        # plt.ylim(0, 5000)
        plt.xlim(-1, 150)
        plt.legend()
        plt.show()
    #Mstar
    if True:
        plt.figure(figsize=(9,3))
        plt.title("Orbit period distribution with F-star host")
        detected = dt[dt["detected"]]
        detected = detected[detected["stype"]==1]
        min_period = min(detected["p_orb"])
        max_period = max(detected["p_orb"])
        p_bins = np.arange(min_period, max_period, 1)
        _, _, patches = plt.hist(detected["p_orb"], p_bins, alpha=0.5,
                                 color="blue", density=True)
        for i in range(len(patches)):
            if i < 100:
                patches[i].set_fc("yellow")
            if i < 50:
                patches[i].set_fc("orange")
            if i < 4:
                patches[i].set_fc("r")
        plt.xlabel("Orbit Period in Days")
        plt.ylabel("Percent of planets")
        plt.axvline(x=4.17, label="1h Instrument: " +
                    str(inper(len(detected[detected.p_orb < 4.17])/len(detected))) + "%", color="red")
        plt.axvline(x=50, label="12h Instrument:" +
                    str(inper(len(detected[detected.p_orb < 50])/len(detected))) + "%", color="orange")
        plt.axvline(x=100, label="24h Instrument:" +
                    str(inper(len(detected[detected.p_orb < 100])/len(detected))) + "%", color="yellow")
        # plt.ylim(0, 5000)
        plt.xlim(-1, 150)
        plt.legend()
        plt.tight_layout()
        plt.show()
    # Average SNR per period histogram
    if False:
        plt.figure()
        p_bins = np.arange(min_period, max_period, 5)
        avg_snr = []
        for small, large in zip(p_bins[:-1], p_bins[1:]):
            avg_snr.append(
                np.mean(dt[(dt["p_orb"] < large) & (dt["p_orb"] > small)]["snr_1h"]))
        avg_snr = np.array(avg_snr)
        # plt.bar(p_bins[1:],avg_snr,5,label="time")
        plt.plot(p_bins[1:], avg_snr, label="time")

        avg_snr = []
        for small, large in zip(p_bins[:-1], p_bins[1:]):
            avg_snr.append(
                np.mean(do[(do["p_orb"] < large) & (do["p_orb"] > small)]["snr_1h"]))
        avg_snr = np.array(avg_snr)
        # plt.bar(p_bins[1:],avg_snr,5, alpha =0.5,label="og")
        plt.plot(p_bins[1:], avg_snr, label="og")
        plt.legend()
        plt.show()
    #  EDGE_On Face_On SNR Thersholhd Bar plot
    if False:
        plt.figure()
        face_on_a = do_d[(do_d["inc_p"] <= 0.6)]
        face_on_c = do_d[do_d["inc_p"] >= (np.pi - 0.6)]
        edge_on = do_d[(do_d["inc_p"] >= (np.pi/2 - 0.3))
                       & (do_d["inc_p"] >= (np.pi/2 + 0.3))]
        counts = [len(face_on_a), len(face_on_c), len(edge_on)]
        plt.bar([0, 1, 2], counts, 0.25, tick_label=[
                "Anticlockwise", "Clockwise", "Edge_on"], label="original")

        face_on_a = dt_d[(dt_d["inc_p"] <= 0.2)]
        face_on_c = dt_d[dt_d["inc_p"] >= (np.pi - 0.2)]
        edge_on = dt_d[(dt_d["inc_p"] >= (np.pi/2 - 0.02))
                       & (dt_d["inc_p"] >= (np.pi/2 + 0.02))]
        counts = [len(face_on_a), len(face_on_c), len(edge_on)]
        plt.bar([0.25, 1.25, 2.25], counts, 0.25, tick_label=[
                "Anticlockwise", "Clockwise", "Edge_on"], label="Time")
        plt.ylabel("SNR_12h > 7")
        plt.legend()
    # Bar plot improvments in SNR
    if False:
        plt.figure()
        difference = dt.copy()
        difference["snr_1h"] = dt["snr_1h"] - do["snr_1h"]
        face_on_a = difference[(difference["inc_p"] <= 0.6)]
        face_on_c = difference[difference["inc_p"] >= (np.pi - 0.6)]
        edge_on = difference[(difference["inc_p"] >= (np.pi/2 - 0.3))
                             & (difference["inc_p"] >= (np.pi/2 + 0.3))]
        counts = np.array([len(face_on_a[face_on_a["snr_1h"] > 0]), len(
            face_on_c[face_on_c["snr_1h"] > 0]), len(edge_on[edge_on["snr_1h"] > 0])])/len(dt)

        plt.bar([-0.15, 0.85, 1.85], counts, 0.3,
                label="improvments", color="green")
        counts = np.array([len(face_on_a[face_on_a["snr_1h"] < 0]), len(
            face_on_c[face_on_c["snr_1h"] < 0]), len(edge_on[edge_on["snr_1h"] < 0])])/len(dt)

        plt.bar([0.15, 1.15, 2.15], counts, 0.3,
                label="decreases", color="red")
        plt.xticks([0, 1, 2], labels=[
            "Anticlockwise", "Clockwise", "Edge_on"])
        plt.ylabel("Counts in Percent")
        plt.title("Changes in SNR for certain inclinations")
        plt.legend()
        plt.show()
    # SNR Distribution
    if False:
        snr_min = -1  # min((min(dt.snr_1h), min(do.snr_1h)))
        snr_max = 20  # max((max(dt.snr_1h), max(do.snr_1h)))
        snr_bins = np.arange(snr_min, snr_max, 0.01)
        plt.figure()
        plt.hist(do.snr_1h, snr_bins, label="original")
        plt.hist(dt.snr_1h, snr_bins, label="time", alpha=0.3)
        plt.legend()
        plt.show()
    if False:
        rng = 100
        steps = 1
        diff_snr = dt.snr_1h - do.snr_1h
        max_diff = max(diff_snr)
        max_diff = max_diff if max_diff <= 1000 else 1000
        min_diff = min(diff_snr)
        min_diff = min_diff if min_diff >= -1000 else -1000
        rng = max([max_diff,-min_diff])
        snr_bins = np.arange(-rng,rng, steps)
        plt.figure()
        _, _, patches = plt.hist(diff_snr,snr_bins,density=True)
        for i in range(len(patches)):
            pos = i*steps -rng
            if  pos < 0:
                patches[i].set_fc("red")
            else:
                patches[i].set_fc("green")
    if False:
        diff_snr = dt.snr_1h - do.snr_1h
        bins = np.arange(0,180,1)
        counts_p_1 = np.zeros((180)) 
        counts_p_2 = np.zeros((180)) 
        counts_p_10 = np.zeros((180))
        counts_p_50 = np.zeros((180))  
        counts_n_1 = np.zeros((180)) 
        counts_n_2 = np.zeros((180)) 
        counts_n_10 = np.zeros((180))
        counts_n_50 = np.zeros((180))   
        for i in range(len(dt)):
            d_snr = diff_snr[i]
            inc = dt.inc_p[i]/pi*180
            if d_snr>=0:
                if d_snr>=50:
                    counts_p_50[int(inc)] += 1
                elif d_snr>=10:
                    counts_p_10[int(inc)] += 1
                elif d_snr>=2:
                    counts_p_2[int(inc)] += 1
                else:
                    counts_p_1[int(inc)] += 1
                    
            elif d_snr <= 0:
                if d_snr<=-50:
                    counts_n_50[int(inc)] += 1
                elif d_snr<=-10:
                    counts_n_10[int(inc)] += 1
                elif d_snr<=-2:
                    counts_n_2[int(inc)] += 1
                else:
                    counts_n_1[int(inc)] += 1
        tot = counts_p_1+counts_p_2+counts_p_10+counts_p_50+counts_n_1+counts_n_2+counts_n_10+counts_n_50
        plt.figure()
        plt.suptitle("Changes in SNR depending on inclination")
        plt.subplot(1,2,1)
        plt.bar(bins,counts_p_1,1,label = "+ 0-2")
        plt.bar(bins,counts_p_2,1,bottom=counts_p_1,label = "+ 2-10")
        plt.bar(bins,counts_p_10,1,bottom=counts_p_1 + counts_p_2,label = "+ 10-50")
        plt.bar(bins,counts_p_50,1,bottom=counts_p_1 + counts_p_2+counts_p_10,label = "+ 50+")
        plt.bar(bins,-counts_n_1,1,label = "- 0-2")
        plt.bar(bins,-counts_n_2,1,bottom=-counts_n_1,label = "- 2-10")
        plt.bar(bins,-counts_n_10,1,bottom=-counts_n_1 - counts_n_2,label = "- 10-50")
        plt.bar(bins,-counts_n_50,1,bottom=-counts_n_1 - counts_n_2 - counts_n_10,label = "- 50+")
        plt.ylabel("Counts")
        plt.xlabel("Degrees of inclination")
        plt.subplot(1,2,2)
        plt.bar(bins,counts_p_1/tot,1,label = "+ 0-2")
        plt.bar(bins,counts_p_2/tot,1,bottom=counts_p_1/tot,label = "+ 2-10")
        plt.bar(bins,counts_p_10/tot,1,bottom=(counts_p_1 + counts_p_2)/tot,label = "+ 10-50")
        plt.bar(bins,counts_p_50/tot,1,bottom=(counts_p_1 + counts_p_2+counts_p_10)/tot,label = "+ 50+")
        plt.bar(bins,-counts_n_1/tot,1,label = "- 0-2")
        plt.bar(bins,-counts_n_2/tot,1,bottom=-counts_n_1/tot,label = "- 2-10")
        plt.bar(bins,-counts_n_10/tot,1,bottom=-(counts_n_1 + counts_n_2)/tot,label = "- 10-50")
        plt.bar(bins,-counts_n_50/tot,1,bottom=-(counts_n_1 + counts_n_2 + counts_n_10)/tot,label = "- 50+")
        plt.legend(loc="center")
        plt.ylabel("Percentage")
        plt.xlabel("Degrees of inclination")
        #plt.plot(bins, counts_p- counts_n,  "-", color = "gray")
# %%
def get_db(file):
    bus = ls.Bus()
    bus.data.options.set_scenario('baseline')
    bus.data.import_catalog(
        "./output/"+file+".hdf5")
    return bus.data.catalog.copy()
# %%