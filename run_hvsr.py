import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyasdf
import glob

import hvsrpy
from hvsrpy import utils

# Time domain settings
# --------------------------------------------------------------------------- #
windowlength = 1200
filter_bool = True # No filtering
filter_flow = 0.1
filter_fhigh = 9.9
filter_order = 5

# Width of cosine taper
width = 0.1
# --------------------------------------------------------------------------- #

# Frequency domain settings
# --------------------------------------------------------------------------- #
# Konno-Ohmachi smoothing constant
bandwidth = 40

resample_fmin = 0.05
resample_fmax = 10
resample_fnum = 200
resample_type = 'log'

# Upper and lower frequency limits to restrict peak selection. To use the entire range use `None`.
peak_f_lower = 0.1
peak_f_upper = 10
# --------------------------------------------------------------------------- #

# HVSR settings
# --------------------------------------------------------------------------- #
method = "geometric-mean"
azimuth = 0

rejection_bool = True # Cox et al., 2020 frequency domain rejection
n = 2 # Number of standard deviations to consider when performing rejection
max_iterations = 50 # Max number of iterations for rejection
distribution_f0 = "lognormal" # Distribution of f0
distribution_mc = "lognormal" # Distribution of mean curve
# --------------------------------------------------------------------------- #

# Plot settings

# --------------------------------------------------------------------------- #
ymin = 0
ymax = 10
# --------------------------------------------------------------------------- #

# Load data
# --------------------------------------------------------------------------- #

data_path = "/mnt/readynas5/cgoerzen/NZ3DFWI/wfs/CUBES_data/3CASDF/*.h5"
data_files = glob.glob(data_path)
station = "3C_BS04"
tmp_station = station.replace("_", ".")
for data_file in data_files:

    date = data_file.split("/")[-1].split(".")[0]

    try:
        with pyasdf.ASDFDataSet(data_file, mpi=False, mode="r") as ds:
            print(ds)
            if tmp_station in ds.waveforms.list():
                if len(ds.waveforms[station].list()) > 3: # Make sure that all three channels aand metadata is there
                    dt = ds.waveforms[station]["hhz_00"][0].stats.delta

                    hhz = {
                        "amplitude" : ds.waveforms[station]["hhz_00"][0].data + 1e-6,
                        "dt" : dt
                    }

                    hhe = {
                        "amplitude" : ds.waveforms[station]["hhe_00"][0].data + 1e-6,
                        "dt" : dt
                    }

                    hhn = {
                        "amplitude" : ds.waveforms[station]["hhn_00"][0].data + 1e-6,
                        "dt" : dt
                    }

                    data_dict = {
                        "vt" : hhz,
                        "ew" : hhe,
                        "ns" : hhn
                    }
                else:
                    print("Not enough components!")
                    continue
            else:
                print(f"No data on station: {station} for date: {date}")
                continue
            # --------------------------------------------------------------------------- #

        sensor = hvsrpy.Sensor3c.from_dict(data_dict)
        bp_filter = {
                     "flag":filter_bool, "flow":filter_flow,
                     "fhigh":filter_fhigh, "order":filter_order
                }

        resampling = {"minf":resample_fmin, "maxf":resample_fmax,
                      "nf":resample_fnum, "res_type":resample_type}

        hv = sensor.hv(windowlength, bp_filter, width, bandwidth, resampling, method,
                       f_low=peak_f_lower, f_high=peak_f_upper, azimuth=azimuth)

        reliability = utils.sesame_reliability(hv.meta["Window Length"],
                                               len(hv.valid_window_indices), hv.frq,
                                               hv.mean_curve(), hv.std_curve(),
                                               search_limits=(peak_f_lower, peak_f_upper),
                                               verbose=1)

        clarity = utils.sesame_clarity(hv.frq, hv.mean_curve(), hv.std_curve(),
                                       hv.std_f0_frq(distribution="normal"),
                                       search_limits=(peak_f_lower, peak_f_upper),
                                       verbose=1)

        # Do the plotting
        # ------------------------------------------------------------------- #
        fig = plt.figure(figsize=(6,6), dpi=150)
        gs = fig.add_gridspec(nrows=6,ncols=6)

        ax0 = fig.add_subplot(gs[0:2, 0:3])
        ax1 = fig.add_subplot(gs[2:4, 0:3])
        ax2 = fig.add_subplot(gs[4:6, 0:3])

        if rejection_bool:
            ax3 = fig.add_subplot(gs[0:3, 3:6])
            ax4 = fig.add_subplot(gs[3:6, 3:6])
        else:
            ax3 = fig.add_subplot(gs[0:3, 3:6])
            ax4 = False

        individual_width = 0.3
        median_width = 1.3
        for ax, title in zip([ax3, ax4], ["Before Rejection", "After Rejection"]):
            # Rejected Windows
            if title=="After Rejection":
                if len(hv.rejected_window_indices):
                    label = "Rejected"
                    for amp in hv.amp[hv.rejected_window_indices]:
                        ax.plot(hv.frq, amp, color='#00ffff', linewidth=individual_width, zorder=2, label=label)
                        label=None

            # Accepted Windows
            label="Accepted"
            for amp in hv.amp[hv.valid_window_indices]:
                ax.plot(hv.frq, amp, color='#888888', linewidth=individual_width,
                        label = label if title=="Before Rejection" else "")
                label=None

            # Window Peaks
            ax.plot(hv.peak_frq, hv.peak_amp, linestyle="", zorder=2,
                    marker='o', markersize=2.5, markerfacecolor="#ffffff", markeredgewidth=0.5, markeredgecolor='k',
                    label="" if title=="Before Rejection" and rejection_bool else r"$f_{0,i}$")

            # Peak Mean Curve
            ax.plot(hv.mc_peak_frq(distribution_mc), hv.mc_peak_amp(distribution_mc), linestyle="", zorder=4,
                    marker='D', markersize=4, markerfacecolor='#66ff33', markeredgewidth=1, markeredgecolor='k',
                    label = "" if title=="Before Rejection" and rejection_bool else r"$f_{0,mc}$")

            # Mean Curve
            label = r"$LM_{curve}$" if distribution_mc=="lognormal" else "Mean"
            ax.plot(hv.frq, hv.mean_curve(distribution_mc), color='k', linewidth=median_width,
                    label="" if title=="Before Rejection" and rejection_bool else label)

            # Mean +/- Curve
            label = r"$LM_{curve}$"+" ± 1 STD" if distribution_mc=="lognormal" else "Mean ± 1 STD"
            ax.plot(hv.frq, hv.nstd_curve(-1, distribution_mc),
                    color='k', linestyle='--', linewidth=median_width, zorder=3,
                    label = "" if title=="Before Rejection" and rejection_bool else label)
            ax.plot(hv.frq, hv.nstd_curve(+1, distribution_mc),
                    color='k', linestyle='--', linewidth=median_width, zorder=3)

            # f0 +/- STD
            if ymin is not None and ymax is not None:
                ax.set_ylim((ymin, ymax))
            label = r"$LM_{f0}$"+" ± 1 STD" if distribution_f0=="lognormal" else "Mean f0 ± 1 STD"
            _ymin, _ymax = ax.get_ylim()
            ax.plot([hv.mean_f0_frq(distribution_f0)]*2, [_ymin, _ymax], linestyle="-.", color="#000000")
            ax.fill([hv.nstd_f0_frq(-1, distribution_f0)]*2 + [hv.nstd_f0_frq(+1, distribution_f0)]*2, [_ymin, _ymax, _ymax, _ymin], 
                    color = "#ff8080",
                    label="" if title=="Before Rejection" and rejection_bool else label)
            ax.set_ylim((_ymin, _ymax))

            ax.set_xscale('log')
            ax.set_xlabel("Frequency (Hz)")
            ax.set_ylabel("HVSR Amplitude")
            if rejection_bool:
                if title=="Before Rejection":
                    print("\nStatistics before rejection:")
                    #hv.print_stats(distribution_f0)
                    c_iter = hv.reject_windows(n, max_iterations=max_iterations,
                                               distribution_f0=distribution_f0, distribution_mc=distribution_mc)
                elif title=="After Rejection":
                    fig.legend(ncol=4, loc='lower center', bbox_to_anchor=(0.51, 0), columnspacing=2)

                    print("\nAnalysis summary:")
                    print(pd.DataFrame(columns=[""], index=["Window length", "No. of windows", "Number of iterations to convergence", "No. of rejected windows"],
                            data=[f"{windowlength}s", str(sensor.ns.nseries), f"{c_iter} of {max_iterations} allowed", str(sum(hv.rejected_window_indices))]))
                    print("\nStatistics after rejection:")
                    #hv.print_stats(distribution_f0)
            else:
                print(pd.DataFrame(columns=[""], index=["Window length", "No. of windows"],
                                 data=[f"{windowlength}s", str(sensor.ns.nseries)]))
                #hv.print_stats(distribution_f0)
                fig.legend(loc="upper center", bbox_to_anchor=(0.77, 0.4))
                break
            ax.set_title(title)

        norm_factor = sensor.normalization_factor
        for ax, timerecord, name in zip([ax0,ax1,ax2], [sensor.ns, sensor.ew, sensor.vt], ["NS", "EW", "VT"]):
            ctime = timerecord.time
            amp = timerecord.amp/norm_factor
            ax.plot(ctime.T, amp.T, linewidth=0.2, color='#888888')
            ax.set_title(f"Time Records ({name})")
            ax.set_yticks([-1, -0.5, 0, 0.5, 1])
            ax.set_xlim(0, windowlength*timerecord.nseries)
            ax.set_ylim(-1, 1)
            ax.set_xlabel('Time (s)')
            ax.set_ylabel('Normalized Amplitude')
            ax.plot(ctime[hv.rejected_window_indices].T, amp[hv.rejected_window_indices].T, linewidth=0.2, color="cyan")

        if rejection_bool:
            axs = [ax0, ax3, ax1, ax4, ax2]
        else:
            axs = [ax0, ax3, ax1, ax2]

        for ax, letter in zip(axs, list("abcde")):
            ax.text(0.97, 0.97, f"({letter})", ha="right", va="top", transform=ax.transAxes, fontsize=12)
            for spine in ["top", "right"]:
                ax.spines[spine].set_visible(False)


        fig.tight_layout(h_pad=1, w_pad=2, rect=(0,0.08,1,1))
        plt.show()
        # ------------------------------------------------------------------- #


        figure_name_out = "plots/" + date + "_" + station + ".png"

        fig.savefig(figure_name_out, dpi=300, bbox_inches='tight')
        plt.close()
        print("Figure saved successfully!")

        file_name_out = "results/" + date + "_" + station + ".hv"
        hv.to_file(file_name_out, distribution_f0, distribution_mc, data_format="hvsrpy")
        print("Results saved successfully! " + date)
    except OSError:
        print("Cannot open file: " + data_file)
