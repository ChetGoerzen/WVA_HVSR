import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyasdf

import hvsrpy
from hvsrpy import utils

# Time domain settings
# --------------------------------------------------------------------------- #
windowlength = 120
filter_bool = False # No filtering
filter_flow = 0.1
filter_fhigh = 10
filter_order = 5

# Width of cosine taper
width = 0.1
# --------------------------------------------------------------------------- #

# Frequency domain settings
# --------------------------------------------------------------------------- #
# Konno-Ohmachi smoothing constant
bandwidth = 40

resample_fmin = 0.1
resample_fmax = 10
resample_fnum = 200
resample_type = 'log'

# Upper and lower frequency limits to restrict peak selection. To use the entire range use `None`.
peak_f_lower = None
peak_f_upper = None
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

# Load data
# --------------------------------------------------------------------------- #

data_path = "/mnt/readynas5/cgoerzen/NZ3DFWI/wfs/CUBES_data/3CASDF/2018_07_01_00_00_00T2018_07_02_00_00_00.h5"
station = "3C_BS20"

with pyasdf.ASDFDataSet(data_path, mpi=False, mode="r") as ds:

    dt = ds.waveforms[station]["hhz_00"][0].stats.delta

    hhz = {
        "amplitude" : ds.waveforms[station]["hhz_00"][0].data,
        "dt" : dt
    }

    hhe = {
        "amplitude" : ds.waveforms[station]["hhe_00"][0].data,
        "dt" : dt
    }

    hhn = {
        "amplitude" : ds.waveforms[station]["hhn_00"][0].data,
        "dt" : dt
    }
    data_dict = {
        "vt" : hhz,
        "ew" : hhe,
        "ns" : hhn
    }
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
file_name_out = "test.hv"
hv.to_file(file_name_out, distribution_f0, distribution_mc, data_format="hvsrpy")
print("Results saved successfully!")
