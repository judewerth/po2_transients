"""
po2_transients.py: DATA layer of the event-triggered pO2 transient pipeline.

This module turns `pO2_all_synced_runXX.mat` files into pandas DataFrames, and it holds every computation
done on those tables (averages, binning, mixed models). It does NOT plot: see po2_plots.py for figures and
po2_gui.py for the interactive viewers. Dependencies run one way only:

    po2_transients (tables)  <--  po2_plots (figures)  <--  po2_gui (widgets)

Data flow for one run
=====================

    .mat file
       |  STAGE 1  load_synced_run         raw event-source signals (pupil and/or locomotion) + pO2 matrix
       v
    sig, po2                                sig = {t_po2, fs_po2, t_<src>, <src>, fs_<src>}, po2 = (n_time, n_points)
       |  STAGE 2  preprocess_event_source EVENT_METHODS[method].preprocess -> filtered signals (ext)
       v                                    pupil: moving median -> low-pass -> derivative;  locom: |x - mean| -> Gaussian
    ext                                     ext = {t_s, s} for every signal s of the method
       |  STAGE 3  detect_events           peaks of ext[detect_on] above a threshold + EVENT_METRICS per event
       v
    events dict                             event_time, event_bounds, event_height, event_halfwidth, event_auc
       |  STAGE 4  extract_event_segments  fixed window around every event, from every signal and every pO2 point
       v
    segments                                (n_events, n_window), NaN-padded where the window leaves the recording
       |  STAGE 5  pO2 metrics             POINT_METRICS (whole run, per point), PO2_EVENT_METRICS (per event x point)
       v
    TransientData                           the five tables below (built by analyze_run)

`load_run` runs stages 1-2 (slow: file I/O) and `analyze_run` runs stages 3-5 (fast), so the GUI can cache the
first half and re-run the second whenever a threshold changes. `process_run` = both, and `process_all_runs`
stacks every run into one TransientData.

Tables (pandas DataFrames; array-valued cells hold whole traces)
================================================================
Every table starts with RUN_KEYS = [session, mouse_id, date, run]. Join tables on TABLE_KEYS[table].
Below, s is each signal in EVENT_METHODS[method].signals (pupil: pupil, dxdt_pupil; locom: locom).

    table            one row per          key columns         other columns
    ---------------  -------------------  ------------------  --------------------------------------------------
    dataset          run                  RUN_KEYS            t_po2 (n_time,), fs_po2, raw t_<src>, <src>, fs_<src>
    po2_signals      run x pO2 point      POINT_KEYS          t_po2, PO2_VARIANTS traces (po2, abs_po2, norm_po2),
                                                              POINT_METRICS scalars (rms, cv, baseline, nan_frac)
    extract_signals  run                  RUN_KEYS            filtered t_s, s for every signal s; event_threshold
    events           event                EVENT_KEYS          event_time (s), event_bounds [start, end] (s),
                                                              EVENT_METRICS scalars, and per signal s:
                                                              t_s_event (time rel. to event, s), s_event_trace
    po2_events       event x pO2 point    EVENT_POINT_KEYS    t_event (rel. time, s), event_time, event_trace (raw),
                                                              event_trace_norm ((trace - rms) / rms),
                                                              PO2_EVENT_METRICS scalars (event_peak, ...)

pO2 "points" are measurement locations along the vessel (columns of pO2OverTime). A point index only has a
meaning within one run: point 3 of two runs is not the same vessel location.

Two pO2 normalizations exist, and they are not the same:
    norm_po2         (po2_signals) = po2 / baseline          (baseline = the point's full-run mean)
    event_trace_norm (po2_events)  = (trace - rms) / rms     (rms = the point's full-run RMS); event_peak uses it.
    Because RMS >= mean, (trace - rms) / rms sits slightly below 0 on average (by 1 - mean/rms).

How to extend (registries)
==========================
Adding an entry to a registry is enough: the tables, BIN_METRICS, the plots and the GUI dropdowns all read
from these dicts.

* New event-source signal in the .mat file -> add to RAW_SIGNALS: {name: (time variable, value variable)}.
* New event detection method ->
      1. write preprocess_<name>(sig, params) -> {"t_<s>": t, "<s>": x, ...} for each signal it produces
         (put any new parameters in EventParams),
      2. register EVENT_METHODS["<name>"] = EventMethod(raw=[...], signals=[...], detect_on="<s>",
         preprocess=preprocess_<name>, describe=lambda p: "...").
* New per-event metric (one scalar per event, computed from the detection signal) ->
      write fn(ctx: DetectionContext) -> (n_events,) and add EVENT_METRICS["<name>"] = fn.
* New per-point metric (one scalar per run x point, from the whole pO2 trace) ->
      write fn(po2 (n_time, n_points)) -> (n_points,) and add POINT_METRICS["<name>"] = fn.
* New per-event x point pO2 metric (from the normalized event window) ->
      write fn(t_rel (n_window,), norm (n_events, n_window)) -> (n_events,) and add PO2_EVENT_METRICS["<name>"] = fn.
* New run-level pO2 trace variant -> PO2_VARIANTS["<name>"] = fn(po2, point_metrics) -> (n_time, n_points).

Known issues (left unchanged on purpose so results stay reproducible; see the NOTE comments)
=============================================================================================
A. extract_event_segments labels sample t.searchsorted(t_event) as t_rel = 0, but that sample lies up to one
   sample period AFTER the event. pO2 windows (fs as low as 0.33 Hz) can therefore look up to ~3 s early.
B. DEFAULT_BINS [[0, 33], [34, 67], [68, 100]] leaves the gaps P33-34 and P67-68 unbinned.
C. The default min_distance_s (30 s) is shorter than the window (60 s), so neighbouring event windows overlap
   and can share the same pO2 peak sample.
"""

import os
import re
import warnings
from dataclasses import dataclass, asdict, fields, field
from typing import Callable

import h5py
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from patsy import dmatrix
from scipy.io import loadmat
from scipy.ndimage import median_filter, gaussian_filter1d
from scipy.signal import butter, filtfilt, find_peaks
from scipy.stats import norm


# =============================================================================
# 0. Where the data lives
# =============================================================================
# Session folders (e.g. "26-10-05_B6_318") are looked up inside DATA_DIR. By default that is the data/ folder
# next to this file. To keep the data somewhere else (lab drive, external disk), set PO2_DATA_DIR in a .env file
# at the repo root (copy .env.example) or as an environment variable. A relative path is taken relative to the
# repo root. This only changes WHERE files are read from, never what is computed.

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
try:
    from dotenv import load_dotenv           # python-dotenv: copies KEY=value lines from .env into os.environ
    load_dotenv(os.path.join(REPO_ROOT, ".env"))
except ImportError:                          # .env support is optional; environment variables still work
    pass
DATA_DIR = os.path.join(REPO_ROOT, os.path.expanduser(os.environ.get("PO2_DATA_DIR") or "data"))


# =============================================================================
# 1. Keys: the columns that identify rows and join the tables
# =============================================================================

RUN_KEYS = ["session", "mouse_id", "date", "run"]   # one recording run
EVENT_KEYS = RUN_KEYS + ["event"]                   # one detected event within a run (event = 0, 1, ...)
POINT_KEYS = RUN_KEYS + ["point"]                   # one pO2 measurement point within a run
EVENT_POINT_KEYS = EVENT_KEYS + ["point"]           # one event x one pO2 point

# table -> the key columns that identify one of its rows (every table has one row per key combination)
TABLE_KEYS = {
    "dataset": RUN_KEYS,
    "po2_signals": POINT_KEYS,
    "extract_signals": RUN_KEYS,
    "events": EVENT_KEYS,
    "po2_events": EVENT_POINT_KEYS,
}


# =============================================================================
# 2. Parameters
# =============================================================================

@dataclass
class EventParams:
    """
    All tunable parameters of the pipeline, in one place. The stage that uses each field is in its comment.
    Notebooks build one of these. The GUI copies it with `dataclasses.replace` when a control changes.
    """
    event_method: str = "pupil"            # key of EVENT_METHODS: which signal events are detected on (stage 2-3)
    rel_to_rms: bool = True                # stage 3: True -> threshold = height_thresh * RMS(detection signal)
    locom_sigma_s: float = 2.0             # stage 2 (locom): Gaussian smoothing sigma (s)
    movmedian_len_s: float = 6.0           # stage 2 (pupil): moving-median window on the raw pupil (s)
    lpf_cutoff_hz: float = 0.1             # stage 2 (pupil): 4th-order Butterworth low-pass cutoff (Hz)
    neg_peaks: bool = False                  # stage 3: whether to detect negative peaks as well
    height_thresh: float = 1.5             # stage 3: threshold, in x RMS (rel_to_rms) or absolute units
    min_distance_s: float = 30.0           # stage 3: minimum spacing between events (s); NOTE (known issue C)
    extract_window_s: tuple = (-30, 30)    # stage 4: window cut around each event (s, relative to the event)

    def describe(self):
        """One-line human summary of the parameters (used in figure titles and the GUI status line)."""
        filt = EVENT_METHODS[self.event_method].describe(self)
        thresh = f"{self.height_thresh}xRMS" if self.rel_to_rms else f"{self.height_thresh}"
        return (f"{filt} -> Peak thresh: {thresh}, min sep {self.min_distance_s} s -> "
                f"Extract {self.extract_window_s} s")

    def as_dict(self):
        return asdict(self)


# =============================================================================
# 3. Data container
# =============================================================================

@dataclass
class TransientData:
    """
    The five analysis tables (schemas in the module docstring). Every table begins with RUN_KEYS.
    Built by `analyze_run` (one run) or `process_all_runs` (stacked runs).
    """
    dataset: pd.DataFrame
    po2_signals: pd.DataFrame
    extract_signals: pd.DataFrame
    events: pd.DataFrame
    po2_events: pd.DataFrame

    @classmethod
    def concat(cls, parts):
        """Stack the tables of several TransientData (e.g. one per run) into one."""
        return cls(**{f.name: pd.concat([getattr(p, f.name) for p in parts], ignore_index=True)
                      for f in fields(cls)})

    def select_run(self, session, run):
        """A TransientData restricted to one run (used by the single-run plots)."""
        def pick(df):
            return df[(df["session"] == session) & (df["run"] == run)].reset_index(drop=True)
        return TransientData(**{f.name: pick(getattr(self, f.name)) for f in fields(self)})

    def summary(self):
        """Row count and grain of each table."""
        grain = {"dataset": "run", "po2_signals": "run x point", "extract_signals": "run",
                 "events": "event", "po2_events": "event x point"}
        return pd.DataFrame([{"table": f.name, "rows": len(getattr(self, f.name)), "one row per": grain[f.name]}
                             for f in fields(self)])


# =============================================================================
# 4. STAGE 1: load (MATLAB file -> raw numpy arrays)
# =============================================================================

# event-source signal name -> (time variable, value variable) in pO2_all_synced_runXX.mat
RAW_SIGNALS = {
    "pupil": ("tPupilShort", "pupilShort"),
    "locom": ("tAccelShort", "AccelShort"),
    "sleep_state": ("t_img", "frame_state")
}

# --- generic .mat readers (scipy for pre-v7.3 files, h5py for v7.3+) --------------------------------------

def _ndarray_to_list(arr):
    if not isinstance(arr, np.ndarray):
        return arr
    if arr.ndim == 0:
        return _to_python(arr.item())
    return [_ndarray_to_list(sub) for sub in arr]


def _to_python(obj):
    """Convert scipy-loaded MATLAB objects (pre-v7.3) into Python types."""
    if hasattr(obj, "_fieldnames"):  # MATLAB struct
        return {name: _to_python(getattr(obj, name)) for name in obj._fieldnames}
    if isinstance(obj, np.ndarray):
        if obj.dtype == object:       # cell/object arrays
            return _ndarray_to_list(obj)
        if obj.ndim == 0:
            return _to_python(obj.item())
        return obj
    return obj


def _is_matlab_char_dataset(dset):
    cls = dset.attrs.get("MATLAB_class", b"")
    if isinstance(cls, bytes):
        cls = cls.decode("utf-8", errors="ignore")
    elif isinstance(cls, np.ndarray) and cls.size == 1:
        cls = cls.item()
        if isinstance(cls, bytes):
            cls = cls.decode("utf-8", errors="ignore")
    return cls == "char"


def _h5_to_python(node, h5file):
    """Convert HDF5 MATLAB objects (v7.3+) into Python types."""
    if isinstance(node, h5py.Group):
        return {k: _h5_to_python(node[k], h5file) for k in node.keys() if k != "#refs#"}

    if isinstance(node, h5py.Dataset):
        data = node[()]
        # MATLAB char arrays are stored as uint16
        if isinstance(data, np.ndarray) and data.dtype == np.uint16 and _is_matlab_char_dataset(node):
            return "".join(map(chr, data.flatten(order="F"))).rstrip("\x00")
        # object references (cells / struct links)
        if isinstance(data, np.ndarray) and data.dtype.kind == "O":
            out = np.empty(data.shape, dtype=object)
            it = np.nditer(data, flags=["multi_index", "refs_ok"], op_flags=["readonly"])
            for ref in it:
                r = ref.item()
                out[it.multi_index] = _h5_to_python(h5file[r], h5file) if r else None
            return _ndarray_to_list(out)
        if isinstance(data, h5py.Reference):
            return _h5_to_python(h5file[data], h5file) if data else None
        if isinstance(data, bytes):
            return data.decode("utf-8", errors="ignore")
        if isinstance(data, np.ndarray) and data.ndim == 0:
            return data.item()
        return data
    return node


def loadmat_nested_dict(path):
    """Load a MATLAB .mat file (any version) into nested Python dicts / lists / numpy arrays."""
    try:
        raw = loadmat(path, struct_as_record=False, squeeze_me=True)
        raw = {k: v for k, v in raw.items() if not k.startswith("__")}
        return {k: _to_python(v) for k, v in raw.items()}
    except NotImplementedError as e:
        if "v7.3" not in str(e):
            raise
        with h5py.File(path, "r") as f:
            return {k: _h5_to_python(v, f) for k, v in f.items() if k != "#refs#"}


# --- pipeline-specific loading -----------------------------------------------------------------------------

def parse_session(folder):
    """Split a session folder name 'YY-MM-DD_MOUSE' into (date, mouse_id)."""
    m = re.fullmatch(r"(\d{2}-\d{2}-\d{2})_(.+)", os.path.basename(os.path.normpath(folder)))
    if not m:
        raise ValueError(f"Folder name {folder} does not match expected pattern 'YY-MM-DD_MOUSE'")
    return m.group(1), m.group(2)


def run_keys(folder, run):
    """The RUN_KEYS values of one run, as a dict (prepended to every table row of that run)."""
    date, mouse_id = parse_session(folder)
    return {"session": folder, "mouse_id": mouse_id, "date": date, "run": run}


def synced_run_path(folder, run):
    """Path to the synced run file, e.g. <DATA_DIR>/<folder>/pO2_all_synced_run04.mat (an absolute folder is used as is)."""
    path = os.path.join(DATA_DIR, folder, f"pO2_all_synced_run{run:02d}.mat")
    if not os.path.exists(path):
        raise FileNotFoundError(f"File {path} does not exist")
    return path


def load_synced_run(folder, run, raw_names):
    """
    STAGE 1. Load the pO2 time series and the requested raw event-source signals of one run.

    Args:
        raw_names : RAW_SIGNALS keys to load (normally EVENT_METHODS[method].raw)
    Returns:
        sig (dict) : t_po2 (n_time,), fs_po2 (Hz), and for each src in raw_names: t_<src>, <src>, fs_<src>
                     -> stored as the `dataset` table row; input of stage 2
        po2 (n_time, n_points) : pO2 of every measurement point (time base sig["t_po2"]) -> stages 4-5
    """
    raw = loadmat_nested_dict(synced_run_path(folder, run))

    t_po2 = np.asarray(raw["tImaging"]).flatten()
    po2 = np.asarray(raw["pO2OverTime"])
    if po2.ndim == 1:                  # single-point runs come back 1-D
        po2 = po2[:, np.newaxis]
    sig = {"t_po2": t_po2, "fs_po2": 1 / np.mean(np.diff(t_po2))}

    for src in RAW_SIGNALS:            # registry order, so column order does not depend on raw_names order
        if src not in raw_names:
            continue
        if src == "sleep_state":
            raw = loadmat_nested_dict(os.path.join(DATA_DIR, folder, f"sleepVikram_run{run:02d}.mat"))["VS"]["classification"]
        t_var, x_var = RAW_SIGNALS[src]
        t = np.asarray(raw[t_var]).flatten()
        sig.update({f"t_{src}": t, src: np.asarray(raw[x_var]).flatten(), f"fs_{src}": 1 / np.mean(np.diff(t))})
    return sig, po2


# =============================================================================
# 5. STAGE 2: preprocess the event source (raw signal -> filtered signals used for detection)
# =============================================================================

def filter_pupil(pupil, fs_pupil, movmedian_len_s, lpf_cutoff_hz):
    """
    Pupil cleaning: a moving median removes blink / tracking spikes, then a zero-phase 4th-order
    Butterworth low-pass smooths what is left. Returns (n,) like `pupil`.
    """
    pupil_med = median_filter(pupil, size=int(movmedian_len_s * fs_pupil))
    b, a = butter(4, lpf_cutoff_hz / (fs_pupil / 2), btype="low")   # cutoff normalized to Nyquist
    return filtfilt(b, a, pupil_med)                                # filtfilt = forward + backward, no phase lag


def pupil_derivative(t_pupil, pupil_filt):
    """
    First difference of the filtered pupil. Positive peaks are dilation onsets.
    Returns (t[1:], dxdt), both (n - 1,). dxdt is per sample, not per second.
    """
    return t_pupil[1:], np.diff(pupil_filt)


def filter_locom(locom, fs_locom, sigma_s):
    """
    Locomotion envelope: remove the mean (a crude high-pass), rectify, then Gaussian-smooth with sigma_s seconds.
    Returns (n,).
    """
    locom = locom - np.mean(locom)
    return gaussian_filter1d(np.abs(locom), sigma=sigma_s * fs_locom)


def filter_sleep_state(sleep_state, fs_sleep_state, min_distance_s):
    """
    Sleep-state cleaning (0 = unknown, 1 = awake, 2 = NREM, 3 = REM).

    A "bout" is a stretch of one nonzero state; 0s inside it do not break it but do not count towards its
    duration (duration = number of nonzero samples / fs).
      - Bouts shorter than min_distance_s -> 0 (unknown).
      - 0s embedded in a kept bout -> that bout's state.
      - 0s between two kept bouts -> the state of the later bout.
      - 0s next to only one kept bout (the other side is a dropped bout or the recording edge) stay 0.
    Single pass: bouts are not re-merged after short bouts between them are dropped.
    Returns (n,) int like `sleep_state`.
    """
    state = np.asarray(sleep_state).astype(int)
    filtered = np.zeros_like(state)
    idx = np.flatnonzero(state)                         # nonzero samples; 0s are transparent
    if idx.size == 0:
        return filtered

    # bouts = runs of equal value in the nonzero samples (0s between equal values are skipped over)
    vals = state[idx]
    starts = np.r_[0, np.flatnonzero(np.diff(vals)) + 1]
    ends = np.r_[starts[1:], vals.size]                 # exclusive, in idx coordinates
    keep = (ends - starts) >= min_distance_s * fs_sleep_state

    prev_kept_end = None                                # sample index just after the last kept bout's last nonzero sample
    for k, s, e in zip(keep, starts, ends):
        first, last = idx[s], idx[e - 1]
        if k:
            # embedded 0s take the bout's state; the 0 gap back to a directly preceding kept bout goes to this (later) bout
            fill_from = prev_kept_end if prev_kept_end is not None else first
            filtered[fill_from:last + 1] = vals[s]
            prev_kept_end = last + 1
        else:
            prev_kept_end = None 
                           # a dropped bout separates the 0s from earlier kept bouts
    return filtered

def sleep_state_derivative(t_sleep_state, sleep_state_filt):
    sleep_state_filt[sleep_state_filt == 3] = 4 # important to distinguish all 6 posibile connections
    dxdt = np.diff(sleep_state_filt)

    # ensure dxdt transitions are not to an unknown state
    for i in np.where(dxdt != 0)[0]:
        if np.any(sleep_state_filt[max(i - 1, 0):i + 2] == 0):
            dxdt[i] = 0
    return t_sleep_state[1:], dxdt


def preprocess_pupil(sig, params):
    """EVENT_METHODS["pupil"].preprocess: filtered pupil and its derivative (events are detected on dxdt_pupil)."""
    pupil_filt = filter_pupil(sig["pupil"], sig["fs_pupil"], params.movmedian_len_s, params.lpf_cutoff_hz)
    t_dxdt, dxdt = pupil_derivative(sig["t_pupil"], pupil_filt)
    return {"t_pupil": sig["t_pupil"], "pupil": pupil_filt, "t_dxdt_pupil": t_dxdt, "dxdt_pupil": dxdt}


def preprocess_locom(sig, params):
    """EVENT_METHODS["locom"].preprocess: smoothed locomotion envelope (events are detected on it directly)."""
    return {"t_locom": sig["t_locom"], "locom": filter_locom(sig["locom"], sig["fs_locom"], params.locom_sigma_s)}

def preprocess_sleep_state(sig, params):
    """EVENT_METHODS["sleep_state"].preprocess: filtered sleep state and its derivative (events are detected on dxdt_sleep_state)."""
    sleep_state_filt = filter_sleep_state(sig["sleep_state"], sig["fs_sleep_state"], params.min_distance_s)
    t_dxdt, dxdt = sleep_state_derivative(sig["t_sleep_state"], sleep_state_filt)
    return {"t_sleep_state": sig["t_sleep_state"], "sleep_state": sleep_state_filt, "t_dxdt_sleep_state": t_dxdt, "dxdt_sleep_state": dxdt}


@dataclass(frozen=True)
class EventMethod:
    """
    One way to detect events (an EVENT_METHODS entry).

    raw        : RAW_SIGNALS loaded for it (stage 1)
    signals    : filtered signals `preprocess` returns; each one is cut around every event (stage 4)
    detect_on  : the signal (one of `signals`) whose peaks are the events (stage 3)
    preprocess : fn(sig, params) -> {"t_<s>": t, "<s>": x for s in signals} (stage 2)
    fs_from    : the raw signal whose sampling rate is the detection rate (sig[f"fs_{fs_from}"])
    describe   : fn(params) -> short text describing the filtering (used by EventParams.describe)
    """
    raw: list
    signals: list
    detect_on: str
    preprocess: Callable
    fs_from: str
    describe: Callable = field(default=lambda params: "")


EVENT_METHODS = {
    "pupil": EventMethod(
        raw=["pupil"], signals=["pupil", "dxdt_pupil"], detect_on="dxdt_pupil",
        preprocess=preprocess_pupil, fs_from="pupil",
        describe=lambda p: f"Movmedian: {p.movmedian_len_s} s -> LPF: {p.lpf_cutoff_hz} Hz"),
    "locom": EventMethod(
        raw=["locom"], signals=["locom"], detect_on="locom",
        preprocess=preprocess_locom, fs_from="locom",
        describe=lambda p: f"|Locom| -> Gaussian sigma: {p.locom_sigma_s} s"),
    "sleep_state": EventMethod(
        raw=["sleep_state"], signals=["sleep_state", "dxdt_sleep_state"], detect_on="dxdt_sleep_state",
        preprocess=preprocess_sleep_state, fs_from="sleep_state",
        describe=lambda p: f"Bouted with min distance: {p.min_distance_s} s (custom)"),
}


def preprocess_event_source(sig, params):
    """
    STAGE 2. Run the preprocessing of `params.event_method` on the raw signals of `load_synced_run`.

    Returns:
        ext (dict) : t_s, s for every signal s of the method -> `extract_signals` table; input of stages 3-4
        fs (float) : sampling rate of the detection signal (Hz), turns seconds into samples in stage 3
    """
    method = EVENT_METHODS[params.event_method]
    return method.preprocess(sig, params), sig[f"fs_{method.fs_from}"]


# =============================================================================
# 6. STAGE 3: detect events (peaks of the detection signal) and measure each one
# =============================================================================

@dataclass
class DetectionContext:
    """
    Everything known about one run's detected peaks. EVENT_METRICS functions read it.

    t, x       : (n,) time (s) and values of the detection signal
    fs         : its sampling rate (Hz)
    peaks      : (n_events,) sample index of each peak
    bounds     : (n_events, 2) sample indices [start, end] of the positive lobe around each peak
    props      : scipy.signal.find_peaks properties (peak_heights, widths, ...)
    """
    t: np.ndarray
    x: np.ndarray
    fs: float
    peaks: np.ndarray
    bounds: np.ndarray
    props: dict


def event_threshold(signal, height_thresh):
    """Absolute threshold = height_thresh * RMS(signal), so it adapts to each run's signal amplitude."""
    return height_thresh * np.sqrt(np.nanmean(signal ** 2))


def positive_lobe_bounds(signal, peaks):
    """
    Sample indices [start, end] of the positive lobe of `signal` that contains each peak.

    start / end are the last sample <= 0 before the peak and the first sample <= 0 after it, or the recording
    edges if there is none. Returns an (n_peaks, 2) int array.
    """
    # every non-positive sample is a possible lobe end, and the recording edges are too. Peaks are > 0, so a
    # peak is never in this list and searchsorted puts it between its two lobe ends.
    edges = np.concatenate([[0], np.flatnonzero(~(signal > 0)), [len(signal) - 1]])
    i = np.searchsorted(edges, peaks)
    return np.column_stack([edges[i - 1], edges[i]]).astype(int)


def find_event_peaks(t, signal, fs, min_height, min_distance_s, rel_to_rms=False, neg_peaks=False):
    """
    STAGE 3a. Locate the events: peaks of `signal` above the threshold, at least min_distance_s apart.

    Returns:
        ctx (DetectionContext) : peaks, lobe bounds and find_peaks properties -> EVENT_METRICS
        threshold (float)      : the absolute threshold used (stored in `extract_signals.event_threshold`)
    """
    if rel_to_rms:
        min_height = event_threshold(signal, min_height)
    # width=0 makes find_peaks compute the width at rel_height (0.5 = half prominence) for event_halfwidth
    peaks, props = find_peaks(signal, height=min_height, distance=int(min_distance_s * fs), width=0, rel_height=0.5)
    if neg_peaks:
        peaks_neg, props_neg = find_peaks(-signal, height=min_height, distance=int(min_distance_s * fs), width=0, rel_height=0.5)
        peaks = np.concatenate([peaks, peaks_neg])

        # ensure neg props are flipped back to original signal
        for k in props:
            if k == "peak_heights":
                props[k] = np.concatenate([props[k], -props_neg[k]])
            else:
                props[k] = np.concatenate([props[k], props_neg[k]])
        sort_idx = np.argsort(peaks)
        peaks = peaks[sort_idx]
        for k in props:
            props[k] = props[k][sort_idx]

    bounds = positive_lobe_bounds(signal, peaks)
    # a lobe can last minutes on a slow signal; cap each side at min_distance_s from the peak
    max_distance = int(min_distance_s * fs)
    bounds[:, 0] = np.maximum(bounds[:, 0], peaks - max_distance)
    bounds[:, 1] = np.minimum(bounds[:, 1], peaks + max_distance)
    return DetectionContext(t=t, x=signal, fs=fs, peaks=peaks, bounds=bounds, props=props), min_height


# --- per-event metrics (one scalar per event; each becomes an `events` column) -----------------------------

def event_height(ctx):
    """Detection-signal value at the peak (e.g. dPupil/dt per sample, or smoothed |accel|)."""
    return ctx.props["peak_heights"]


def event_halfwidth(ctx):
    """Peak width (s) at half its prominence."""
    return ctx.props["widths"] / ctx.fs


def event_auc(ctx):
    """
    Area of the detection signal over the event's positive lobe, in sample units. For dxdt_pupil this equals
    the pupil change over the lobe.
    """
    return np.array([np.trapezoid(ctx.x[s:e + 1]) for s, e in ctx.bounds])


# name -> fn(DetectionContext) -> (n_events,). The order here is the column order in `events`.
EVENT_METRICS = {
    "event_height": event_height,
    "event_halfwidth": event_halfwidth,
    "event_auc": event_auc,
}


def detect_events(t, signal, fs, min_height, min_distance_s, rel_to_rms=False, neg_peaks=False):
    """
    STAGE 3. Find the events (find_event_peaks) and measure each with every EVENT_METRICS entry.

    Returns:
        events (dict of (n_events,) arrays; event_bounds is (n_events, 2)):
            event_time   : peak time (s)
            event_bounds : [start, end] time (s) of the positive lobe containing the peak
            <metric>     : one entry per EVENT_METRICS metric
        threshold (float) : absolute threshold used on `signal`
    """
    ctx, threshold = find_event_peaks(t, signal, fs, min_height, min_distance_s, rel_to_rms, neg_peaks=neg_peaks)
    events = {"event_time": t[ctx.peaks], "event_bounds": t[ctx.bounds].reshape(-1, 2)}
    events.update({name: fn(ctx) for name, fn in EVENT_METRICS.items()})
    return events, threshold


# =============================================================================
# 7. STAGE 4: extract a fixed window around every event
# =============================================================================

def extract_event_segments(t, data, t_events, extract_window_s):
    """
    Cut a fixed-length window around each event from one 1-D signal.

    Samples falling outside the recording are NaN, so every segment has the same length.

    NOTE (known issue A): sample t.searchsorted(t_event) is the first sample AT OR AFTER the event, but it is
    labeled t_rel = 0. The true time of each sample is t_rel + (t[idx0] - t_event), with a lag in [0, 1/fs).
    A fix would interpolate `data` at t_event + offsets / fs instead.

    Args:
        t (n_samples,), data (n_samples,), t_events (n_events,) in s, extract_window_s (w0, w1) in s
    Returns:
        t_rel    (n_window,)           : time relative to the event (s)
        segments (n_events, n_window)
    """
    if data.ndim != 1:
        raise ValueError("data must be 1D")

    fs = 1 / np.median(np.diff(t))
    offsets = np.arange(int(extract_window_s[0] * fs), int(extract_window_s[-1] * fs))   # sample offsets

    segments = np.full((len(t_events), len(offsets)), np.nan)
    for i, t_event in enumerate(t_events):
        idx = t.searchsorted(t_event) + offsets
        valid = (idx >= 0) & (idx < len(t))           # inside the recording; the rest stay NaN
        segments[i, valid] = data[idx[valid]]
    return offsets / fs, segments


# =============================================================================
# 8. STAGE 5: pO2 metrics
# =============================================================================

# --- per-point metrics over the whole run (one scalar per run x point; `po2_signals` columns) --------------
# Each takes po2 (n_time, n_points) and returns (n_points,). All-NaN points give NaN (warnings are silenced
# by po2_point_qc).

def point_rms(po2):
    """Magnitude: sqrt(mean(po2^2)). Used to normalize event windows (event_trace_norm)."""
    return np.sqrt(np.nanmean(po2 ** 2, axis=0))


def point_cv(po2):
    """Noise: coefficient of variation std / mean (NaN where the mean is 0). Used for the point QC cutoff."""
    mean = np.nanmean(po2, axis=0)
    return np.nanstd(po2, axis=0) / np.where(mean != 0, mean, np.nan)


def point_baseline(po2):
    """Mean pO2 of the point over the run. Used to normalize norm_po2 = po2 / baseline."""
    return np.nanmean(po2, axis=0)


def point_nan_frac(po2):
    """Fraction of NaN samples."""
    return np.isnan(po2).mean(axis=0)


# name -> fn(po2) -> (n_points,). The order here is the column order in `po2_signals`.
POINT_METRICS = {
    "rms": point_rms,
    "cv": point_cv,
    "baseline": point_baseline,
    "nan_frac": point_nan_frac,
}


def po2_point_qc(po2):
    """Every POINT_METRICS metric of every point: {name: (n_points,)}."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN points
        return {name: fn(po2) for name, fn in POINT_METRICS.items()}


# --- run-level pO2 trace variants (`po2_signals` trace columns) ---------------------------------------------
# Each takes po2 (n_time, n_points) and qc {POINT_METRICS name: (n_points,)} and returns (n_time, n_points).

def variant_raw(po2, qc):
    """pO2 as recorded."""
    return po2


def variant_abs(po2, qc):
    """|pO2 / baseline (the point's full-run mean)|."""
    return np.abs(po2 - qc["baseline"])


def variant_norm(po2, qc):
    """pO2 / baseline (the point's full-run mean), so 1 = the point's average level."""
    return po2 / qc["baseline"]


# name -> fn(po2, qc). The order here is the column order in `po2_signals`.
PO2_VARIANTS = {
    "po2": variant_raw,
    "abs_po2": variant_abs,
    "norm_po2": variant_norm,
}


def po2_variants(po2, qc):
    """Every PO2_VARIANTS trace: {name: (n_time, n_points)}."""
    return {name: fn(po2, qc) for name, fn in PO2_VARIANTS.items()}


# --- per-event x point metrics (one scalar per `po2_events` row) ----------------------------------------------

def normalize_event_segments(segments, rms):
    """
    (segments - rms) / rms: a point's event windows relative to its full-run RMS (-> event_trace_norm).
    Dividing by each point's own RMS puts points with different absolute pO2 on one scale.
    """
    return (segments - rms) / rms


def _abs_peak(t_rel, norm):
    """
    Signed value and time of the largest |norm| in each window. Both are NaN for all-NaN windows.
    Shared by the event_peak* metrics so they always describe the same sample.
    """
    peak = np.full(len(norm), np.nan)
    peak_time = np.full(len(norm), np.nan)
    ok = ~np.all(np.isnan(norm), axis=1)
    idx = np.nanargmax(np.abs(norm[ok]), axis=1)
    peak[ok] = norm[ok][np.arange(len(idx)), idx]
    peak_time[ok] = t_rel[idx]
    return peak, peak_time


def po2_peak(t_rel, norm):
    """Signed largest deviation of (pO2 - rms) / rms within the window (positive = above the point's RMS)."""
    return _abs_peak(t_rel, norm)[0]


def po2_peak_abs_height(t_rel, norm):
    """|event_peak|: size of the largest deviation, whatever its sign."""
    return np.abs(_abs_peak(t_rel, norm)[0])


def po2_peak_time(t_rel, norm):
    """Time (s, relative to the event) of event_peak. NOTE (known issue A): up to one pO2 sample early."""
    return _abs_peak(t_rel, norm)[1]


def po2_abs_peak_time(t_rel, norm):
    """|event_peak_time|: how far (s) from the event the peak falls, before or after."""
    return np.abs(_abs_peak(t_rel, norm)[1])


def po2_rms_over_baseline(t_rel, norm):
    """
    rms(event window) / rms(whole-run pO2) of the point. Since norm = (window - rms) / rms, the window over the
    run RMS is norm + 1, so this is sqrt(mean((norm + 1)^2)). 1 = as large as the run average. NaN if all-NaN.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN windows
        return np.sqrt(np.nanmean((norm + 1) ** 2, axis=1))


# name -> fn(t_rel (n_window,), norm (n_events, n_window)) -> (n_events,).
# The order here is the column order in `po2_events`. norm = normalize_event_segments(...) of one point.
PO2_EVENT_METRICS = {
    "event_peak": po2_peak,
    "peak_abs_height": po2_peak_abs_height,
    "event_peak_time": po2_peak_time,
    "abs_event_peak_time": po2_abs_peak_time,
    "rms_over_baseline": po2_rms_over_baseline,
}


def po2_event_metrics(t_rel, norm):
    """Every PO2_EVENT_METRICS metric of one point's events: {name: (n_events,)}."""
    return {name: fn(t_rel, norm) for name, fn in PO2_EVENT_METRICS.items()}


# =============================================================================
# 9. Per-run drivers: stages -> tables
# =============================================================================

@dataclass
class LoadedRun:
    """
    Output of stages 1-2 for one run (file I/O + filtering). It is cached by the GUI, so a threshold change
    only re-runs `analyze_run`.

    keys : RUN_KEYS values;  sig, po2 : load_synced_run output;  ext, fs : preprocess_event_source output
    event_method : the method `ext` was computed for
    """
    keys: dict
    sig: dict
    po2: np.ndarray
    ext: dict
    fs: float
    event_method: str


def load_run(folder, run, params):
    """STAGES 1-2 for one run -> LoadedRun. Warns if more than half of the pO2 samples are NaN."""
    sig, po2 = load_synced_run(folder, run, EVENT_METHODS[params.event_method].raw)
    nan_frac = np.isnan(po2).mean()
    if nan_frac > 0.5:
        warnings.warn(f"{folder} run {run}: {nan_frac:.0%} of pO2 samples are NaN")
    ext, fs = preprocess_event_source(sig, params)
    return LoadedRun(keys=run_keys(folder, run), sig=sig, po2=po2, ext=ext, fs=fs, event_method=params.event_method)


def _with_keys(df, keys):
    """Prepend constant key columns (e.g. the run identity) to a DataFrame, in place. Returns df."""
    for i, (name, value) in enumerate(keys.items()):
        df.insert(i, name, value)
    return df


def _po2_signals_table(loaded):
    """`po2_signals`: one row per pO2 point with its PO2_VARIANTS traces and POINT_METRICS scalars."""
    qc = po2_point_qc(loaded.po2)
    variants = po2_variants(loaded.po2, qc)
    return pd.DataFrame([{
        **loaded.keys,
        "point": point,
        "t_po2": loaded.sig["t_po2"],
        **{name: values[:, point] for name, values in variants.items()},
        **{name: values[point] for name, values in qc.items()},
    } for point in range(loaded.po2.shape[1])])


def _events_table(loaded, ev, params):
    """
    `events`: one row per event with its time, lobe, EVENT_METRICS and its window of every method signal.
    Built column-wise so that a run without events still has every column.
    """
    t_events = ev["event_time"]
    n_events = len(t_events)
    events = {"event": np.arange(n_events), "event_time": t_events, "event_bounds": list(ev["event_bounds"]),
              **{name: ev[name] for name in EVENT_METRICS}}
    for signal in EVENT_METHODS[params.event_method].signals:
        t_rel, segments = extract_event_segments(loaded.ext["t_" + signal], loaded.ext[signal], t_events,
                                                 params.extract_window_s)
        events[f"t_{signal}_event"] = [t_rel] * n_events
        events[f"{signal}_event_trace"] = list(segments)
    return _with_keys(pd.DataFrame(events), loaded.keys)


def _po2_events_table(loaded, po2_signals, t_events, params):
    """
    `po2_events`: one row per event x point with the raw and RMS-normalized pO2 window and PO2_EVENT_METRICS.
    """
    n_events = len(t_events)
    parts = []
    for row in po2_signals.itertuples():
        t_rel, segments = extract_event_segments(row.t_po2, row.po2, t_events, params.extract_window_s)
        norm = normalize_event_segments(segments, row.rms)
        parts.append(pd.DataFrame({
            "event": np.arange(n_events),
            "point": row.point,
            "t_event": [t_rel] * n_events,
            "event_time": t_events,
            "event_trace": list(segments),
            "event_trace_norm": list(norm),
            **po2_event_metrics(t_rel, norm),
        }))
    return _with_keys(pd.concat(parts, ignore_index=True), loaded.keys)


def analyze_run(loaded, params):
    """
    STAGES 3-5 for one loaded run -> TransientData with that run's five tables.

    `params.event_method` must match `loaded.event_method`, because the filtered signals depend on it.
    The other fields (threshold, spacing, window) can change freely; the GUI relies on this.
    """
    if loaded.event_method != params.event_method:
        raise ValueError(f"run was preprocessed for {loaded.event_method!r}, params ask for {params.event_method!r}")
    method = EVENT_METHODS[params.event_method]

    dataset = pd.DataFrame([{**loaded.keys, **loaded.sig}])
    po2_signals = _po2_signals_table(loaded)

    # stage 3: detect on the method's detection signal
    detect_on = method.detect_on
    ev, threshold = detect_events(loaded.ext["t_" + detect_on], loaded.ext[detect_on], loaded.fs,
                                  params.height_thresh, params.min_distance_s, rel_to_rms=params.rel_to_rms, neg_peaks=params.neg_peaks)
    extract_signals = pd.DataFrame([{**loaded.keys, **loaded.ext, "event_threshold": threshold}])

    # stages 4-5: event windows of the method signals and of every pO2 point
    events = _events_table(loaded, ev, params)
    po2_events = _po2_events_table(loaded, po2_signals, ev["event_time"], params)
    return TransientData(dataset=dataset, po2_signals=po2_signals, extract_signals=extract_signals,
                         events=events, po2_events=po2_events)


def process_run(folder, run, params):
    """STAGES 1-5 for one run (load_run + analyze_run) -> TransientData."""
    return analyze_run(load_run(folder, run, params), params)


def process_all_runs(runs, params, verbose=True):
    """Apply process_run to every {session folder: [run numbers]} entry and stack the tables."""
    parts = []
    for folder, run_list in runs.items():
        if verbose:
            date, mouse_id = parse_session(folder)
            print(f"Processing mouse {mouse_id} on date {date}: runs {run_list}")
        parts.extend(process_run(folder, run, params) for run in run_list)
    return TransientData.concat(parts)


# =============================================================================
# 10. Event-locked traces and averages (one method: interpolate onto a shared dt grid)
# =============================================================================
# Runs differ in sampling rate (pO2: 0.33-0.92 Hz), so the windows of different runs have different sample
# times. Every average therefore interpolates each window onto np.arange(w0, w1 + dt/2, dt) first
# (interpolate_event). Then it works with a long table (event_responses) and groups it (group_average).

def signal_columns(signal):
    """
    (table, time column, trace column) holding one signal's event windows.

    "po2" / "po2_norm" are the raw and (trace - rms) / rms pO2 windows in `po2_events`. Any other signal s
    (e.g. "pupil", "dxdt_pupil", "locom") is the t_s_event / s_event_trace pair in `events`.
    """
    if signal == "po2":
        return "po2_events", "t_event", "event_trace"
    if signal == "po2_norm":
        return "po2_events", "t_event", "event_trace_norm"
    return "events", f"t_{signal}_event", f"{signal}_event_trace"


def _keep_points(df, points):
    """Rows of a point-level table whose (run, point) is in `points` (a `point_qc` subset); all rows if None."""
    if points is None:
        return df
    return df.merge(points[POINT_KEYS], on=POINT_KEYS)


def event_traces(data, signal, points=None):
    """
    Yield (row, t_rel, trace) for every event of one signal. trace is (n_window,) for event-source signals and
    (n_window, n_points) for pO2, whose per-point `po2_events` rows are stacked back into one event.

    Args:
        points : optional `point_qc` subset (e.g. QC-passing points). Only those points are kept from pO2
                 windows, and events with none left are skipped. Ignored for other signals.
    """
    table, t_col, trace_col = signal_columns(signal)
    df = getattr(data, table)
    if table != "po2_events":
        for row in df.itertuples():
            yield row, getattr(row, t_col), getattr(row, trace_col)
        return

    df = _keep_points(df, points)
    for _, g in df.groupby(EVENT_KEYS, sort=False):
        trace = np.column_stack(g[trace_col].to_list())
        row = next(g.itertuples())
        yield row, getattr(row, t_col), trace


def interpolate_event(t_rel, trace, extract_window_s, dt):
    """
    Resample one event window onto the regular grid np.arange(w0, w1 + dt/2, dt), so runs with different
    sampling rates can be averaged sample by sample.

    2-D traces (n_window, n_points) are averaged over points first. Grid times outside the finite samples
    (e.g. the NaN padding at recording edges) are NaN: values are never extrapolated.
    Returns (t_grid, values), both (n_grid,).
    """
    t_grid = np.arange(extract_window_s[0], extract_window_s[1] + dt / 2, dt)
    if trace.ndim == 2:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)  # samples where every point is NaN
            trace = np.nanmean(trace, axis=1)
    values = np.full(len(t_grid), np.nan)
    ok = np.isfinite(trace)
    if ok.any():
        t_ok = t_rel[ok]
        inside = (t_grid >= t_ok[0]) & (t_grid <= t_ok[-1])
        values[inside] = np.interp(t_grid[inside], t_ok, trace[ok])
    return t_grid, values


def event_responses(data, signal, extract_window_s, dt=1.0, level="event", points=None):
    """
    Every event window of `signal` interpolated onto the shared dt grid, as one long table.

    Args:
        signal : an event-source signal (e.g. "pupil", "locom") or "po2" / "po2_norm" (see signal_columns)
        level  : "event" -> one trace per event (pO2 averaged over the event's points first, so each event counts
                            once however many points its run has);
                 "point" -> one trace per event x pO2 point (`po2_events` row); pO2 signals only
        points : optional `point_qc` subset; see `event_traces`
    Returns DataFrame: EVENT_KEYS (+ point if level="point") + [t, value]
        -> group_average, and the GUI bin averages
    """
    parts = []
    if level == "event":
        key_cols = EVENT_KEYS
        for row, t_rel, trace in event_traces(data, signal, points):
            t_grid, values = interpolate_event(t_rel, trace, extract_window_s, dt)
            parts.append(_with_keys(pd.DataFrame({"t": t_grid, "value": values}),
                                    {k: getattr(row, k) for k in key_cols}))
    elif level == "point":
        table, t_col, trace_col = signal_columns(signal)
        if table != "po2_events":
            raise ValueError(f"level='point' needs a pO2 signal ('po2' or 'po2_norm'), got {signal!r}")
        key_cols = EVENT_POINT_KEYS
        for row in _keep_points(data.po2_events, points).itertuples():
            t_grid, values = interpolate_event(getattr(row, t_col), getattr(row, trace_col), extract_window_s, dt)
            parts.append(_with_keys(pd.DataFrame({"t": t_grid, "value": values}),
                                    {k: getattr(row, k) for k in key_cols}))
    else:
        raise ValueError(f"level must be 'event' or 'point', got {level!r}")
    if not parts:
        return pd.DataFrame(columns=key_cols + ["t", "value"])
    return pd.concat(parts, ignore_index=True)


def group_average(responses, labels=None, by=None):
    """
    Mean, std and count of `event_responses` output at each grid time, per group.

    Args:
        responses : long table with key columns + [t, value] (from event_responses)
        labels    : optional DataFrame with key columns + the `by` columns that `responses` lacks (e.g.
                    data.events[EVENT_KEYS + ["event_class"]], or pt.bin_labels output). It is merged with
                    `responses` on every column they share, so labels can be per event while responses are per
                    event x point. Not needed when every `by` column is already in `responses` (e.g. "mouse_id",
                    "point"). Rows whose group is NaN (e.g. events in no bin) are dropped.
        by        : None (one overall average), a column name, or a list of columns
    Returns DataFrame: [by...] + [t, mean, std, n]   (n = rows with a finite value at that t)
    """
    by = [by] if isinstance(by, str) else list(by or [])
    if labels is not None:
        on = [c for c in labels.columns if c in responses.columns]
        responses = responses.merge(labels[on + [c for c in by if c not in responses.columns]], on=on)
    if by:
        responses = responses.dropna(subset=by)
    return (responses.groupby(by + ["t"], observed=True)["value"]
            .agg(mean="mean", std="std", n="count")
            .reset_index())


def trace_mean_std(traces):
    """
    Mean, standard deviation (ddof=1) and count across equal-length traces at each sample, ignoring NaNs.
    For wide arrays (e.g. one run's windows in the GUI); long tables use group_average.

    Args:
        traces : sequence of (n,) arrays
    Returns (mean, std, n), each (n,); NaN where fewer than 1 (mean) or 2 (std) traces have data
    """
    traces = np.vstack(list(traces))
    n = np.sum(np.isfinite(traces), axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # samples with no or one finite trace
        mean = np.nanmean(traces, axis=0)
        std = np.nanstd(traces, axis=0, ddof=1)
    return mean, std, n


def po2_event_traces(data, variant="po2"):
    """
    One pO2 window per `po2_events` row in a chosen normalization (at the run's own sampling rate).

    Args:
        variant : "po2" (raw event_trace), "po2_norm" ((trace - rms) / rms = event_trace_norm), or a
                  PO2_VARIANTS name applied to the window (e.g. "norm_po2" = trace / baseline, "abs_po2")
    Returns DataFrame: EVENT_POINT_KEYS + [t_event, trace]   (same rows and order as data.po2_events)
        -> the GUI's single-run event panel
    """
    pe = data.po2_events
    df = pe[EVENT_POINT_KEYS + ["t_event"]].copy()
    if variant == "po2":
        df["trace"] = pe["event_trace"].to_list()
    elif variant == "po2_norm":
        df["trace"] = pe["event_trace_norm"].to_list()
    elif variant in PO2_VARIANTS:
        # the variant functions take (po2 (n, n_points), qc); apply them to one window with that point's metrics
        qc_cols = list(POINT_METRICS)
        qc = pe[POINT_KEYS].merge(data.po2_signals[POINT_KEYS + qc_cols], on=POINT_KEYS, how="left")
        df["trace"] = [PO2_VARIANTS[variant](trace[:, None], {c: np.array([q[c]]) for c in qc_cols})[:, 0]
                       for trace, (_, q) in zip(pe["event_trace"], qc.iterrows())]
    else:
        raise ValueError(f"Unknown pO2 variant {variant!r}")
    return df


# =============================================================================
# 11. Binning: classify rows by a metric, and line up metrics of different tables
# =============================================================================

# metrics that must not appear as bin / bar / correlation choices
NOT_BINNABLE = {"nan_frac"}

# key columns that rows can also be binned by (as numbers) -> the table they identify rows of.
# "point" bins put every po2_events row of a pO2 point into that point's bin, like the po2_signals metrics.
BIN_KEYS = {"point": "po2_signals"}

# every scalar metric that rows can be binned by -> the TransientData table holding it.
# Derived from the registries, so a new metric becomes binnable without editing this.
BIN_METRICS = {
    **{m: "events" for m in EVENT_METRICS},
    **{m: "po2_signals" for m in POINT_METRICS if m not in NOT_BINNABLE},
    **{m: "po2_events" for m in PO2_EVENT_METRICS},
    **BIN_KEYS,
}

# NOTE (known issue B): these percentile bins leave P33-34 and P67-68 unbinned. Contiguous bins such as
# [[0, 33.3], [33.3, 66.7], [66.7, 100]] would avoid that.
DEFAULT_BINS = [[0, 33.3], [33.3, 66.7], [66.7, 100]]

CORR_X_KEYS = ["session", "mouse_id", "date", "run", "point"]   # key columns that can be a correlation x axis
NO_BIN, MULTI_BIN = -1, -2   # correlation_values color_idx: no binned row / binned rows in different bins


def bin_edges(values, bins, percentile=True):
    """
    [lo, hi] edges of each bin in units of `values`.

    Args:
        bins       : (n_bins, 2) [lo, hi] bounds (inclusive)
        percentile : the bounds are percentiles of `values` (otherwise they are already values)
    Returns (n_bins, 2) float array
    """
    bins = np.asarray(bins, dtype=float).reshape(-1, 2)
    if not percentile or len(bins) == 0:
        return bins
    values = np.asarray(values, dtype=float)
    if not np.isfinite(values).any():
        return np.full(bins.shape, np.nan)
    return np.nanpercentile(values, np.clip(bins, 0, 100))


def values_to_percentiles(bounds, values):
    """
    Bounds in metric units -> percentiles of `values` (the inverse of np.nanpercentile's linear interpolation).
    The GUI uses it to convert typed bounds when "Percentiles" is toggled.
    """
    values = np.asarray(values, dtype=float)
    values = np.sort(values[np.isfinite(values)])
    bounds = np.asarray(bounds, dtype=float).reshape(-1, 2)
    if len(values) < 2:
        return np.full(bounds.shape, 0.0 if len(values) else np.nan)
    return np.interp(bounds, values, np.linspace(0, 100, len(values)))


def bin_index(values, edges):
    """
    Bin of each value given (n_bins, 2) [lo, hi] `edges` (inclusive), or -1 if the value is in no bin (or NaN).
    Where bins overlap, the LATER bin wins, because each bin overwrites the ones before it.
    """
    values = np.asarray(values, dtype=float)
    idx = np.full(len(values), -1)
    for k, (lo, hi) in enumerate(np.asarray(edges, dtype=float).reshape(-1, 2)):
        idx[(values >= lo) & (values <= hi)] = k
    return idx


def classify_events(df, column, bins, percentile=True):
    """
    Assign each row of `df` to a bin of `column` (any BIN_METRICS metric, with df its table).

    Args:
        bins       : (n_bins, 2) array of [lo, hi] bounds (inclusive), e.g. [[0, 33], [34, 67], [68, 100]]
        percentile : the bounds are percentiles of `column` (default) or values of `column`
    Returns:
        ordered categorical Series on df.index, with labels like "event_auc P0-33" (percentiles) or
        "event_auc 0.5-2" (values). Rows in a gap between bins (or with a NaN value) are NaN.
    """
    bins = np.asarray(bins, dtype=float).reshape(-1, 2)
    values = df[column].to_numpy(dtype=float)
    edges = bin_edges(values, bins, percentile)
    prefix = "P" if percentile else ""
    labels = [f"{column} {prefix}{lo:g}-{hi:g}" for lo, hi in bins]
    if len(set(labels)) < len(labels):  # identical bins -> number them so the categories stay unique
        labels = [f"{label} (bin {k + 1})" if labels.count(label) > 1 else label for k, label in enumerate(labels)]

    idx = bin_index(values, edges)
    classes = np.full(len(values), None, dtype=object)
    for k, label in enumerate(labels):
        classes[idx == k] = label

    n_missing = int(pd.isna(classes).sum())
    if n_missing:
        warnings.warn(f"classify_events: {n_missing}/{len(values)} rows fall outside every {column} bin")
    return pd.Series(pd.Categorical(classes, categories=labels, ordered=True), index=df.index)


def metric_table(data, metric):
    """The TransientData table holding a BIN_METRICS metric."""
    return getattr(data, BIN_METRICS[metric])


def metric_frame(data, name):
    """
    Key columns of `name`'s table + its values as "value". `name` is a BIN_METRICS metric (including "point").
    This is the common format that lets correlation_values pair metrics from different tables.
    """
    if name == "point":
        df = data.po2_signals[POINT_KEYS]
        return df.assign(value=df["point"])
    table = BIN_METRICS[name]
    return getattr(data, table)[TABLE_KEYS[table] + [name]].rename(columns={name: "value"})


def _quiet_classify(data, metric, bins, percentile):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)  # rows outside every bin are expected here
        return classify_events(metric_table(data, metric), metric, bins, percentile)


def bin_labels(data, metric, bins, percentile=True):
    """
    The bin of every event or every event x point row, ready to merge with event_responses / group_average.

    Returns (labels, per_point):
        metric in `events`      -> labels = EVENT_KEYS + [bin], per_point False (one bin per event)
        metric in `po2_events`  -> labels = EVENT_POINT_KEYS + [bin], per_point True
        metric in `po2_signals` -> every po2_events row takes its point's bin (EVENT_POINT_KEYS + [bin]), True
    "bin" is an ordered categorical (classify_events labels); NaN means the row is in no bin.
    """
    classes = _quiet_classify(data, metric, bins, percentile)
    table = BIN_METRICS[metric]
    if table == "events":
        return data.events[EVENT_KEYS].assign(bin=classes), False
    point_rows = data.po2_events[EVENT_POINT_KEYS]
    if table == "po2_events":
        return point_rows.assign(bin=classes), True
    run_points = data.po2_signals[POINT_KEYS].assign(bin=classes)
    return point_rows.merge(run_points, on=POINT_KEYS, how="left"), True


def bar_values(data, metric, labels, per_point):
    """
    `metric` (an events or po2_events metric) of every binned row, with its "bin" (from bin_labels).

    Rows follow the metric's table: one per event, or one per event x point. A po2_events metric with event bins
    takes its event's bin. An events metric with point bins keeps one row per event and bin it falls in, not one
    per point.
    Returns DataFrame: key columns + [bin, metric]
    """
    table_name = BIN_METRICS[metric]
    table = getattr(data, table_name)
    on = [c for c in labels.columns if c != "bin" and c in table.columns]
    df = labels.merge(table[on + [metric]], on=on)
    if table_name == "events" and per_point:
        df = df.drop_duplicates(EVENT_KEYS + ["bin"])
    return df[df["bin"].notna()].reset_index(drop=True)


def correlation_values(data, x, y, bin_metric, bins, percentile=True):
    """
    Rows of a y-vs-x scatter, colored by the bins of `bin_metric`.

    x and y are BIN_METRICS metrics (x may also be a CORR_X_KEYS key). They are merged on the key columns their
    tables share; for example rms vs event_auc pairs every event with every point of its run. Run keys other
    than point come from y's rows.
    color_idx is the bin of the binned rows each scatter row matches (merged on shared keys): bin k if they all
    fall in bin k, NO_BIN if none is binned, MULTI_BIN if they differ (e.g. per-point rms vs bins of the
    per-event x point event_peak).
    Returns DataFrame: key columns + [y, x, color_idx]
    """
    df = metric_frame(data, y).rename(columns={"value": "y"})
    if x in CORR_X_KEYS and x != "point":
        df["x"] = df[x]
    else:
        x_frame = metric_frame(data, x).rename(columns={"value": "x"})
        df = df.merge(x_frame, on=[c for c in x_frame.columns if c in df.columns])
    df = df[df["x"].notna() & np.isfinite(df["y"].to_numpy(dtype=float))].reset_index(drop=True)

    classes = _quiet_classify(data, bin_metric, bins, percentile)
    bin_table = BIN_METRICS[bin_metric]
    bin_frame = getattr(data, bin_table)[TABLE_KEYS[bin_table]].assign(bin=classes.cat.codes)
    on = [c for c in bin_frame.columns if c != "bin" and c in df.columns]
    matched = df[on].assign(_row=np.arange(len(df))).merge(bin_frame, on=on, how="left")
    per_row = matched.groupby("_row")["bin"]
    n_distinct = per_row.nunique(dropna=False).reindex(range(len(df)), fill_value=0)
    first = per_row.first().reindex(range(len(df)))
    color_idx = np.where(n_distinct > 1, MULTI_BIN, first.fillna(NO_BIN)).astype(int)
    return df.assign(color_idx=color_idx)


# =============================================================================
# 12. Statistics: point QC and linear mixed models
# =============================================================================

def point_qc(data):
    """Run x pO2 point quality table (the scalar columns of `po2_signals`): POINT_KEYS + POINT_METRICS."""
    return data.po2_signals[POINT_KEYS + list(POINT_METRICS)].copy()


def event_point_table(data, points=None):
    """
    The scalar columns of `po2_events` (one row per event x pO2 point), for mixed models.

    run_id, event_id and point_id are unique across sessions (run / event / point numbers repeat), so they can be
    nested random effects. point_id is one pO2 location within one run.

    Args:
        points : optional `point_qc` subset (e.g. QC-passing points); only those points are kept
    Returns DataFrame: EVENT_KEYS + [run_id, event_id, point_id, point] + PO2_EVENT_METRICS
    """
    df = _keep_points(data.po2_events[EVENT_KEYS + ["point"] + list(PO2_EVENT_METRICS)], points)
    df.insert(len(EVENT_KEYS), "run_id", df["date"] + "_run" + df["run"].astype(str))
    df.insert(len(EVENT_KEYS) + 1, "event_id", df["run_id"] + "_event" + df["event"].astype(str))
    df.insert(len(EVENT_KEYS) + 2, "point_id", df["run_id"] + "_point" + df["point"].astype(str))
    return df


def fit_mixed_model(df, outcome, fixed_effects, random_effects, groups, references=None, formula=None):
    """
    Fit a linear mixed model with categorical fixed effects, a random intercept per `groups`, and each of
    `random_effects` as a variance component nested within `groups`.

    Args:
        references : optional {fixed effect: reference level} for treatment coding
        formula    : optional full formula (overrides the one built from fixed_effects)
    Returns:
        (fit result, rows used in the fit, formula)
    """
    references = references or {}
    cols = list(dict.fromkeys([outcome] + list(fixed_effects) + list(random_effects) + [groups]))
    model_df = df[cols].copy()
    model_df[outcome] = pd.to_numeric(model_df[outcome], errors="coerce")
    model_df = model_df.dropna()
    model_df = model_df[np.isfinite(model_df[outcome])].reset_index(drop=True)
    print(f"Fitting on {len(model_df)}/{len(df)} rows with a finite {outcome}")

    if formula is None:
        terms = [f'C({fe}, Treatment(reference="{references[fe]}"))' if fe in references else f"C({fe})"
                 for fe in fixed_effects]
        formula = f"{outcome} ~ " + (" + ".join(terms) or "1")

    model = smf.mixedlm(formula, data=model_df, groups=model_df[groups], re_formula="1",
                        vc_formula={re: f"0 + C({re})" for re in random_effects})
    fit = model.fit(reml=True, method="lbfgs", maxiter=1000)
    return fit, model_df, formula


def extract_fixed_effects(fit, alpha=0.05, drop_intercept=True):
    """Fixed-effect estimates with Wald CIs: [term, estimate, se, pvalue, ci_lower, ci_upper, significant]."""
    z = norm.ppf(1 - alpha / 2)
    fe = pd.Series(fit.fe_params, dtype=float)
    effects = pd.DataFrame({
        "term": fe.index,
        "estimate": fe.values,
        "se": pd.Series(fit.bse_fe, index=fe.index).values,
        "pvalue": pd.Series(fit.pvalues, index=fit.params.index).reindex(fe.index).values,
    })
    effects["ci_lower"] = effects["estimate"] - z * effects["se"]
    effects["ci_upper"] = effects["estimate"] + z * effects["se"]
    effects["significant"] = effects["pvalue"] < alpha
    if drop_intercept:
        effects = effects[effects["term"] != "Intercept"]
    return effects.reset_index(drop=True)


def extract_variance_decomposition(fit, model_df, formula, group_name="group"):
    """
    Split the model variance into fixed-effect terms, random effects and the residual.

    Fixed terms are sequential (Type-I): each term's variance is how much the variance of the fixed-effect
    prediction grows when that term is added, so it depends on the term order in `formula`.
    Random effects are the group random-intercept variance and each variance component.
    Returns DataFrame: [component, component_type, variance, percent_total]
    """
    x = dmatrix(formula.split("~", 1)[1].strip(), model_df, return_type="dataframe")
    term_slices = x.design_info.term_name_slices
    x = x.reindex(columns=fit.fe_params.index, fill_value=0.0)

    rows, included, explained = [], [], 0.0
    for term, cols in term_slices.items():
        if term == "Intercept":
            continue
        included += [c for c in x.columns[cols] if c in fit.fe_params.index]
        partial = ["Intercept"] + included
        var = float(np.var(x[partial].to_numpy() @ fit.fe_params[partial].to_numpy(), ddof=1))
        rows.append({"component": term, "component_type": "fixed", "variance": max(var - explained, 0.0)})
        explained = max(var, explained)

    group_var = float(np.diag(fit.cov_re).sum()) if fit.cov_re is not None and fit.cov_re.size else 0.0
    rows.append({"component": group_name, "component_type": "random", "variance": max(group_var, 0.0)})
    for name, var in zip(getattr(fit.model.exog_vc, "names", []), np.asarray(fit.vcomp, dtype=float)):
        rows.append({"component": str(name), "component_type": "random", "variance": max(float(var), 0.0)})
    rows.append({"component": "Residual", "component_type": "residual", "variance": max(float(fit.scale), 0.0)})

    variance = pd.DataFrame(rows)
    total = variance["variance"].sum()
    variance["percent_total"] = 100 * variance["variance"] / total if total > 0 else np.nan
    return variance
