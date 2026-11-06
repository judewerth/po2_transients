"""
presentation_figures.py: publication-style figures for the data presentation.

One shared style for every figure (fonts, palette, sizes, spines, legend, export), so the deck looks
consistent. Every figure is drawn from the same example run (EXAMPLE, same as notebooks/01_pipeline_walkthrough.ipynb) with the same
parameters as the analysis notebooks (PARAMS). Data and filtering come from po2_transients (pt); nothing is
recomputed here.

Run `python presentation_figures.py` to (re)build every figure into figures/ (or PO2_FIG_DIR) as PNG (300 dpi) and SVG.

Sections
--------
1. Style      palette, rcParams, title / legend / save helpers
2. Data       example run + parameters
3. Figures    one `fig_*` function per figure (returns the Figure)
"""

import os
import warnings

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import po2_transients as pt


# =============================================================================
# 1. Style
# =============================================================================

# Validated categorical order (blue, orange, aqua, ...). Assign in this order, never cycle.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]

COLORS = {
    "raw": "#b5b4ae",          # unprocessed signal: recessive neutral, sits behind the processed one
    "filt": SERIES[0],         # processed / primary signal (filtered pupil)
    "dxdt": SERIES[1],         # detection signal (dPupil/dt)
    "event": "#0b0b0b",        # detected events (markers, labels)
    "threshold": "#52514e",    # detection threshold (dashed)
    "window": "#d9d7cf",       # extraction-window shading (context, not data)
    "trace": "#b5b4ae",        # individual traces (e.g. one pO2 point) behind their mean
    "mean": "#0b0b0b",         # average of the individual traces
    "text": "#0b0b0b",
    "text_muted": "#52514e",
    "axis": "#52514e",
    "grid": "#e6e5e0",
    "surface": "#ffffff",
}

# Ordered bins (low -> high metric) are a magnitude, so they get one sequential hue, light -> dark.
# bin_colors(n) picks n evenly spaced steps of this ramp; the same bin is the same color in every panel.
BIN_RAMP = ["mediumorchid", "limegreen", "tomato"]


def bin_colors(n_bins):
    """n_bins colors from BIN_RAMP, lightest = lowest bin."""
    if n_bins == 1:
        return [BIN_RAMP[2]]
    return [BIN_RAMP[round(k * (len(BIN_RAMP) - 1) / (n_bins - 1))] for k in range(n_bins)]


# Axis labels shared by every figure (signal names as in pt.EVENT_METHODS / pt.signal_columns)
SIGNAL_LABELS = {
    "pupil": "Pupil diameter (norm.)",
    "dxdt_pupil": "dPupil/dt (per sample)",
    "locom": "Locomotion (a.u.)",
    "sleep_state": "Sleep state",
    "dxdt_sleep_state": "dSleepState/dt",
    "po2": "pO$_2$ (mmHg)",
    "po2_norm": "Normalized pO$_2$",
}
SIGNAL_NAMES = {"pupil": "Pupil diameter", "dxdt_pupil": "dPupil/dt", "locom": "Locomotion",
                "sleep_state": "Sleep state", "dxdt_sleep_state": "dSleepState/dt",
                "po2": "pO$_2$", "po2_norm": "Normalized pO$_2$"}
SIGNAL_UNITS = {"pupil": "norm.", "dxdt_pupil": "per sample", "locom": "a.u."}
# coded signals: value -> name for the y ticks (sleep state as stored by pt.preprocess_sleep_state: REM 3 -> 4)
SIGNAL_TICKS = {"sleep_state": {0: "Unknown", 1: "Awake", 2: "NREM", 4: "REM"}}

# What an "event" is for each pt.EVENT_METHODS entry, so labels name it in context
# (e.g. "dilation" / "Time from peak dilation rate" instead of "event" / "Time relative to event").
EVENT_CONTEXT = {
    "pupil": {"event": "dilation", "events": "Events", "time": "Time (s) from pupil event"},
    "locom": {"event": "event", "events": "events", "time": "Time from peak locomotion (s)"},
    "sleep_state": {"event": "transition", "events": "transitions", "time": "Time from state transition (s)"},
}

# Event metrics (pt.EVENT_METRICS) are measured on the detection signal: "{sig}" is its name.
METRIC_TEMPLATES = {"event_auc": "{sig} area under the curve", "event_height": "{sig} peak height",
                    "event_halfwidth": "{sig} peak half-width"}
# metric units: (metric, detection signal) first, then metric alone, else the detection signal's own units
METRIC_UNITS = {("event_auc", "dxdt_pupil"): "Δ pupil diameter, norm.", "event_halfwidth": "s"}


def upper_first(text):
    """Uppercase only the first character (str.capitalize also lowercases the rest). For plain words only:
    signal names such as 'dPupil/dt' are kept as written."""
    return text[:1].upper() + text[1:]


def event_context(params):
    """EVENT_CONTEXT entry of params.event_method (generic 'event' wording if it has none)."""
    return EVENT_CONTEXT.get(params.event_method,
                             {"event": "event", "events": "events", "time": "Time relative to event (s)"})


def metric_label(metric, params, units=True):
    """e.g. 'dPupil/dt area under the curve (Δ pupil diameter, norm.)' for event_auc of pupil events."""
    detect_on = pt.EVENT_METHODS[params.event_method].detect_on
    sig = SIGNAL_NAMES.get(detect_on, detect_on)
    name = METRIC_TEMPLATES.get(metric, f"{sig} {metric}").format(sig=sig)
    unit = METRIC_UNITS.get((metric, detect_on), METRIC_UNITS.get(metric, SIGNAL_UNITS.get(detect_on)))
    return f"{name} ({unit})" if units and unit else name

FIGSIZE = {
    "wide": (10, 3.4),         # full slide width, one row (time series)
    "wide2": (10, 6.2),        # full slide width, two stacked rows sharing the time axis
    "wide3": (10, 9.0),        # full slide width, three stacked rows
    "half": (5, 3.4),          # half slide
    "square": (4.2, 4.2),
}

LW = {"raw": 0.8, "filt": 2.0, "threshold": 1.2, "trace": 1.0, "mean": 2.5, "event_line": 1.0}
MARKER = {"s": 55, "color": COLORS["event"], "edgecolor": COLORS["surface"], "linewidth": 1.5, "zorder": 5}
WINDOW_ALPHA = 0.45
FIG_DIR = os.path.join(pt.REPO_ROOT, os.path.expanduser(os.environ.get("PO2_FIG_DIR") or "figures"))  # set PO2_FIG_DIR in .env to change


def apply_style():
    """Set the shared matplotlib style. Called on import, so every figure gets it."""
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.labelsize": 11,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "text.color": COLORS["text"],
        "axes.labelcolor": COLORS["text"],
        "axes.edgecolor": COLORS["axis"],
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.grid.axis": "y",
        "axes.axisbelow": True,
        "grid.color": COLORS["grid"],
        "grid.linewidth": 0.8,
        "xtick.color": COLORS["axis"],
        "ytick.color": COLORS["axis"],
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "legend.frameon": False,
        "lines.solid_capstyle": "round",
        "figure.facecolor": COLORS["surface"],
        "axes.facecolor": COLORS["surface"],
        "figure.dpi": 110,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
        "svg.fonttype": "none",       # keep text editable in SVG
        "pdf.fonttype": 42,
    })


apply_style()


def title(ax, text, subtitle=None):
    """Left-aligned bold title with an optional muted subtitle line underneath."""
    ax.set_title(text, loc="left", fontweight="bold", pad=22 if subtitle else 8)
    if subtitle:
        ax.text(0, 1.02, subtitle, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=10, color=COLORS["text_muted"])


def legend(ax, **kw):
    """Frameless legend, top right above the plot area, in one row."""
    kw = {"loc": "lower right", "bbox_to_anchor": (1, 1.0), "ncol": 4, "handlelength": 1.6,
          "borderaxespad": 0.2, **kw}
    return ax.legend(**kw)


def save(fig, name):
    """Write figures/<name>.png (300 dpi) and figures/<name>.svg."""
    os.makedirs(FIG_DIR, exist_ok=True)
    paths = [os.path.join(FIG_DIR, f"{name}.{ext}") for ext in ("png", "svg")]
    for p in paths:
        fig.savefig(p)
    return paths


# =============================================================================
# 2. Data: example run + parameters (match notebooks/01_pipeline_walkthrough.ipynb / 03_analyze_pupil_po2_transients.ipynb)
# =============================================================================

EXAMPLE = ("26-10-05_B6_318", 6)   # (session folder, run) used for every single-run figure
ALL_RUNS = {
    "26-10-05_B6_318": [5, 6, 8, 9, 11, 12, 14],
    "26-10-07_B6_322": [6, 12],
    "26-10-21_B6_322": [2, 3],
}
# One parameter set per event source. Pupil figures default to PUPIL_PARAMS (PARAMS); pass LOCOM_PARAMS for locomotion.
PUPIL_PARAMS = pt.EventParams(
    event_method="pupil",
    rel_to_rms=False,
    movmedian_len_s=6,           # s
    lpf_cutoff_hz=0.1,           # Hz
    height_thresh=0.004,
    min_distance_s=30,           # s
    extract_window_s=(-30, 30),  # s
)
LOCOM_PARAMS = pt.EventParams(
    event_method="locom",
    rel_to_rms=False,
    movmedian_len_s=6,           # s
    lpf_cutoff_hz=0.1,           # Hz
    height_thresh=0.0025,
    min_distance_s=30,           # s
    extract_window_s=(-30, 30),  # s
)
SLEEP_STATE_PARAMS = pt.EventParams(
    event_method="sleep_state",
    rel_to_rms=False,      
    height_thresh=0.1,
    min_distance_s=10,           # s
    extract_window_s=(-30, 30),  # s
    neg_peaks=True
)
PARAMS = PUPIL_PARAMS            # default of every pupil figure


def example_label(example=EXAMPLE):
    folder, run = example
    date, mouse = pt.parse_session(folder)
    return f"Mouse {mouse} · {date} · run {run:02d}"


_RUN_CACHE = {}


def load_example_run(example=EXAMPLE, params=PARAMS):
    """Stages 1-5 of the pipeline (pt.process_run) for one run -> one-run TransientData. Cached per params."""
    key = (example, repr(params))
    if key not in _RUN_CACHE:
        _RUN_CACHE[key] = pt.process_run(*example, params)
    return _RUN_CACHE[key]


_ALL_CACHE = {}

DT_S = 1.0                                   # s, shared grid that pooled event windows are interpolated onto
TERTILE_BINS = [[0, 33], [33, 66], [66, 100]]   # percentile bins; contiguous, so no event falls between bins
AUC_BINS = TERTILE_BINS
# dSleepState/dt peak height (value bins): -1 = NREM -> awake, +1 = awake -> NREM (awake 1, NREM 2, REM 4)
SLEEP_TRANSITION_BINS = [[-1, 0], [0, 1]]
SLEEP_TRANSITION_NAMES = ["NREM → awake", "Awake → NREM"]


def load_all_runs(runs=ALL_RUNS, params=PARAMS):
    """Stages 1-5 for every run (pt.process_all_runs) -> stacked TransientData. Cached per (runs, params)."""
    key = (repr(runs), repr(params))
    if key not in _ALL_CACHE:
        _ALL_CACHE[key] = pt.process_all_runs(runs, params, verbose=False)
    return _ALL_CACHE[key]


def runs_label(data, events="events"):
    """e.g. '3 sessions · 2 mice · 11 runs · 42 dilations'."""
    ds = data.dataset
    return (f"{ds['session'].nunique()} sessions · {ds['mouse_id'].nunique()} mice · {len(ds)} runs · "
            f"{len(data.events)} {events}")


def run_series(run, name):
    """(t, values) of a signal of a one-run TransientData: raw_<src> from `dataset`, else `extract_signals`."""
    if name.startswith("raw_"):
        row, name = run.dataset.iloc[0], name[len("raw_"):]
    else:
        row = run.extract_signals.iloc[0]
    return row["t_" + name], row[name]


# =============================================================================
# 3. Shared panels (draw onto an existing axis, so every figure shows a signal the same way)
# =============================================================================

def filter_label(params=PARAMS):
    return f"{params.movmedian_len_s:g} s median + {params.lpf_cutoff_hz:g} Hz low-pass"


def panel_pupil(ax, run, params=PARAMS, detailed_label=True):
    """Raw pupil (thin gray) under the filtered pupil (blue), on the shared 0-1 axis."""
    t_raw, raw = run_series(run, "raw_pupil")
    t, filt = run_series(run, "pupil")
    ax.plot(t_raw, raw, color=COLORS["raw"], lw=LW["raw"], label="Raw")
    ax.plot(t, filt, color=COLORS["filt"], lw=LW["filt"],
            label=f"Filtered ({filter_label(params)})" if detailed_label else "Filtered")
    ax.set_xlim(0, t[-1])
    ax.set_ylim(-0.05, 1.05)          # filtfilt can overshoot the 0-1 range slightly; don't clip it
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_ylabel(SIGNAL_LABELS["pupil"])


def panel_event_windows(ax, run, params=PARAMS, label="Extraction window"):
    """Shade each event's extraction window [t + w0, t + w1]. Only the first span carries the legend label."""
    w0, w1 = params.extract_window_s
    for i, t_event in enumerate(run.events["event_time"]):
        ax.axvspan(t_event + w0, t_event + w1, color=COLORS["window"], alpha=WINDOW_ALPHA, lw=0, zorder=0,
                   label=label if i == 0 else None)


def panel_event_markers(ax, run, signal, label="Event", number=True):
    """Marker at each event on `signal` (interpolated to the event time); optional event number above it."""
    t, x = run_series(run, signal)
    t_ev = run.events["event_time"].to_numpy()
    y_ev = np.interp(t_ev, t, x)
    ax.scatter(t_ev, y_ev, label=label, **MARKER)
    if number:
        for i, (te, ye) in enumerate(zip(t_ev, y_ev), start=1):     # 1-based for the audience
            ax.annotate(str(i), (te, ye), xytext=(0, 9), textcoords="offset points", ha="center",
                        va="bottom", fontsize=10, fontweight="bold", color=COLORS["event"])


def panel_detection(ax, run, params=PARAMS):
    """Detection signal (dPupil/dt) with the threshold and each detected peak."""
    detect_on = pt.EVENT_METHODS[params.event_method].detect_on
    t, x = run_series(run, detect_on)
    threshold = run.extract_signals["event_threshold"].iloc[0]
    ax.axhline(0, color=COLORS["axis"], lw=0.8, zorder=1)
    ax.plot(t, x, color=COLORS["dxdt"], lw=LW["filt"], label=SIGNAL_NAMES[detect_on])
    ax.axhline(threshold, color=COLORS["threshold"], lw=LW["threshold"], ls=(0, (5, 3)), zorder=2,
               label=f"Threshold ({threshold:g})")
    ax.scatter(run.events["event_time"], run.events["event_height"],
               label=f"Detected {event_context(params)['event']}", **MARKER)
    ax.set_ylabel(SIGNAL_LABELS[detect_on])


# =============================================================================
# 4. Figures
# =============================================================================

def fig_pupil_raw_filtered(example=EXAMPLE, params=PARAMS):
    """Raw pupil and the filtered pupil used for event detection, full run."""
    run = load_example_run(example, params)
    fig, ax = plt.subplots(figsize=FIGSIZE["wide"], layout="constrained")
    panel_pupil(ax, run, params)
    ax.set_xlabel("Time (s)")
    title(ax, "Pupil diameter", example_label(example))
    legend(ax)
    return fig


def fig_pupil_event_detection(example=EXAMPLE, params=PARAMS):
    """
    Top: raw + filtered pupil with each event's extraction window and an event marker.
    Bottom: dPupil/dt (detection signal) with the threshold and the detected peaks.
    """
    run = load_example_run(example, params)
    fig, (ax_p, ax_d) = plt.subplots(2, 1, figsize=FIGSIZE["wide2"], sharex=True, layout="constrained",
                                     height_ratios=[1.25, 1])

    ctx = event_context(params)
    panel_event_windows(ax_p, run, params)
    panel_pupil(ax_p, run, params, detailed_label=False)
    panel_event_markers(ax_p, run, "pupil", label=upper_first(ctx["event"]))
    w0, w1 = params.extract_window_s
    title(ax_p, "Detecting pupil dilations from dPupil/dt",
          f"{example_label(example)} · {len(run.events)} {ctx['events']} · window −{abs(w0):g} to +{w1:g} s")
    legend(ax_p)

    panel_event_windows(ax_d, run, params, label=None)
    panel_detection(ax_d, run, params)
    ymax = np.nanmax(np.abs(run_series(run, "dxdt_pupil")[1]))
    ax_d.set_ylim(-1.15 * ymax, 1.15 * ymax)
    ax_d.yaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=4, symmetric=True))
    ax_d.set_xlabel("Time (s)")
    legend(ax_d)
    return fig


def event_row(run, event_number):
    """Row of run.events for a 1-based event number (the numbers shown in fig_pupil_event_detection)."""
    n = len(run.events)
    if not 1 <= event_number <= n:
        raise ValueError(f"event_number must be 1..{n}, got {event_number}")
    return run.events.iloc[event_number - 1]


def fig_event_pupil_po2(event_number=1, example=EXAMPLE, params=PARAMS):
    """
    One event, on time relative to the event (the extraction window):
    Top: filtered pupil (left axis) and dPupil/dt (right axis) with the threshold and the event peak.
    Bottom: raw pO2 of every measurement point (gray) and their mean across points (black).
    Uses the windows the pipeline extracted (events / po2_events tables), so samples outside the recording are NaN.
    """
    run = load_example_run(example, params)
    ev = event_row(run, event_number)
    w0, w1 = params.extract_window_s
    threshold = run.extract_signals["event_threshold"].iloc[0]

    fig, (ax_p, ax_o) = plt.subplots(2, 1, figsize=FIGSIZE["wide2"], sharex=True, layout="constrained",
                                     height_ratios=[1, 1.25])

    # --- top: filtered pupil + dPupil/dt (two units -> two y-axes, each label in its line's color)
    ax_d = ax_p.twinx()
    ax_d.spines["right"].set_visible(True)
    ax_d.grid(False)
    for ax in (ax_p, ax_d):
        ax.axvline(0, color=COLORS["axis"], lw=LW["event_line"], ls=":", zorder=1)
    ax_p.plot(ev["t_pupil_event"], ev["pupil_event_trace"], color=COLORS["filt"], lw=LW["filt"],
              label="Filtered pupil")
    ax_d.axhline(threshold, color=COLORS["threshold"], lw=LW["threshold"], ls=(0, (5, 3)),
                 label=f"Threshold ({threshold:g})")
    ax_d.plot(ev["t_dxdt_pupil_event"], ev["dxdt_pupil_event_trace"], color=COLORS["dxdt"], lw=LW["filt"],
              label="dPupil/dt")
    ax_d.scatter([0], [ev["event_height"]], label="Peak dilation rate", **MARKER)

    ax_p.set_ylim(-0.05, 1.05)
    ax_p.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax_p.set_ylabel(SIGNAL_LABELS["pupil"], color=COLORS["filt"])
    ymax = 1.15 * np.nanmax(np.abs(ev["dxdt_pupil_event_trace"]))
    ax_d.set_ylim(-ymax, ymax)
    ax_d.yaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=4, symmetric=True))
    ax_d.set_ylabel(SIGNAL_LABELS["dxdt_pupil"], color=COLORS["dxdt"])
    ax_p.set_zorder(ax_d.get_zorder() + 1)    # pupil line drawn above the derivative axis
    ax_p.patch.set_visible(False)

    ctx = event_context(params)
    title(ax_p, f"{upper_first(ctx['event'])} {event_number}: pupil and pO$_2$ response",
          f"{example_label(example)} · peak dilation rate at t = {ev['event_time']:.0f} s")
    h1, l1 = ax_p.get_legend_handles_labels()
    h2, l2 = ax_d.get_legend_handles_labels()
    legend(ax_p, handles=h1 + h2, labels=l1 + l2, ncol=2)

    # --- bottom: raw pO2 of every point + mean across points
    po2 = run.po2_events[run.po2_events["event"] == ev["event"]]
    t_rel = po2["t_event"].iloc[0]
    traces = np.vstack(po2["event_trace"].to_numpy())             # (n_points, n_time)
    ax_o.axvline(0, color=COLORS["axis"], lw=LW["event_line"], ls=":", zorder=1)
    for i, trace in enumerate(traces):
        ax_o.plot(t_rel, trace, color=COLORS["trace"], lw=LW["trace"],
                  label=f"Individual pO$_2$ points (n = {len(traces)})" if i == 0 else None)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)          # samples outside the recording: all NaN
        mean = np.nanmean(traces, axis=0)
    ax_o.plot(t_rel, mean, color=COLORS["mean"], lw=LW["mean"], label="Mean of pO$_2$ points")
    ax_o.set_ylim(bottom=0)
    ax_o.set_ylabel(SIGNAL_LABELS["po2"])
    ax_o.set_xlabel(ctx["time"])
    ax_o.set_xlim(w0, w1)
    legend(ax_o)

    # parts of the window outside the recording are NaN: say so instead of leaving them blank
    t_rec = run.dataset["t_pupil"].iloc[0]
    for lo, hi, text in [(w0, t_rec[0] - ev["event_time"], "Before recording start"),
                         (t_rec[-1] - ev["event_time"], w1, "After recording end")]:
        if hi > lo:
            narrow = (hi - lo) < 0.3 * (w1 - w0)          # too narrow for horizontal text: run it vertically
            for ax in (ax_p, ax_o):
                ax.axvspan(lo, hi, color=COLORS["window"], alpha=WINDOW_ALPHA, lw=0, zorder=0)
                ax.text((lo + hi) / 2, 0.5, text, transform=ax.get_xaxis_transform(), ha="center",
                        va="center", fontsize=10 if not narrow else 9, color=COLORS["text_muted"],
                        rotation=90 if narrow else 0)
    return fig


def panel_bin_traces(ax, avg, colors, ylabel):
    """Mean ± std per bin (pt.group_average by "bin"), bin k in colors[k]. Returns the drawn line handles."""
    ax.axvline(0, color=COLORS["axis"], lw=LW["event_line"], ls=":", zorder=1)
    for label, color in zip(avg["bin"].cat.categories, colors):
        g = avg[avg["bin"] == label]
        if g.empty:
            continue
        std = g["std"].fillna(0)
        lo, hi = g["mean"] - std, g["mean"] + std
        ax.fill_between(g["t"], lo, hi, color=color, alpha=0.10, lw=0, zorder=2)
        for edge in (lo, hi):                       # thin edges keep each band readable where bands overlap
            ax.plot(g["t"], edge, color=color, lw=0.7, alpha=0.7, zorder=2)
        ax.plot(g["t"], g["mean"], color=color, lw=LW["filt"], zorder=3)
    ax.set_ylabel(ylabel)


def panel_bin_histogram(ax, values, edges, colors, xlabel, n_hist=20, ylabel="Events"):
    """Stacked histogram of a metric, each bar segment in its bin's color (pt.bin_index); unbinned values gray."""
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    idx = pt.bin_index(values, edges)
    stacks = [values[idx == k] for k in range(len(edges))]
    cols = list(colors)
    if np.any(finite & (idx == -1)):
        stacks.append(values[finite & (idx == -1)])
        cols.append(COLORS["raw"])
    v = values[finite]
    discrete = len(v) and np.all(v == np.round(v)) and len(np.unique(v)) <= n_hist
    if discrete:                       # integer-valued metric (e.g. sleep-state transitions): one bar per value
        hist_edges = np.arange(v.min() - 0.5, v.max() + 1.5)
        ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
    else:
        hist_edges = np.histogram_bin_edges(v, bins=n_hist)
    ax.hist(stacks, bins=hist_edges, stacked=True, color=cols, edgecolor=COLORS["surface"], linewidth=1,
            rwidth=0.6 if discrete else None)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))


def fig_binned_event_responses(signal=None, po2_signal="po2", metric="event_auc", bins=AUC_BINS,
                               percentile=True, runs=ALL_RUNS, params=PARAMS, dt=DT_S, title_text=None,
                               bin_names=None):
    """
    Event-triggered responses pooled over `runs`, grouped by bins of an event metric:
    Top: mean ± std of the behavior `signal` per bin.  Middle: same for `po2_signal` (each event's pO2 is first
    averaged over its points, so every event counts once).  Bottom: histogram of `metric`, colored by bin.

    signal defaults to the filtered event source of params.event_method (e.g. "pupil", "locom"), so the same
    figure works for any behavior: pass params with another event_method.
    bin_names : optional legend name per bin (e.g. ["NREM → awake", "Awake → NREM"]); default shows the bounds.
    """
    method = pt.EVENT_METHODS[params.event_method]
    signal = signal or method.signals[0]
    data = load_all_runs(runs, params)

    # bins (pt): one label per event, then mean ± std per bin of each signal's interpolated windows
    labels, per_point = pt.bin_labels(data, metric, bins, percentile)
    if per_point:
        raise ValueError(f"{metric!r} is not an event metric; this figure bins events")
    values = pt.metric_table(data, metric)[metric].to_numpy(dtype=float)
    edges = pt.bin_edges(values, bins, percentile)
    avgs = {s: pt.group_average(pt.event_responses(data, s, params.extract_window_s, dt, level="event"),
                                labels, "bin")
            for s in (signal, po2_signal)}
    colors = bin_colors(len(edges))
    counts = labels["bin"].value_counts(sort=False)

    fig = plt.figure(figsize=FIGSIZE["wide3"], layout="constrained")
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 0.85])
    ax_b = fig.add_subplot(gs[0])
    ax_o = fig.add_subplot(gs[1], sharex=ax_b)
    ax_h = fig.add_subplot(gs[2])

    ctx = event_context(params)
    sig_name, po2_name = SIGNAL_NAMES.get(signal, signal), SIGNAL_NAMES.get(po2_signal, po2_signal)
    metric_name = metric_label(metric, params, units=False)

    panel_bin_traces(ax_b, avgs[signal], colors, SIGNAL_LABELS.get(signal, signal))
    panel_bin_traces(ax_o, avgs[po2_signal], colors, SIGNAL_LABELS.get(po2_signal, po2_signal))
    if signal in SIGNAL_TICKS:         # coded signal (e.g. sleep state): name the codes inside the y range
        lo, hi = ax_b.get_ylim()
        ticks = {v: n for v, n in SIGNAL_TICKS[signal].items() if lo <= v <= hi}
        ax_b.set_yticks(list(ticks), list(ticks.values()))
    ax_b.tick_params(labelbottom=False)
    ax_o.set_xlim(*params.extract_window_s)
    ax_o.set_xlabel(ctx["time"])
    panel_bin_histogram(ax_h, values, edges, colors, metric_label(metric, params), ylabel=upper_first(ctx["events"]))

    # one legend for the whole figure: bin k = color k in every panel
    prefix = "P" if percentile else ""
    names = bin_names or [f"{prefix}{b[0]:g}–{b[1]:g}" for b in np.asarray(bins).reshape(-1, 2)]
    handles = [mpl.patches.Patch(color=c, label=f"{name} (n = {counts.iloc[k]})")
               for k, (c, name) in enumerate(zip(colors, names))]
    legend(ax_b, handles=handles, ncol=len(handles), handlelength=1.0)
    title(ax_b, title_text or f"{sig_name} and {po2_name} around {ctx['events']}, by {metric_name}",
          f"{runs_label(data, ctx['events'])} · mean ± SD")
    ax_o.set_title(f"{po2_name} (each {ctx['event']} averaged over pO$_2$ points)", loc="left",
                   fontweight="bold")
    ax_h.set_title(f"{metric_name} "
                   f"per {ctx['event']}, colored by bin", loc="left", fontweight="bold")
    return fig


FIGURES = {
    "pupil_raw_filtered": fig_pupil_raw_filtered,
    "pupil_event_detection": fig_pupil_event_detection,
    "event1_pupil_po2": lambda: fig_event_pupil_po2(event_number=1),
    "pupil_po2_by_auc": fig_binned_event_responses,
    "locom_po2_by_height": lambda: fig_binned_event_responses(
        params=LOCOM_PARAMS, metric="event_height", bins=TERTILE_BINS,
        title_text="Locomotion and pO$_2$ by locomotion peak height"),
    "sleep_po2_by_transition": lambda: fig_binned_event_responses(
        params=SLEEP_STATE_PARAMS, metric="event_height", bins=SLEEP_TRANSITION_BINS, percentile=False,
        bin_names=SLEEP_TRANSITION_NAMES, title_text="Sleep state and pO$_2$ by transition type"),
}


if __name__ == "__main__":
    for name, make in FIGURES.items():
        fig = make()
        print(name, "->", ", ".join(save(fig, name)))
        plt.close(fig)
