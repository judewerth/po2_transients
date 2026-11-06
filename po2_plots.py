"""
po2_plots.py: FIGURE layer of the event-triggered pO2 transient pipeline.

Every function draws tables built by po2_transients (pt) onto matplotlib axes. The same functions serve the
analysis notebooks (static figures) and po2_gui (redrawn on every widget change):

    po2_transients (tables)  <--  po2_plots (figures)  <--  po2_gui (widgets)

Conventions
-----------
* Functions named `plot_*` that take `ax` (or `ax_*`) draw into existing axes and return nothing (or the twin
  axis they created). The GUI clears and redraws those axes.
* Functions that build a whole figure take `axs=None` / `ax=None`, create the figure when needed, and return it.
* No analysis is done here: every number comes from a pt table or a pt function (event_responses,
  group_average, interpolate_event, ...). To change what is computed, edit po2_transients.py.

Sections
--------
1. Style        colors, labels, color helpers
2. Single run   detection panel, per-signal panels, pO2 point bars / traces / event windows (RunGUI panels)
3. Pooled       QC histograms, pooled / group / class / per-point averages, interpolation check
4. Bins         histogram, bin bars, bin-average traces, correlation scatter (AllRunsGUI panels)
5. Statistics   mixed-model forest plot and variance decomposition
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import po2_transients as pt


# =============================================================================
# 1. Style
# =============================================================================

COLORS = {
    "raw": "slategrey",        # raw event-source signal
    "filt": "dodgerblue",      # filtered signal that is not the detection signal
    "dxdt": "firebrick",       # detection signal
    "event_peak": "black",     # event markers
    "window": "bisque",        # extraction window shading
    "lobe": "salmon",          # positive lobe of an event (event_bounds)
    "threshold": "forestgreen",
    "trace": "gray",           # individual event windows
    "mean": "black",           # averages
}

# axis labels of signals (filtered event-source signals, pO2 variants and the event-window pO2 signals)
SIGNAL_LABELS = {"po2": "pO2", "abs_po2": "|pO2|", "norm_po2": "pO2 / point mean",
                 "po2_norm": "(pO2 - point RMS) / point RMS", "pupil": "Filtered Pupil",
                 "dxdt_pupil": "dPupil/dt", "locom": "Smoothed Locomotion",
                 "sleep_state": "Filtered Sleep State", "dxdt_sleep_state": "dSleepState/dt"}
# short titles of the pO2 signals offered in the GUIs
PO2_TITLES = {"po2": "Raw pO2", "norm_po2": "Normalized pO2", "po2_norm": "Normalized pO2"}

RAW_PREFIX = "raw_"   # "raw_<src>" names the unfiltered event source (dataset table) instead of the filtered one

# table -> what one of its rows is called (titles and dropdown labels)
TABLE_LABELS = {"events": "event", "po2_signals": "po2_signal", "po2_events": "po2_event"}
TABLE_ROW_UNITS = {"events": "events", "po2_signals": "run x point rows", "po2_events": "event x point rows"}

UNBINNED_COLOR = "lightgray"   # values / dots in no bin
MULTI_BIN_COLOR = "black"      # correlation dots that match binned rows of different bins
MAX_INTEGER_HIST_BARS = 50     # integer-valued histograms with at most this range get one bar per integer
BAR_JITTER = 0.25              # half-width of the horizontal dot jitter in bar / category plots
EVENT_LINE_COLORS = {True: "black", False: "darkgray"}   # included / excluded event lines (RunGUI)


def signal_label(signal):
    """Axis / legend label of a filtered signal, a pO2 signal, or a raw_<src> signal."""
    if signal.startswith(RAW_PREFIX):
        return f"Raw {signal[len(RAW_PREFIX):]}"
    return SIGNAL_LABELS.get(signal, signal)


def event_xlabel(params):
    """x-axis label for event-locked time."""
    return f"Time relative to {params.event_method} event (s)"


def default_signals(params, extra=("po2",)):
    """Detection signal first, then the method's other signals, then `extra` (e.g. pO2)."""
    method = pt.EVENT_METHODS[params.event_method]
    return [method.detect_on] + [s for s in method.signals if s != method.detect_on] + list(extra)


def point_colors(n_points):
    """One color per pO2 point, shared by the RunGUI bars, traces and event windows."""
    cmap = plt.get_cmap("tab10" if n_points <= 10 else "tab20")
    return [cmap(i % cmap.N) for i in range(n_points)]


def bin_colors(n_bins):
    """One color per bin (tab10), shared by every AllRunsGUI panel and the bound inputs."""
    cmap = plt.get_cmap("tab10")
    return [cmap(i % cmap.N) for i in range(n_bins)]


def _mean_band(ax, t, mean, std, label, color=None, linewidth=2):
    """Line at `mean` with a ±std band in the same color. Returns the line."""
    line, = ax.plot(t, mean, linewidth=linewidth, label=label, color=color)
    ax.fill_between(t, mean - std, mean + std, color=line.get_color(), alpha=0.2)
    return line


def _empty(ax, text):
    """Gray centered hint for an axis with nothing to show."""
    ax.text(0.5, 0.5, text, ha="center", va="center", transform=ax.transAxes, color="gray")


# =============================================================================
# 2. Single run (inputs: `data.select_run(session, run)`, i.e. a TransientData holding one run)
# =============================================================================

def run_signal(run, signal):
    """
    (t, values) of one signal of a one-run TransientData: raw_<src> comes from `dataset`, anything else from
    `extract_signals` (the filtered signals of the event method).
    """
    if signal.startswith(RAW_PREFIX):
        src = signal[len(RAW_PREFIX):]
        ds = run.dataset.iloc[0]
        return ds["t_" + src], ds[src]
    ext = run.extract_signals.iloc[0]
    return ext["t_" + signal], ext[signal]


def plot_event_windows(ax, event_times, extract_window_s):
    """Shade the extraction window [t_event + w0, t_event + w1] around each event."""
    for t_event in event_times:
        ax.axvspan(t_event + extract_window_s[0], t_event + extract_window_s[-1], color=COLORS["window"], alpha=0.3)


def plot_signal_panel(ax, run, signal, params, show_raw=True):
    """
    One method signal of one run (not the detection signal): the raw source (if `show_raw` and loaded), the
    filtered signal, each event's extraction window and a marker at each event time.
    """
    t, x = run_signal(run, signal)
    ds = run.dataset.iloc[0]
    if show_raw and signal in ds.index:
        ax.plot(ds["t_" + signal], ds[signal], label=f"Raw {signal}", color=COLORS["raw"], alpha=0.6, linewidth=2)
    ax.plot(t, x, label=signal_label(signal), color=COLORS["filt"], alpha=0.6, linewidth=2)
    plot_event_windows(ax, run.events["event_time"], params.extract_window_s)
    ax.scatter(run.events["event_time"], np.interp(run.events["event_time"], t, x),
               color=COLORS["event_peak"], zorder=5)
    ax.set_title(f"Detected events on {signal_label(signal)}", fontweight="bold")
    ax.set_ylabel(signal_label(signal))
    ax.legend(loc="upper right", ncol=3)


def plot_detection_panel(ax, run, params, second=None, show_raw=False, show_lobes=True):
    """
    The detection signal of one run (EVENT_METHODS[method].detect_on) with the threshold, each event's peak,
    its positive lobe (event_bounds) and its extraction window. Used by plot_run_events and RunGUI [0, 0].

    Args:
        second     : optional other signal (a method signal, or raw_<src>) drawn on a twin y-axis
        show_raw   : also draw the raw source when the detection signal is a raw signal filtered in place
                     (e.g. locom)
    Returns the twin axis (None without `second`); the GUI removes it before redrawing.
    """
    detect_on = pt.EVENT_METHODS[params.event_method].detect_on
    t, x = run_signal(run, detect_on)
    ds, events = run.dataset.iloc[0], run.events
    label = signal_label(detect_on)

    if show_raw and detect_on in ds.index:
        ax2 = ax.twinx()
        ax2.plot(ds["t_" + detect_on], ds[detect_on], label=f"Raw {detect_on}", color=COLORS["raw"],
                 alpha=0.6, linewidth=2)
        ax2.set_ylabel(f"Raw {detect_on}", color=COLORS["raw"])
    ax.plot(t, x, label=label, color=COLORS["dxdt"], alpha=0.6, linewidth=2)
    plot_event_windows(ax, events["event_time"], params.extract_window_s)
    if show_lobes:
        for bounds in events["event_bounds"]:
            ax.axvspan(*bounds, color=COLORS["lobe"], alpha=0.3)
    ax.scatter(events["event_time"], events["event_height"], color=COLORS["event_peak"], zorder=5, label="Events")
    ax.axhline(run.extract_signals["event_threshold"].iloc[0], color=COLORS["threshold"], linestyle="--",
               linewidth=1, label="Threshold")
    ax.set_ylabel(label, color=COLORS["dxdt"])
    ax.set_xlabel("Time (s)")
    ax.set_xlim(t[0], t[-1])

    twin = None
    handles, labels = ax.get_legend_handles_labels()
    if second:
        twin = ax.twinx()
        color = COLORS["raw" if second.startswith(RAW_PREFIX) else "filt"]
        t2, values = run_signal(run, second)
        twin.plot(t2, values, label=signal_label(second), color=color, alpha=0.6, linewidth=2)
        twin.set_ylabel(signal_label(second), color=color)
        ax.set_zorder(twin.get_zorder() + 1)   # keep the events and threshold above the second signal
        ax.patch.set_visible(False)
        h2, l2 = twin.get_legend_handles_labels()
        handles, labels = handles + h2, labels + l2
    ax.legend(handles, labels, loc="best", ncol=2, framealpha=0.7, fontsize=8)
    return twin


def plot_run_events(data, session, run, params, axs=None):
    """
    Event-detection QC figure for one run: one panel per method signal (EVENT_METHODS[method].signals).
    The detection-signal panel adds the threshold, peaks and lobes (plot_detection_panel).
    """
    one = data.select_run(session, run)
    method = pt.EVENT_METHODS[params.event_method]
    if axs is None:
        fig, axs = plt.subplots(len(method.signals), 1, figsize=(10, 3.5 * len(method.signals) + 0.5),
                                sharex=True, squeeze=False)
        axs = axs[:, 0]
    fig = axs[0].figure

    for ax, signal in zip(axs, method.signals):
        if signal == method.detect_on:
            plot_detection_panel(ax, one, params, show_raw=True)
            ax.set_xlabel("")
            ax.set_title(f"Events detected from peaks of {signal_label(signal)}", fontweight="bold")
        else:
            plot_signal_panel(ax, one, signal, params)
    t_detect, _ = run_signal(one, method.detect_on)
    axs[-1].set_xlabel("Time (s)")
    axs[-1].set_xlim(t_detect[0], t_detect[-1])

    ds = one.dataset.iloc[0]
    fig.suptitle(f"Filtered {params.event_method} signal for event detection\n({params.describe()})\n"
                 f"Mouse: {ds['mouse_id']}, Date: {ds['date']}, Run: {ds['run']}")
    fig.tight_layout()
    return fig


def plot_point_bars(ax_rms, ax_cv, po2_signals, selected, colors):
    """
    RMS (top) and CV (bottom) of every pO2 point of one run (`po2_signals` rows). Selected points are drawn in
    their point color with a black edge and a bold colored tick label; the rest are faded gray.
    """
    points = po2_signals["point"].to_numpy()
    is_sel = np.isin(points, list(selected))
    face = [colors[i] if s else "lightgray" for i, s in zip(points, is_sel)]
    edge = ["black" if s else "none" for s in is_sel]
    for ax, metric, title in ((ax_rms, "rms", "Point RMS"), (ax_cv, "cv", "Point CV")):
        ax.bar(points, po2_signals[metric], color=face, edgecolor=edge, linewidth=1.5)
        ax.set_ylabel(metric.upper())
        ax.set_title(title, fontsize=10)
    ax_cv.set_xticks(points)
    ax_cv.set_xlim(-0.6, len(points) - 0.4)
    ax_cv.set_xlabel("pO2 point (click to select)")
    for i, label in zip(points, ax_cv.get_xticklabels()):
        if is_sel[i]:
            label.set_color(colors[i])
            label.set_fontweight("bold")
    ax_rms.tick_params(labelbottom=False)


def plot_po2_traces(ax, po2_signals, variant, selected, colors, event_times, included, extract_window_s):
    """
    Full-run pO2 (PO2_VARIANTS `variant` column of one run's `po2_signals`) of each selected point in its point
    color, over thin event lines (black if `included`, gray if excluded) and the included events' windows.
    """
    event_times = np.asarray(event_times)
    included = np.asarray(included, dtype=bool)
    t_po2 = po2_signals["t_po2"].iloc[0]
    values = np.column_stack(po2_signals[variant].to_list())   # (n_time, n_points)

    plot_event_windows(ax, event_times[included], extract_window_s)
    for t_event, inc in zip(event_times, included):
        ax.axvline(t_event, color=EVENT_LINE_COLORS[bool(inc)], linewidth=0.8)
    if not selected:
        _empty(ax, "Click a point in the bar plots to show its pO2")
    for i in sorted(selected):
        ax.plot(t_po2, values[:, i], color=colors[i], linewidth=1)
    if selected:
        # y-limits from the traces alone, so the event shading / lines never stretch them
        sel = values[:, sorted(selected)]
        lo, hi = np.nanmin(sel), np.nanmax(sel)
        pad = 0.05 * (hi - lo) if hi > lo else 1
        ax.set_ylim(lo - pad, hi + pad)
    ax.set_xlim(t_po2[0], t_po2[-1])
    ax.set_xlabel("Time (s)")
    ax.set_ylabel(SIGNAL_LABELS[variant])
    ax.set_title(f"{PO2_TITLES[variant]} ({len(selected)}/{values.shape[1]} points, "
                 f"{included.sum()}/{len(event_times)} events; click an event line to toggle)", fontweight="bold")


def plot_po2_events(ax, traces, variant, colors, n_events, n_points, n_detected, extract_window_s, event_method):
    """
    pO2 event windows of one run (rows of pt.po2_event_traces, already restricted to the included events and
    selected points), each in its point color, plus their mean (thick black, pt.trace_mean_std).

    Args:
        n_events, n_points : how many events / points `traces` covers (title)
        n_detected         : number of detected events, to pick the hint when nothing is shown
    """
    ax.axvline(0, color="black", linewidth=0.8)
    for row in traces.itertuples():
        ax.plot(row.t_event, row.trace, color=colors[row.point], linewidth=1, alpha=0.6)
    if len(traces):
        mean, _, _ = pt.trace_mean_std(traces["trace"])
        ax.plot(traces["t_event"].iloc[0], mean, color="black", linewidth=2.5, label="Mean")
        ax.legend(loc="upper right", fontsize=8, framealpha=0.7)
    elif n_detected == 0:
        _empty(ax, "No events detected")
    elif n_events == 0:
        _empty(ax, "All events excluded (click a gray event line in the pO2 traces to include it)")
    else:
        _empty(ax, "Click a point in the bar plots to show its pO2 events")
    # after clear(), the x-only axvline keeps the previous y data limits; recompute them from scratch
    ax.relim()
    ax.autoscale_view()
    ax.set_xlim(extract_window_s[0], extract_window_s[-1])
    ax.set_xlabel(f"Time relative to {event_method} event (s)")
    ax.set_ylabel(SIGNAL_LABELS[variant])
    ax.set_title(f"{PO2_TITLES[variant]} event responses ({n_events} events x {n_points} points)",
                 fontweight="bold")


# =============================================================================
# 3. Pooled over runs
# =============================================================================

def plot_qc_histograms(points, events, max_cv, n_bins=100, event_label="dPupil/dt"):
    """
    Distributions of pO2 point CV (with the QC cutoff), of pO2 RMS for points passing QC, and of event peak
    heights of the detection signal (`event_label`, e.g. signal_label(EVENT_METHODS[method].detect_on)).
    `points` is pt.point_qc(data).
    """
    fig, axs = plt.subplots(3, 1, figsize=(12, 8))

    axs[0].hist(points["cv"], bins=n_bins, color="skyblue", edgecolor="black")
    axs[0].axvline(max_cv, color="red", linestyle="--", label=f"Max CV: ({max_cv})")
    axs[0].set_title("pO2 Error Across All Runs")
    axs[0].set_xlabel("Error (CV)")
    axs[0].legend()

    axs[1].hist(points.loc[points["cv"] <= max_cv, "rms"], bins=n_bins, color="lightgreen", edgecolor="black")
    axs[1].set_title(f"pO2 Magnitude Across All Runs (points with CV <= {max_cv})")
    axs[1].set_xlabel("Magnitude (RMS)")

    axs[2].hist(events["event_height"], bins=max(n_bins // 5, 1), color="lightcoral", edgecolor="black")
    axs[2].set_title(f"{event_label} Event Peak Across All Runs")
    axs[2].set_xlabel(f"Peak {event_label} (event_height)")

    for ax in axs:
        ax.set_ylabel("Frequency")
    fig.tight_layout()
    return fig


def plot_pooled_events(data, params, signals=None, dt=1.0, points=None):
    """
    Every event window (gray; every point for pO2) of each signal, overlaid with the event average (black,
    ±SD band) after interpolation onto the dt grid (pt.event_responses -> pt.group_average).
    signals defaults to the method's signals + pO2. `points`: optional pt.point_qc subset.
    """
    signals = signals or default_signals(params)
    fig, axs = plt.subplots(len(signals), 1, figsize=(10, 4 * len(signals)), sharex=True, squeeze=False)
    for ax, signal in zip(axs[:, 0], signals):
        for _, t_rel, trace in pt.event_traces(data, signal, points):
            ax.plot(t_rel, trace, color=COLORS["trace"], alpha=0.3)
        avg = pt.group_average(pt.event_responses(data, signal, params.extract_window_s, dt, points=points))
        _mean_band(ax, avg["t"], avg["mean"], avg["std"], color=COLORS["mean"], linewidth=1,
                   label=f"Mean of {avg['n'].max()} events (interpolated to {dt:g} s)")
        ax.axvline(0, color="k", linestyle=":", linewidth=1)
        ax.set_title(f"{signal_label(signal)} Events")
        ax.set_ylabel(signal_label(signal))
        ax.legend(loc="upper right")
    axs[-1, 0].set_xlabel(event_xlabel(params))
    fig.tight_layout()
    return fig


def plot_group_average(avg, by, signal, ax=None, xlabel="Time relative to event (s)"):
    """
    Mean ± SD per group from pt.group_average(..., by=by), e.g. by="mouse_id":
        avg = pt.group_average(pt.event_responses(data, "po2", window), data.events[pt.EVENT_KEYS + ["mouse_id"]],
                               by="mouse_id")
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 4))
    for group, g in avg.groupby(by, observed=True):
        label = group if isinstance(group, str) else ", ".join(map(str, np.atleast_1d(group)))
        _mean_band(ax, g["t"], g["mean"], g["std"], linewidth=1.5, label=f"{label} (n={g['n'].max()} events)")
    ax.axvline(0, color="k", linestyle=":", linewidth=1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(signal_label(signal))
    ax.set_title(f"{signal_label(signal)} event response by {by}")
    ax.legend()
    ax.figure.tight_layout()
    return ax.figure


def plot_class_averages(data, class_col, params, dt=1.0, signals=None, points=None):
    """
    Mean ± SD event response per class of `class_col` (a data.events column, e.g. from pt.classify_events),
    overlaid, one subplot per signal. Each event counts once. signals defaults to the method's signals +
    pO2 + normalized pO2.
    """
    signals = signals or default_signals(params, extra=("po2", "po2_norm"))
    labels = data.events[pt.EVENT_KEYS + [class_col]]
    fig, axs = plt.subplots(len(signals), 1, figsize=(10, 3.5 * len(signals)), sharex=True, squeeze=False)
    for ax, signal in zip(axs[:, 0], signals):
        responses = pt.event_responses(data, signal, params.extract_window_s, dt, points=points)
        avg = pt.group_average(responses, labels, class_col)
        for cls, g in avg.groupby(class_col, observed=True):
            _mean_band(ax, g["t"], g["mean"], g["std"], label=f"{cls} (n={g['n'].max()} events)")
        ax.axvline(0, color="k", linestyle=":", linewidth=1)
        ax.set_title(f"{signal_label(signal)} event response by {class_col}")
        ax.set_ylabel(signal_label(signal))
        ax.legend(loc="upper right", ncol=3)
    axs[-1, 0].set_xlabel(f"{event_xlabel(params)}, interpolated to {dt} s")
    fig.tight_layout()
    return fig


def plot_point_class_averages(data, class_col, params, dt=1.0, normalize=True, points=None):
    """
    pO2 response per class of `class_col` (mean ± SD across events, overlaid), one subplot per pO2 point index,
    pooling every mouse, date and run (point i of every run is averaged together).
    Computed as pt.group_average(pt.event_responses(level="point"), events labels, by=["point", class_col]).
    """
    signal = "po2_norm" if normalize else "po2"
    responses = pt.event_responses(data, signal, params.extract_window_s, dt, level="point", points=points)
    avg = pt.group_average(responses, data.events[pt.EVENT_KEYS + [class_col]], by=["point", class_col])
    point_ids = np.sort(avg["point"].unique())
    ylabel = SIGNAL_LABELS[signal]

    fig, axs = plt.subplots(len(point_ids), 1, figsize=(10, 2.2 * len(point_ids) + 1), sharex=True, sharey=True,
                            squeeze=False, layout="constrained")   # keeps the suptitle clear of the first panel
    for ax, point in zip(axs[:, 0], point_ids):
        for cls, g in avg[avg["point"] == point].groupby(class_col, observed=True):
            _mean_band(ax, g["t"], g["mean"], g["std"], label=f"{cls} (n={g['n'].max()} events)")
        ax.axvline(0, color="k", linestyle=":", linewidth=1)
        ax.set_title(f"Point {point}", loc="left", fontsize=10)
        ax.set_ylabel(ylabel)
    axs[0, 0].legend(loc="upper right", ncol=3, fontsize=8)
    axs[-1, 0].set_xlabel(f"{event_xlabel(params)}, interpolated to {dt} s")
    fig.suptitle(f"{ylabel} event response per point by {class_col}\n"
                 f"Pooled over {data.dataset['mouse_id'].nunique()} mice, {data.dataset['date'].nunique()} dates, "
                 f"{len(data.dataset)} runs")
    return fig


def plot_interpolation_examples(data, class_col, params, dt=1.0, normalize=True, points=None, seed=None):
    """
    4x4 check of pt.interpolate_event: 16 random event x point rows, each showing the observed pO2 window
    (`po2_events`, at the run's own fs) overlaid with its interpolated trace on the shared dt grid.

    Args:
        points : optional pt.point_qc subset; only those points are sampled
        seed   : random seed, so the same 16 rows can be redrawn
    """
    trace_col = "event_trace_norm" if normalize else "event_trace"
    rows = data.po2_events.merge(data.events[pt.EVENT_KEYS + [class_col]], on=pt.EVENT_KEYS, how="left")
    rows = pt._keep_points(rows, points)
    interp = [pt.interpolate_event(r.t_event, getattr(r, trace_col), params.extract_window_s, dt)
              for r in rows.itertuples()]
    rows = rows.assign(t_grid=[t for t, _ in interp], interpolated=[v for _, v in interp])
    rows = rows[[np.isfinite(v).any() for v in rows["interpolated"]]]   # skip all-NaN windows
    rows = rows.iloc[np.sort(np.random.default_rng(seed).choice(len(rows), min(16, len(rows)), replace=False))]
    ylabel = SIGNAL_LABELS["po2_norm" if normalize else "po2"]

    fig, axs = plt.subplots(4, 4, figsize=(16, 12), sharex=True, layout="constrained")
    for ax, row in zip(axs.flat, rows.itertuples()):
        ax.plot(row.t_event, getattr(row, trace_col), "o-", color=COLORS["raw"], markersize=4, linewidth=1,
                alpha=0.7, label="Observed (run fs)")
        ax.plot(row.t_grid, row.interpolated, ".-", color=COLORS["filt"], markersize=4, linewidth=1.5,
                label=f"Interpolated ({dt} s)")
        ax.axvline(0, color="k", linestyle=":", linewidth=1)
        ax.set_title(f"{row.mouse_id} {row.date} run {row.run}\nevent {row.event}, point {row.point}"
                     f" ({getattr(row, class_col)})", fontsize=9)
    for ax in axs.flat[len(rows):]:
        ax.set_visible(False)
    for ax in axs[:, 0]:
        ax.set_ylabel(ylabel)
    for ax in axs[-1, :]:
        ax.set_xlabel(event_xlabel(params))
    fig.legend(*axs.flat[0].get_legend_handles_labels(), loc="outside upper right", ncol=2)
    fig.suptitle(f"Observed vs. interpolated pO2 for {len(rows)} random event x point windows")
    return fig


# =============================================================================
# 4. Bins (the AllRunsGUI panels; inputs from pt.bin_edges / bin_labels / bar_values / correlation_values)
# =============================================================================

def sample_title(data, params):
    """
    Two-line figure title: the pooled sample (runs, events, pO2 points; each mouse with its dates and run count)
    and the event detection method with its parameters.
    """
    runs = data.dataset[pt.RUN_KEYS]
    mice = [f"{mouse} ({len(g)} runs: {', '.join(sorted(g['date'].astype(str).unique()))})"
            for mouse, g in runs.groupby("mouse_id", sort=True)]
    sample = (f"{len(runs)} runs, {runs['date'].nunique()} dates, {runs['mouse_id'].nunique()} mice | "
              f"{len(data.events)} events, {len(data.po2_signals)} pO2 points | " + "; ".join(mice))
    detect_on = pt.EVENT_METHODS[params.event_method].detect_on
    method = f"Events: {params.event_method} method, peaks of {signal_label(detect_on)} | {params.describe()}"
    return f"{sample}\n{method}"


def plot_bin_histogram(ax, values, edges, colors, metric, n_runs, n_hist=20):
    """
    Stacked histogram of `values` (one BIN_METRICS metric pooled over runs), each value colored by its bin
    (edges from pt.bin_edges; unbinned values gray), with dashed lines at every bin's [lo, hi].
    """
    table = pt.BIN_METRICS[metric]
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    ax.set_title(f"{metric} ({TABLE_LABELS[table]}): {finite.sum()} {TABLE_ROW_UNITS[table]} across {n_runs} runs",
                 fontweight="bold")
    ax.set_xlabel(metric)
    ax.set_ylabel("Count")
    if not finite.any():
        _empty(ax, f"No finite {metric} values")
        return

    v = values[finite]
    if np.all(v == np.round(v)) and np.ptp(v) <= MAX_INTEGER_HIST_BARS:   # integers (e.g. point): one bar each
        hist_edges = np.arange(v.min() - 0.5, v.max() + 1.5)
    else:
        hist_edges = np.histogram_bin_edges(v, bins=n_hist)
    idx = pt.bin_index(values, edges)
    stacks, cols, labels = [], [], []
    for k, ((lo, hi), color) in enumerate(zip(edges, colors)):
        stacks.append(values[idx == k])
        cols.append(color)
        labels.append(f"Bin {k + 1}: {lo:.3g} to {hi:.3g} (n={np.sum(idx == k)})")
    unbinned = values[finite & (idx == -1)]
    if len(unbinned) or not stacks:
        stacks.append(unbinned)
        cols.append(UNBINNED_COLOR)
        labels.append(f"No bin (n={len(unbinned)})")
    ax.hist(stacks, bins=hist_edges, stacked=True, color=cols, label=labels, edgecolor="white", linewidth=0.3)

    for (lo, hi), color in zip(edges, colors):
        for edge in (lo, hi):
            if np.isfinite(edge):
                ax.axvline(edge, color=color, linestyle="--", linewidth=1.2)
    ax.set_xlim(hist_edges[0], hist_edges[-1])
    ax.legend(loc="upper right", fontsize=8, framealpha=0.7)


def plot_bin_traces(ax, avg, colors, ylabel, unit="events"):
    """
    Mean ± std event-locked trace per bin (pt.group_average with by="bin", an ordered categorical), bin k in
    colors[k] so the bands match the histogram. `unit` names what n counts ("events" or "rows").
    """
    ax.axvline(0, color="black", linestyle=":", linewidth=1)
    bounds = []   # mean ± std of every drawn bin, for the y-limits
    for k, (label, color) in enumerate(zip(avg["bin"].cat.categories, colors)):
        g = avg[avg["bin"] == label]
        if g.empty or not np.isfinite(g["mean"]).any():
            continue
        lo, hi = g["mean"] - g["std"].fillna(0), g["mean"] + g["std"].fillna(0)
        ax.plot(g["t"], g["mean"], color=color, linewidth=2, label=f"Bin {k + 1}: {label} (n={int(g['n'].max())} {unit})")
        ax.fill_between(g["t"], lo, hi, color=color, alpha=0.2, linewidth=0)
        bounds.extend([np.nanmin(lo), np.nanmax(hi)])
    if bounds:
        # explicit limits: after clear() the axvline keeps stale y data limits, and relim ignores the bands
        lo, hi = min(bounds), max(bounds)
        pad = 0.05 * (hi - lo) if hi > lo else 1
        ax.set_ylim(lo - pad, hi + pad)
    else:
        _empty(ax, "No binned traces")
    ax.set_ylabel(ylabel, fontsize=9)


def plot_bin_bars(ax, values, colors, metric, bin_metric, unit="events"):
    """
    One bar per bin (`values["bin"]`, an ordered categorical, from pt.bar_values) of `metric`: height = mean,
    error bar = std, dots = the individual values (fixed jitter, so they do not move between redraws).
    """
    ax.set_title(f"{metric} by {bin_metric} bin", fontweight="bold")
    ax.set_ylabel(metric)
    categories = values["bin"].cat.categories
    if len(categories) == 0:
        _empty(ax, "No bins")
        return
    rng = np.random.default_rng(0)
    ticks, finite_all = [], []
    for k, (label, color) in enumerate(zip(categories, colors)):
        v = values.loc[values["bin"] == label, metric].to_numpy(dtype=float)
        v = v[np.isfinite(v)]
        ticks.append(f"Bin {k + 1}\n(n={len(v)})")
        if not len(v):
            continue
        std = v.std(ddof=1) if len(v) > 1 else 0.0
        ax.bar(k, v.mean(), yerr=std, capsize=4, width=0.8, color=color, alpha=0.6, edgecolor=color,
               error_kw={"elinewidth": 1.5, "ecolor": "black"})
        ax.scatter(k + rng.uniform(-BAR_JITTER, BAR_JITTER, len(v)), v, s=8, color="black", alpha=0.5,
                   linewidths=0, zorder=3)
        finite_all.append(v)
    ax.set_xticks(range(len(categories)), ticks, fontsize=8)
    ax.set_xlim(-0.6, len(categories) - 0.4)
    if not finite_all:
        _empty(ax, f"No finite {metric} values in any bin")
        return
    finite_all = np.concatenate(finite_all)
    if finite_all.min() < 0 < finite_all.max():
        ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xlabel(f"{bin_metric} bin (dots: {unit})")


def plot_correlation(ax, df, colors, x, y, bin_metric):
    """
    Scatter of df["y"] vs df["x"] (pt.correlation_values), each dot colored by df["color_idx"] (bin k ->
    colors[k], pt.NO_BIN gray, pt.MULTI_BIN black), with a least-squares line and r² over all dots when x is
    numeric. A string x is drawn as jittered categories (no fit).
    """
    ax.set_title(f"{y} vs {x} (colored by {bin_metric} bin)", fontweight="bold")
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    if df.empty:
        _empty(ax, f"No rows with both {x} and {y}")
        return

    numeric = pd.api.types.is_numeric_dtype(df["x"])
    if numeric:
        xs = df["x"].to_numpy(dtype=float)
    else:
        categories = sorted(df["x"].astype(str).unique())
        codes = pd.Categorical(df["x"].astype(str), categories=categories).codes
        xs = codes + np.random.default_rng(0).uniform(-BAR_JITTER, BAR_JITTER, len(df))
        many = len(categories) > 3
        ax.set_xticks(range(len(categories)), categories, fontsize=8, rotation=30 if many else 0,
                      ha="right" if many else "center")
    ys = df["y"].to_numpy(dtype=float)
    idx = df["color_idx"].to_numpy()

    groups = [(k, colors[k], f"Bin {k + 1}") for k in range(len(colors))]
    groups += [(pt.MULTI_BIN, MULTI_BIN_COLOR, "Multiple bins"), (pt.NO_BIN, UNBINNED_COLOR, "No bin")]
    for k, color, label in groups:
        sel = idx == k
        if sel.any():
            ax.scatter(xs[sel], ys[sel], s=10, color=color, alpha=0.6, linewidths=0, label=f"{label} (n={sel.sum()})")

    if numeric and len(df) >= 3 and np.ptp(xs) > 0:
        slope, intercept = np.polyfit(xs, ys, 1)
        r = np.corrcoef(xs, ys)[0, 1]
        x_line = np.array([xs.min(), xs.max()])
        ax.plot(x_line, slope * x_line + intercept, color="dimgray", linestyle="--", linewidth=1.5,
                label=f"Fit: r² = {r ** 2:.3f}, n = {len(xs)}")
    ax.legend(loc="best", fontsize=7, framealpha=0.7)


# =============================================================================
# 5. Statistics (inputs from pt.extract_fixed_effects / pt.extract_variance_decomposition)
# =============================================================================

def plot_fixed_effects(effects, title="Fixed-effect estimates", xlabel="Estimated effect", alpha=0.05):
    """Forest plot: estimate ± CI per term, red if p < alpha, grey otherwise."""
    effects = effects.sort_values("estimate").reset_index(drop=True)
    colors = np.where(effects["pvalue"] < alpha, "#C62828", "#666666")
    y = np.arange(len(effects))

    fig, ax = plt.subplots(figsize=(10, max(3.0, 0.55 * len(effects) + 1.5)))
    ax.hlines(y, effects["ci_lower"], effects["ci_upper"], color=colors, linewidth=2.2, zorder=2)
    ax.scatter(effects["estimate"], y, color=colors, s=65, zorder=3)
    for i, row in effects.iterrows():
        if np.isnan(row["pvalue"]):   # singular fit: no standard error
            p_text = "p = n/a"
        else:
            p_text = "p < 0.001" if row["pvalue"] < 0.001 else f"p = {row['pvalue']:.3f}"
        ax.annotate(p_text, xy=(np.nan_to_num(row["ci_upper"], nan=row["estimate"]), i), xytext=(6, 0),
                    textcoords="offset points", va="center", fontsize=9, color=colors[i])

    ax.axvline(0, color="k", linestyle="--", linewidth=1, alpha=0.65)
    ax.set_yticks(y)
    ax.set_yticklabels(effects["term"])
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontweight="bold")
    ax.grid(axis="x", linestyle=":", alpha=0.35)
    ax.spines[["top", "right", "left"]].set_visible(False)
    fig.tight_layout()
    return fig


def plot_variance_decomposition(variance, title="Variance decomposition"):
    """Percent of the total model variance per component."""
    variance = variance.sort_values("percent_total").reset_index(drop=True)
    palette = {"fixed": "#377EB8", "random": "#984EA3", "residual": "#9E9E9E"}

    fig, ax = plt.subplots(figsize=(9, max(3.0, 0.55 * len(variance) + 1.25)))
    ax.barh(variance["component"], variance["percent_total"], color=variance["component_type"].map(palette),
            edgecolor="black", linewidth=0.5)
    for i, value in enumerate(variance["percent_total"]):
        ax.text(value + 0.7, i, f"{value:.1f}%", va="center", fontsize=9)

    ax.set_xlim(0, max(100, variance["percent_total"].max() * 1.22))
    ax.set_xlabel("Percent of model variance")
    ax.set_title(title, fontweight="bold")
    ax.grid(axis="x", linestyle=":", alpha=0.35)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=color, label=kind.title())
                       for kind, color in palette.items() if kind in variance["component_type"].values],
              frameon=False, loc="lower right")
    fig.tight_layout()
    return fig
