"""
po2_gui.py: INTERACTIVE layer of the pO2 transient pipeline (ipywidgets + ipympl).

This module only holds widgets, GUI state (selections, bins) and callbacks. Every number comes from
po2_transients (pt), and every panel is drawn by po2_plots (pp):

    po2_transients (tables)  <--  po2_plots (figures)  <--  po2_gui (widgets)

Use in a notebook with `%matplotlib widget`:

    gui = RunGUI(runs, params)                  # one run at a time, from the .mat files
    gui.display()

    data = pt.process_all_runs(runs, params)
    all_gui = AllRunsGUI(data, params)          # every processed run pooled
    all_gui.display()

RunGUI figure layout (2x2)                     data (pt)                         drawing (pp)
-------------------------------------------------------------------------------------------------------------
[0, 0] detection panel: detection signal,       pt.load_run (cached per run and   pp.plot_detection_panel
       threshold, event peaks, positive lobes,  method) -> pt.analyze_run
       extraction windows; optionally a second  (re-run on every threshold
       signal on a twin y-axis                  change)
[1, 0] pO2 traces of the selected points, with  run.po2_signals[variant]          pp.plot_po2_traces
       event lines; click a line to exclude
       (gray) / include (black) that event
[1, 1] RMS (top) and CV (bottom) per point;     run.po2_signals[rms, cv]          pp.plot_point_bars
       click a column to (de)select the point
[0, 1] pO2 windows of the included events x     pt.po2_event_traces(run, variant) pp.plot_po2_events
       selected points, and their mean

Controls: run, event method, second signal, pO2 signal, threshold, rel_to_rms; select / deselect all
points buttons sit in the figure above [1, 1]. Changing the threshold or rel_to_rms only re-runs
pt.analyze_run on the cached pt.LoadedRun, so only a new run / method reloads the .mat file. Selecting points
or events only redraws [0, 1], [1, 0] and [1, 1]. Excluded events are kept (by event time) across
re-detections and reset when the run changes.

AllRunsGUI figure layout (2x2)                 data (pt)                         drawing (pp)
-------------------------------------------------------------------------------------------------------------
[0, 0] histogram of the "Bin by" metric, each   pt.metric_table, pt.bin_edges     pp.plot_bin_histogram
       value colored by its bin
[0, 1] mean ± std (+ dots) of the bar metric    pt.bin_labels -> pt.bar_values    pp.plot_bin_bars
       per bin
[1, 0] mean ± std per bin of a behavior signal  pt.event_responses (cached) ->    pp.plot_bin_traces
       (top) and of pO2 (bottom), interpolated  pt.group_average(labels, "bin")
       onto a dt grid. Event metric: each event
       counts once; point metric: each event x
       point row counts once
[1, 1] y vs x scatter colored by bin (black:    pt.correlation_values             pp.plot_correlation
       rows of several bins), linear fit + r²

Controls (above the figure): "Bins:" + adds a bin row, - removes the bottom one; each row's lower/upper bound
is labeled with its bin color. "Bin by" picks the metric (or "point": bins of pO2 point numbers). "Percentiles" switches the bounds between
percentiles of the metric and metric values, and converts the typed bounds so the binning is unchanged.
`AllRunsGUI.classify()` returns pt.classify_events with the GUI's bins, for use downstream.
"""

from dataclasses import replace
from functools import wraps

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from matplotlib.widgets import Button
import ipywidgets as widgets
from IPython.display import display

import po2_transients as pt
import po2_plots as pp


NO_SIGNAL = "None"   # "Second signal" radio value for no twin-axis signal
EVENT_PICK_PX = 5    # click distance (pixels) to an event line in [1, 0] that toggles it
RUN_PO2_VARIANTS = ["po2", "norm_po2"]            # RunGUI "pO2 signal" choices (pt.PO2_VARIANTS names)
EVENT_PO2_SIGNALS = ["po2", "po2_norm"]           # AllRunsGUI "pO2" choices (pt.signal_columns names)
INTERP_DT_S = 1.0    # grid step (s) the AllRunsGUI [1, 0] averages are interpolated onto
# metrics the [0, 1] bars can show: one value per event or per event x point (not per run x point)
BAR_METRICS = [m for m, table in pt.BIN_METRICS.items() if table in ("events", "po2_events")]


def logged(fn):
    """
    Run a GUI method inside its `self.log_w` Output widget, so an exception in a widget / mpl callback shows its
    traceback under the controls instead of silently leaving the figure blank or frozen.
    """
    @wraps(fn)
    def wrapper(self, *args, **kwargs):
        with self.log_w:
            return fn(self, *args, **kwargs)
    return wrapper


def second_signal_options(event_method):
    """Signals that can share the detection panel: the method's other filtered signals, then the raw sources."""
    method = pt.EVENT_METHODS[event_method]
    others = [s for s in method.signals if s != method.detect_on]
    raw = [pp.RAW_PREFIX + s for s in method.raw]
    return [NO_SIGNAL] + others + raw


def metric_option_label(metric):
    """Dropdown label of a pt.BIN_METRICS metric with the table it comes from, e.g. "event_auc  (event)"."""
    return f"{metric}  ({pp.TABLE_LABELS[pt.BIN_METRICS[metric]]})"


# =============================================================================
# Single-run GUI
# =============================================================================

class RunGUI:
    """
    2x2 figure + controls for one run at a time.

    Args:
        runs   : {session folder: [run numbers]}
        params : pt.EventParams; supplies everything not set in the GUI (filters, min distance, window)

    State:
        self.run_data : pt.TransientData of the current run (pt.analyze_run output), redrawn from
        self.selected : set of selected pO2 point indices
        self.excluded : set of excluded event times
    """

    def __init__(self, runs, params, figsize=(11.9, 6.8)):
        self.params = params
        self._cache = {}       # (folder, run, event_method) -> pt.LoadedRun
        self._twin = None      # twin axis of [0, 0] (removed before each redraw)
        self._busy = False     # True while the GUI itself changes widgets, to suppress redraws
        self.selected = set()
        self.excluded = set()

        run_options = [(f"{folder}  run {run:02d}", (folder, run)) for folder, run_list in runs.items() for run in run_list]
        style = {"description_width": "initial"}
        self.run_w = widgets.Dropdown(options=run_options, description="Run", style=style)
        self.method_w = widgets.Dropdown(options=list(pt.EVENT_METHODS), value=params.event_method,
                                         description="Event method", style=style)
        self.second_w = widgets.RadioButtons(options=self._second_options(params.event_method),
                                             value=NO_SIGNAL, description="Second signal", style=style)
        self.thresh_w = widgets.FloatText(value=params.height_thresh, step=0.1, style=style,
                                          layout=widgets.Layout(width="220px"))
        self.rel_w = widgets.Checkbox(value=params.rel_to_rms, description="Threshold relative to RMS", style=style)
        self.po2_w = widgets.RadioButtons(options=[(pp.PO2_TITLES[v], v) for v in RUN_PO2_VARIANTS],
                                          value="po2", description="pO2 signal", style=style)
        self.status_w = widgets.HTML()
        self.log_w = widgets.Output()   # tracebacks of failed callbacks (see `logged`)

        with plt.ioff():
            self.fig, self.axs = plt.subplots(2, 2, figsize=figsize, layout="constrained")
        self.fig.canvas.header_visible = False  # hide ipympl's "Figure 1" title
        # [1, 1] -> a row of select / deselect buttons above two stacked bar plots sharing the point axis
        grid = self.axs[1, 1].get_subplotspec().subgridspec(3, 2, height_ratios=[0.25, 1, 1])
        self.axs[1, 1].remove()
        self.select_all_btn = Button(self.fig.add_subplot(grid[0, 0]), "Select all points",
                                     color="whitesmoke", hovercolor="lightgray")
        self.deselect_all_btn = Button(self.fig.add_subplot(grid[0, 1]), "Deselect all points",
                                       color="whitesmoke", hovercolor="lightgray")
        self.ax_rms = self.fig.add_subplot(grid[1, :])
        self.ax_cv = self.fig.add_subplot(grid[2, :], sharex=self.ax_rms)

        self.method_w.observe(self._on_method, names="value")
        self.run_w.observe(self._on_run, names="value")
        for w in (self.second_w, self.thresh_w, self.rel_w):
            w.observe(self._on_change, names="value")
        self.po2_w.observe(lambda _: self.update_points(), names="value")
        self.select_all_btn.on_clicked(lambda _: self.set_selected(range(self.n_points())))
        self.deselect_all_btn.on_clicked(lambda _: self.set_selected([]))
        self.fig.canvas.mpl_connect("button_press_event", self._on_click)
        self.update()

    # ---- data -------------------------------------------------------------------------------------

    def run_params(self):
        """`params` with the GUI's event method, threshold and rel_to_rms."""
        return replace(self.params, event_method=self.method_w.value,
                       height_thresh=self.thresh_w.value, rel_to_rms=self.rel_w.value)

    def loaded_run(self):
        """pt.load_run (stages 1-2) of the chosen run and method, cached."""
        folder, run = self.run_w.value
        key = (folder, run, self.method_w.value)
        if key not in self._cache:
            self._cache[key] = pt.load_run(folder, run, self.run_params())
        return self._cache[key]

    def n_points(self):
        return self.loaded_run().po2.shape[1]

    def included_events(self):
        """Boolean mask over self.run_data.events: True for events not excluded in [1, 0]."""
        return ~np.isin(self.run_data.events["event_time"], list(self.excluded))

    # ---- callbacks ----------------------------------------------------------------------------------

    @staticmethod
    def _second_options(event_method):
        """(label, value) radio options of second_signal_options."""
        return [(s if s == NO_SIGNAL else pp.signal_label(s), s) for s in second_signal_options(event_method)]

    def _on_method(self, change):
        self._busy = True   # resetting the radio options would otherwise redraw with the old method
        self.second_w.options = self._second_options(change["new"])
        self.second_w.value = NO_SIGNAL
        self._busy = False
        self.update()

    def _on_change(self, change):
        if not self._busy:
            self.update()   # new events -> new windows in [1, 0] too

    def _on_run(self, change):
        self.selected = set()   # point indices belong to the previous run
        self.excluded = set()
        self.update()

    @logged
    def _on_click(self, event):
        """Toggle the clicked point's bar column, or the event line in [1, 0] nearest the click (not while panning)."""
        toolbar = getattr(self.fig.canvas, "toolbar", None)
        if event.xdata is None or (toolbar and toolbar.mode):
            return
        if event.inaxes in (self.ax_rms, self.ax_cv):
            i = int(round(event.xdata))
            if 0 <= i < self.n_points():
                self.set_selected(self.selected ^ {i})
        elif event.inaxes is self.axs[1, 0]:
            self._toggle_event_near(event)

    def _toggle_event_near(self, event):
        """Include / exclude the event line in [1, 0] closest to the click, if within EVENT_PICK_PX."""
        times = self.run_data.events["event_time"].to_numpy()
        if len(times) == 0:
            return
        ax = self.axs[1, 0]
        x_px = ax.transData.transform(np.column_stack([times, np.full(len(times), event.ydata)]))[:, 0]
        i = int(np.argmin(np.abs(x_px - event.x)))
        if abs(x_px[i] - event.x) <= EVENT_PICK_PX:
            self.excluded ^= {times[i]}
            self.update_points()

    def set_selected(self, points):
        """Select exactly `points` and redraw the point panels."""
        self.selected = set(points)
        self.update_points()

    # ---- drawing ------------------------------------------------------------------------------------

    @logged
    def update(self):
        """Re-detect and redraw every panel."""
        self.log_w.clear_output()
        self.update_detection()
        self.update_points()

    def update_detection(self):
        """Re-run pt.analyze_run with the current controls and redraw [0, 0]."""
        params = self.run_params()
        self.run_data = pt.analyze_run(self.loaded_run(), params)
        events = self.run_data.events
        self.excluded &= set(events["event_time"].tolist())   # drop events that are no longer detected
        self.thresh_w.description = "Threshold (x RMS)" if params.rel_to_rms else "Threshold (abs)"

        ax = self.axs[0, 0]
        if self._twin is not None:
            self._twin.remove()
        ax.clear()
        second = None if self.second_w.value == NO_SIGNAL else self.second_w.value
        self._twin = pp.plot_detection_panel(ax, self.run_data, params, second=second)

        keys = self.run_data.dataset.iloc[0]
        detect_on = pt.EVENT_METHODS[params.event_method].detect_on
        threshold = self.run_data.extract_signals["event_threshold"].iloc[0]
        ax.set_title(f"Events detected from peaks of {pp.signal_label(detect_on)}\n"
                     f"Mouse: {keys['mouse_id']}, Date: {keys['date']}, Run: {keys['run']} | "
                     f"{len(events)} events, threshold {threshold:.3g}", fontweight="bold")
        self.status_w.value = f"<small>{params.describe()}</small>"
        self.fig.canvas.draw_idle()

    @logged
    def update_points(self):
        """Redraw the point bars [1, 1], the selected points' traces [1, 0] and their included event windows [0, 1]."""
        run, variant = self.run_data, self.po2_w.value
        colors = pp.point_colors(self.n_points())
        for ax in (self.axs[0, 1], self.axs[1, 0], self.ax_rms, self.ax_cv):
            ax.clear()

        event_times = run.events["event_time"].to_numpy()
        included = self.included_events()
        pp.plot_point_bars(self.ax_rms, self.ax_cv, run.po2_signals, self.selected, colors)
        pp.plot_po2_traces(self.axs[1, 0], run.po2_signals, variant, self.selected, colors,
                           event_times, included, self.params.extract_window_s)

        traces = pt.po2_event_traces(run, variant)
        keep = traces["event"].isin(run.events["event"][included]) & traces["point"].isin(self.selected)
        pp.plot_po2_events(self.axs[0, 1], traces[keep].sort_values("point", kind="stable"), variant, colors,
                           n_events=int(included.sum()), n_points=len(self.selected), n_detected=len(event_times),
                           extract_window_s=self.params.extract_window_s, event_method=self.method_w.value)
        self.fig.canvas.draw_idle()

    def display(self):
        """Show the controls above the figure (needs `%matplotlib widget`)."""
        controls = widgets.HBox([
            widgets.VBox([self.run_w, self.method_w]),
            self.second_w,
            self.po2_w,
            widgets.VBox([self.thresh_w, self.rel_w]),
        ])
        display(widgets.VBox([controls, self.status_w, self.log_w, self.fig.canvas]))


# =============================================================================
# All-runs GUI
# =============================================================================

class AllRunsGUI:
    """
    2x2 figure + controls for every processed run at once (layout in the module docstring).

    Args:
        data   : pt.TransientData from pt.process_all_runs
        params : the pt.EventParams it was processed with
        metric, bins, percentile : initial "Bin by" metric and [lo, hi] bounds (percentiles or values)
        dt     : grid step (s) for the [1, 0] interpolated averages
        bar_metric, corr_x, corr_y : initial [0, 1] / [1, 1] choices

    Downstream: `bins()`, `edges()`, `metric_w.value`, `percentile_w.value`, or `classify()`.
    """

    def __init__(self, data, params, metric="event_auc", bins=pt.DEFAULT_BINS, percentile=True,
                 dt=INTERP_DT_S, bar_metric="event_peak", corr_x="rms", corr_y="event_peak",
                 figsize=(11.9, 6.8)):
        self.data = data
        self.params = params
        self.dt = dt
        self._busy = False
        self.log_w = widgets.Output()   # tracebacks of failed callbacks (see `logged`)
        self._rows = []        # one dict per bin row: box, swatch, lo, hi
        self._responses = {}   # (level, signal) -> pt.event_responses output (slow, so cached)

        style = {"description_width": "initial"}
        small = widgets.Layout(width="34px")
        metric_options = [(metric_option_label(m), m) for m in pt.BIN_METRICS]
        self.add_w = widgets.Button(description="+", tooltip="Add a bin", layout=small)
        self.remove_w = widgets.Button(description="−", tooltip="Remove the bottom bin", layout=small)
        self.metric_w = widgets.Dropdown(options=metric_options, value=metric, description="Bin by", style=style)
        self.percentile_w = widgets.Checkbox(value=percentile, description="Percentiles", style=style,
                                             indent=False, layout=widgets.Layout(width="110px"))
        self.rows_w = widgets.VBox()
        self.behavior_w = widgets.Dropdown(
            options=[(pp.signal_label(s), s) for s in pt.EVENT_METHODS[params.event_method].signals],
            description="Behavior", style=style)
        self.po2_w = widgets.Dropdown(options=[(pp.PO2_TITLES[s], s) for s in EVENT_PO2_SIGNALS],
                                      value="po2", description="pO2", style=style)
        self.bar_metric_w = widgets.Dropdown(options=[(metric_option_label(m), m) for m in BAR_METRICS],
                                             value=bar_metric, description="Metric", style=style)
        key_options = [(k, k) for k in pt.CORR_X_KEYS if k not in pt.BIN_METRICS]   # "point" is already a metric
        self.corr_x_w = widgets.Dropdown(options=key_options + metric_options,
                                         value=corr_x, description="x", style=style)
        self.corr_y_w = widgets.Dropdown(options=metric_options, value=corr_y, description="y", style=style)

        with plt.ioff():
            self.fig, self.axs = plt.subplots(2, 2, figsize=figsize, layout="constrained")
        self.fig.canvas.header_visible = False
        self.fig.suptitle(pp.sample_title(data, params), fontsize=10)
        # [1, 0] -> behavior signal (top) and pO2 (bottom) bin averages sharing the event-time axis
        grid = self.axs[1, 0].get_subplotspec().subgridspec(2, 1)
        self.axs[1, 0].remove()
        self.ax_behavior = self.fig.add_subplot(grid[0])
        self.ax_po2 = self.fig.add_subplot(grid[1], sharex=self.ax_behavior)

        for lo, hi in np.asarray(bins, dtype=float).reshape(-1, 2):
            self._add_row(lo, hi)
        self.add_w.on_click(self._on_add)
        self.remove_w.on_click(self._on_remove)
        self.metric_w.observe(self._on_metric, names="value")
        self.percentile_w.observe(self._on_percentile, names="value")
        for w in (self.behavior_w, self.po2_w):
            w.observe(lambda _: self.update_bin_traces(), names="value")
        self.bar_metric_w.observe(lambda _: self.update_bin_bars(), names="value")
        for w in (self.corr_x_w, self.corr_y_w):
            w.observe(lambda _: self.update_correlation(), names="value")
        self.update_bins()

    # ---- data (thin wrappers around pt) -------------------------------------------------------------

    def metric_table(self, metric=None):
        """The TransientData table holding `metric` (default: the "Bin by" metric)."""
        return pt.metric_table(self.data, metric or self.metric_w.value)

    def metric_values(self, metric=None):
        metric = metric or self.metric_w.value
        return self.metric_table(metric)[metric].to_numpy(dtype=float)

    def bins(self):
        """(n_bins, 2) [lo, hi] bounds as typed (percentiles if `percentile_w`, else metric values)."""
        return np.array([[r["lo"].value, r["hi"].value] for r in self._rows], dtype=float).reshape(-1, 2)

    def edges(self):
        """(n_bins, 2) bin bounds in metric units."""
        return pt.bin_edges(self.metric_values(), self.bins(), self.percentile_w.value)

    def classify(self):
        """pt.classify_events of the "Bin by" metric's table with the GUI's bins (Series on that table's index)."""
        return pt.classify_events(self.metric_table(), self.metric_w.value, self.bins(), self.percentile_w.value)

    def bin_labels(self):
        """pt.bin_labels with the GUI's bins: (labels, per_point)."""
        return pt.bin_labels(self.data, self.metric_w.value, self.bins(), self.percentile_w.value)

    def responses(self, level, signal):
        """pt.event_responses at `level` ("event" / "point"), cached per (level, signal)."""
        key = (level, signal)
        if key not in self._responses:
            self._responses[key] = pt.event_responses(self.data, signal, self.params.extract_window_s, self.dt,
                                                      level=level)
        return self._responses[key]

    # ---- bin rows -----------------------------------------------------------------------------------

    def _add_row(self, lo, hi):
        text = widgets.Layout(width="90px")
        row = {"swatch": widgets.HTML(), "lo": widgets.FloatText(value=lo, layout=text),
               "hi": widgets.FloatText(value=hi, layout=text)}
        row["box"] = widgets.HBox([row["swatch"], row["lo"], widgets.Label("to"), row["hi"]])
        for key in ("lo", "hi"):
            row[key].observe(self._on_bound, names="value")
        self._rows.append(row)
        self._refresh_rows()

    def _refresh_rows(self):
        """Re-show the bin rows with their colors (the colors follow the row order)."""
        for k, (row, color) in enumerate(zip(self._rows, pp.bin_colors(len(self._rows)))):
            row["swatch"].value = (f"<span style='display:inline-block;width:14px;height:14px;"
                                   f"background:{to_hex(color)};border:1px solid #333;vertical-align:middle'>"
                                   f"</span>&nbsp;<b>Bin {k + 1}</b>")
        self.rows_w.children = [r["box"] for r in self._rows]

    def _set_bounds(self, bounds):
        """Write `bounds` into the rows without a redraw per value."""
        self._busy = True
        for row, (lo, hi) in zip(self._rows, bounds):
            row["lo"].value, row["hi"].value = self._round(lo, np.floor), self._round(hi, np.ceil)
        self._busy = False

    def _round(self, v, how=np.round):
        """
        Bound shown in a row: percentiles to 0.1, metric values to 4 significant digits. `how` (np.floor for lower
        bounds, np.ceil for upper ones) widens converted bins so they keep every value they held.
        """
        v = float(v)
        if not np.isfinite(v) or v == 0:
            return v
        decimals = 1 if self.percentile_w.value else 3 - int(np.floor(np.log10(abs(v))))
        return float(round(how(v * 10.0 ** decimals) / 10.0 ** decimals, max(decimals, 0)))

    # ---- callbacks ----------------------------------------------------------------------------------

    def _on_add(self, _):
        """New bin from the previous bin's upper bound to the top (the full range if that would be empty)."""
        if self.percentile_w.value:
            vmin, vmax = 0.0, 100.0
        else:
            values = self.metric_values()
            vmin, vmax = (np.nanmin(values), np.nanmax(values)) if np.isfinite(values).any() else (0.0, 0.0)
        prev_hi = self._rows[-1]["hi"].value if self._rows else None
        lo = prev_hi if prev_hi is not None and vmin <= prev_hi and self._round(prev_hi) < self._round(vmax) else vmin
        self._add_row(self._round(lo, np.floor), self._round(vmax, np.ceil))
        self.update_bins()

    def _on_remove(self, _):
        if self._rows:
            self._rows.pop()
            self._refresh_rows()
            self.update_bins()

    def _on_bound(self, _):
        if not self._busy:
            self.update_bins()

    def _on_percentile(self, change):
        """Convert the typed bounds to the new units so the binning stays the same."""
        values = self.metric_values()
        if change["new"]:
            self._set_bounds(pt.values_to_percentiles(self.bins(), values))
        else:
            self._set_bounds(pt.bin_edges(values, self.bins(), percentile=True))
        self.update_bins()

    def _on_metric(self, change):
        """Percentile bounds carry over; value bounds are re-mapped through their percentile."""
        if not self.percentile_w.value:
            pct = pt.values_to_percentiles(self.bins(), self.metric_values(change["old"]))
            self._set_bounds(pt.bin_edges(self.metric_values(change["new"]), pct, percentile=True))
        self.update_bins()

    # ---- drawing ------------------------------------------------------------------------------------

    @logged
    def update_bins(self):
        """Redraw [0, 0] with the current metric and bins, then every panel that depends on the bins."""
        self.log_w.clear_output()
        ax = self.axs[0, 0]
        ax.clear()
        edges = self.edges()
        pp.plot_bin_histogram(ax, self.metric_values(), edges, pp.bin_colors(len(edges)), self.metric_w.value,
                              n_runs=len(self.data.dataset))
        self.update_bin_bars()
        self.update_correlation()
        self.update_bin_traces()   # draws the canvas

    @logged
    def update_correlation(self):
        """Redraw the [1, 1] scatter of y vs x, colored by [0, 0] bin."""
        ax = self.axs[1, 1]
        ax.clear()
        df = pt.correlation_values(self.data, self.corr_x_w.value, self.corr_y_w.value, self.metric_w.value,
                                   self.bins(), self.percentile_w.value)
        pp.plot_correlation(ax, df, pp.bin_colors(len(self._rows)), self.corr_x_w.value, self.corr_y_w.value,
                            self.metric_w.value)
        self.fig.canvas.draw_idle()

    @logged
    def update_bin_bars(self):
        """Redraw the [0, 1] bars: mean ± std (and the individual values) of the bar metric per bin."""
        ax = self.axs[0, 1]
        ax.clear()
        metric = self.bar_metric_w.value
        labels, per_point = self.bin_labels()
        values = pt.bar_values(self.data, metric, labels, per_point)
        unit = pp.TABLE_ROW_UNITS[pt.BIN_METRICS[metric]]
        pp.plot_bin_bars(ax, values, pp.bin_colors(len(values["bin"].cat.categories)), metric,
                         self.metric_w.value, unit)
        self.fig.canvas.draw_idle()

    @logged
    def update_bin_traces(self):
        """Redraw the [1, 0] bin averages: behavior signal (top) and pO2 (bottom), mean ± std per bin."""
        ax_b, ax_p = self.ax_behavior, self.ax_po2
        ax_b.clear()
        ax_p.clear()
        labels, per_point = self.bin_labels()
        colors = pp.bin_colors(len(labels["bin"].cat.categories))
        unit = "rows" if per_point else "events"
        metric, behavior, po2 = self.metric_w.value, self.behavior_w.value, self.po2_w.value

        # behavior traces exist once per event; pO2 traces per event x point when the bins are per point
        behavior_avg = pt.group_average(self.responses("event", behavior), labels, "bin")
        po2_avg = pt.group_average(self.responses("point" if per_point else "event", po2), labels, "bin")
        pp.plot_bin_traces(ax_b, behavior_avg, colors, pp.signal_label(behavior), unit)
        pp.plot_bin_traces(ax_p, po2_avg, colors, pp.PO2_TITLES[po2], unit)

        ax_b.set_title(f"{pp.signal_label(behavior)} by {metric} bin", fontweight="bold", fontsize=10)
        ax_p.set_title(f"{pp.signal_label(po2)} by {metric} bin", fontweight="bold", fontsize=10)
        ax_b.tick_params(labelbottom=False)
        ax_p.set_xlim(*self.params.extract_window_s)
        ax_p.set_xlabel(f"{pp.event_xlabel(self.params)}, interpolated to {self.dt:g} s")
        self.fig.canvas.draw_idle()

    def display(self):
        """Show the controls above the figure (needs `%matplotlib widget`)."""
        bins_panel = widgets.VBox([
            widgets.HTML("<b>[0, 0] Bins</b>"),
            widgets.HBox([widgets.Label("Bins:"), self.add_w, self.remove_w, self.metric_w, self.percentile_w]),
            self.rows_w,
        ])
        margin = widgets.Layout(margin="0 0 0 24px")
        bars_panel = widgets.VBox([widgets.HTML("<b>[0, 1] Bin bars</b>"), self.bar_metric_w], layout=margin)
        traces_panel = widgets.VBox([widgets.HTML("<b>[1, 0] Bin averages</b>"), self.behavior_w, self.po2_w],
                                    layout=margin)
        corr_panel = widgets.VBox([widgets.HTML("<b>[1, 1] Correlation</b>"), self.corr_x_w, self.corr_y_w],
                                  layout=margin)
        display(widgets.VBox([widgets.HBox([bins_panel, bars_panel, traces_panel, corr_panel]), self.log_w,
                              self.fig.canvas]))
