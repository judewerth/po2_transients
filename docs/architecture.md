# Architecture: how the back end produces what you see

This page explains how the pieces of the code fit together: which file does what, how data flows from a `.mat` file to a figure, and what happens when you click something in an interactive viewer. Read the [README](../README.md) first. To see every function's actual output on real data, run [notebooks/01_pipeline_walkthrough.ipynb](../notebooks/01_pipeline_walkthrough.ipynb) alongside this page.

**Contents**
1. [The big picture](#1-the-big-picture)
2. [How a notebook finds the code and the data](#2-how-a-notebook-finds-the-code-and-the-data)
3. [The pipeline for one run (stages 1–5)](#3-the-pipeline-for-one-run-stages-15)
4. [The five tables](#4-the-five-tables)
5. [After the tables: averaging, binning, statistics](#5-after-the-tables-averaging-binning-statistics)
6. [The figure layer (`po2_plots.py`)](#6-the-figure-layer-po2_plotspy)
7. [The interactive layer (`po2_gui.py`): from a click to a redraw](#7-the-interactive-layer-po2_guipy-from-a-click-to-a-redraw)
8. [Extending the pipeline: registries](#8-extending-the-pipeline-registries)
9. [Map of `po2_transients.py`](#9-map-of-po2_transientspy)

---

## 1. The big picture

```
                    YOU
                     │  edit parameters, run cells, click widgets
                     ▼
┌──────────────────────────────────────────────────────────────────────┐
│ notebooks/*.ipynb                         (front end)                │
│   params = pt.EventParams(...)          ← all the numbers you choose │
│   data = pt.process_all_runs(runs, params)                           │
│   pp.plot_class_averages(data, ...)                                  │
└──────┬─────────────────────────┬────────────────────────┬────────────┘
       │ calls                   │ calls                  │ calls
       ▼                         ▼                        ▼
 po2_gui.py (gui)  ──calls──▶ po2_plots.py (pp) ──calls──▶ po2_transients.py (pt)
 widgets, state,              draws tables onto            reads .mat files, builds
 click callbacks              matplotlib axes              tables, computes everything
       │                                                   ▲
       └──────────────────────calls────────────────────────┘
                                                           │ reads
                                         data/<session>/pO2_all_synced_runXX.mat
```

There are three rules, and they explain most of the design:

1. **Arrows point one way only.** `pt` imports nothing from `pp` or `gui`, and `pp` imports nothing from `gui`. You can therefore use `pt` in a script with no plotting at all, and a plotting bug can never change a number.
2. **`pt` never plots and `pp` never computes.** If a figure shows a mean, `pp` got it by calling a `pt` function (e.g. `pt.group_average`). If you want to know *how* a number was computed, look in `po2_transients.py`.
3. **Parameters live in one object.** `pt.EventParams` is a dataclass (a struct with named fields and defaults). The notebook builds one and passes it down. Nothing deeper in the code has its own hidden settings.

**MATLAB analogy.** Picture three folders of `.m` functions: `+data`, `+plots` and `+gui`. Your Live Script calls them, `+plots` functions call `+data` functions, and `+data` never calls the other two.

---

## 2. How a notebook finds the code and the data

**The code.** Notebooks live in `notebooks/`, but the `.py` files are in the folder above. Normally Python only finds modules next to the running file or in installed packages. The environment's last install step, `pip install -e .`, uses [`pyproject.toml`](../pyproject.toml) to register the four modules (`po2_transients`, `po2_plots`, `po2_gui`, `presentation_figures`) as installed packages that point back at this folder. So `import po2_transients as pt` works from any folder, and because the install is *editable* (`-e`), Python reads your actual files. Edit `po2_transients.py`, and with `%autoreload 2` the next cell uses the new version.

**The data.** At import, `po2_transients.py` runs this logic (section "0. Where the data lives"):

```python
REPO_ROOT = folder that contains po2_transients.py
load .env from REPO_ROOT if it exists          # python-dotenv copies KEY=value lines into the environment
DATA_DIR = REPO_ROOT / (PO2_DATA_DIR or "data")
```

Every file read builds its path as `DATA_DIR / <session folder> / <file>`. Each notebook prints `pt.DATA_DIR` in its first code cell. `presentation_figures.py` does the same for `FIG_DIR` (`PO2_FIG_DIR`, default `figures/`).

---

## 3. The pipeline for one run (stages 1–5)

These stages carry the same numbers in the section headers of `po2_transients.py` and in the walkthrough notebook.

```
pO2_all_synced_run06.mat (+ sleepVikram_run06.mat for sleep state)
  │
  │ STAGE 1  load_synced_run(folder, run, raw_names)
  │   reads the .mat (scipy for old files, h5py for v7.3), picks variables via RAW_SIGNALS
  ▼
sig = {t_po2, fs_po2, t_pupil, pupil, fs_pupil, ...}      po2 = array (n_time × n_points)
  │
  │ STAGE 2  preprocess_event_source(sig, params)  →  EVENT_METHODS[params.event_method].preprocess
  │   pupil:        moving median (movmedian_len_s) → zero-phase low-pass (lpf_cutoff_hz) → first difference
  │   locom:        |accel − mean| → Gaussian smoothing (locom_sigma_s)
  │   sleep_state:  drop bouts shorter than min_distance_s, fill gaps → first difference
  ▼
ext = {t_<s>, <s>} for every filtered signal s;  one of them (detect_on) is the detection signal
  │
  │ STAGE 3  detect_events(t, x, fs, height_thresh, min_distance_s, rel_to_rms, neg_peaks)
  │   find_event_peaks: scipy.signal.find_peaks above the threshold, ≥ min_distance_s apart
  │   → DetectionContext (signal, peaks, lobe bounds, find_peaks properties)
  │   → every EVENT_METRICS function applied to it: event_height, event_halfwidth, event_auc
  ▼
events: event_time, event_bounds, metrics (one value per event)
  │
  │ STAGE 4  extract_event_segments(t, signal, event_times, extract_window_s)
  │   same-length window around each event (NaN where it runs off the recording),
  │   applied to every filtered signal AND every pO2 point, each at its own sample rate
  ▼
segments (n_events × n_window)
  │
  │ STAGE 5  pO2 metrics
  │   POINT_METRICS on the whole-run pO2 of each point: rms, cv, baseline, nan_frac
  │   normalize_event_segments: (window − point rms) / point rms
  │   PO2_EVENT_METRICS on each normalized window: event_peak, peak_abs_height, event_peak_time, ...
  ▼
TransientData (five tables, next section)
```

Three wrapper functions group the stages:

| Function | Stages | Why it exists |
|---|---|---|
| `load_run(folder, run, params)` → `LoadedRun` | 1–2 | The slow part (file reading + filtering). The GUI **caches** its result. |
| `analyze_run(loaded, params)` → `TransientData` | 3–5 | The fast part. It re-runs whenever a threshold changes. |
| `process_run` / `process_all_runs(runs, params)` | 1–5 | One run / every run in the `runs` dict, stacked into one `TransientData`. This is what the analysis notebooks call. |

---

## 4. The five tables

`TransientData` is a small container holding five pandas DataFrames. A DataFrame is like a MATLAB `table`. One difference: here a single cell can hold a whole array (a trace), which pandas shows as `array(600,)`.

| Table | One row per | Key columns | Main other columns |
|---|---|---|---|
| `dataset` | run | `session, mouse_id, date, run` | `t_po2, fs_po2`, raw `t_<src>, <src>, fs_<src>` |
| `po2_signals` | run × pO₂ point | + `point` | `po2, abs_po2, norm_po2` (traces), `rms, cv, baseline, nan_frac` |
| `extract_signals` | run | run keys | filtered `t_<s>, <s>`, `event_threshold` |
| `events` | event | + `event` | `event_time, event_bounds, event_height, event_halfwidth, event_auc`, per signal: `t_<s>_event, <s>_event_trace` |
| `po2_events` | event × pO₂ point | + `event, point` | `t_event, event_trace, event_trace_norm, event_peak, peak_abs_height, event_peak_time, ...` |

**Keys and joins.** Every table starts with the run keys, and `pt.TABLE_KEYS[table]` lists the columns that make a row unique. To attach event information to each event × point row, merge on the event keys, the same idea as MATLAB's `join`:

```python
df = data.po2_events.merge(data.events[pt.EVENT_KEYS + ["event_auc"]], on=pt.EVENT_KEYS)
```

The analysis notebooks do exactly this before fitting the mixed model.

**Point indices only mean something within a run.** Point 3 of run 5 and point 3 of run 6 are not the same vessel location.

**Two pO₂ normalizations.** `norm_po2 = po2 / baseline` (1 = the point's mean) is used for whole-run traces. `event_trace_norm = (trace − rms) / rms` (0 = the point's RMS) is used for event windows and `event_peak`.

---

## 5. After the tables: averaging, binning, statistics

These functions read the tables (sections 10–12 of `po2_transients.py`):

- **Event averages.** Runs have different sample rates, so windows can't be averaged sample by sample. `interpolate_event` puts one window on a regular `dt` grid. `event_responses` does that for every row of a signal. `group_average` gives the mean ± SD per group (e.g. per event class).
- **Binning.** `classify_events(df, column, bins, percentile)` labels each row with a class (e.g. tertiles of `event_auc`). `bin_labels`, `bar_values` and `correlation_values` are what the all-runs viewer uses to line up metrics from different tables.
- **Statistics.** `point_qc` (per-point CV for the QC cutoff), `event_point_table` (one row per event × point, ready for modeling), `fit_mixed_model` (builds a formula and fits `statsmodels` MixedLM, similar to MATLAB's `fitlme`), `extract_fixed_effects`, `extract_variance_decomposition`.

---

## 6. The figure layer (`po2_plots.py`)

Every function takes tables (or the result of a `pt` function) and draws. There are two kinds:

- `plot_*(ax, ...)` draws **into an axes you give it** and returns nothing. The GUI uses these to redraw single panels.
- Whole-figure functions (e.g. `plot_class_averages(data, ...)`, `plot_run_events(...)`) create a figure when you don't pass axes, and **return it**. The notebooks use these.

**Example trace: `pp.plot_class_averages(data, "event_class", params, dt=1.0)`**

```
pp.plot_class_averages
 ├─ for each signal (behavior signals, then pO2):
 │    pt.event_responses(data, signal, params.extract_window_s, dt)   # every window on a 1 s grid
 │    pt.group_average(responses, labels=data.events["event_class"])  # mean ± SD per class
 └─ ax.plot(mean) + ax.fill_between(mean ± sd)                        # matplotlib, like MATLAB plot/fill
```

Colors, labels and axis names are defined once at the top of `po2_plots.py` (`COLORS`, `signal_label`, …), so every figure matches.

---

## 7. The interactive layer (`po2_gui.py`): from a click to a redraw

The viewers are built from **ipywidgets** (dropdowns, text boxes, buttons drawn in the notebook) and **ipympl** (makes a matplotlib figure live in the notebook, so it can receive mouse clicks; enabled by `%matplotlib widget`).

Interactive code is **event-driven**. Instead of running top to bottom once, it *registers callbacks*: functions that run when something happens. MATLAB GUIs work the same way (`ButtonPushedFcn`, `ValueChangedFcn`).

```python
self.run_w.observe(self._on_run, names="value")                     # run dropdown changed → _on_run
self.fig.canvas.mpl_connect("button_press_event", self._on_click)   # click on the figure → _on_click
```

**What happens when you type a new threshold in `RunGUI` and press Enter:**

```
threshold text box value changes
 └─ _on_change  (callback)
     └─ update()
         ├─ update_detection()
         │    params = current params with the new threshold          (dataclasses.replace)
         │    loaded = loaded_run()            ← from the cache {(folder, run, method): LoadedRun};
         │                                       the .mat file is NOT read again
         │    self.run_data = pt.analyze_run(loaded, params)   ← stages 3–5 recomputed
         │    pp.plot_detection_panel(ax[0,0], self.run_data, params)
         └─ update_points()
              pp.plot_point_bars / plot_po2_traces / plot_po2_events  → panels [1,1], [1,0], [0,1]
```

**Clicking a pO₂ bar** only changes `self.selected` (which points are shown) and calls `update_points()`, with no recomputation. **Clicking an event line** toggles it in `self.excluded` and redraws the same panels. **Choosing a new run or event method** is the only action that calls `pt.load_run` (it reads the file), and only on a cache miss.

`AllRunsGUI` works the same way on an already-computed `data`. Changing bins calls `pt.bin_labels` / `pt.group_average` / `pt.correlation_values` and redraws the affected panels. `all_viewer.classify()` returns the GUI's bins as a column you can use downstream.

---

## 8. Extending the pipeline: registries

A *registry* is a dict that maps a name to a function. The tables, `BIN_METRICS`, the plots and the GUI dropdowns all loop over these dicts, so **adding one entry is enough** for a new metric to show up everywhere.

| Registry | Add an entry to… | Function signature |
|---|---|---|
| `RAW_SIGNALS` | read a new variable from the `.mat` | `name: (time_var, value_var)` (not a function) |
| `EVENT_METHODS` | detect events on a new signal | `EventMethod(raw=[...], signals=[...], detect_on=..., preprocess=fn, fs_from=..., describe=...)` |
| `EVENT_METRICS` | add a per-event scalar | `fn(ctx: DetectionContext) -> array (n_events,)` |
| `POINT_METRICS` | add a per-point scalar | `fn(po2 (n_time, n_points)) -> array (n_points,)` |
| `PO2_VARIANTS` | add a whole-run pO₂ trace variant | `fn(po2, point_metrics) -> array (n_time, n_points)` |
| `PO2_EVENT_METRICS` | add a per-event × point scalar | `fn(t_rel (n_window,), norm (n_events, n_window)) -> array (n_events,)` |

**Example: add peak prominence as an event metric.** Put this in `po2_transients.py` next to the other event metrics:

```python
def event_prominence(ctx):
    """Prominence of each event peak (from scipy.signal.find_peaks)."""
    return ctx.props["prominences"]

EVENT_METRICS["event_prominence"] = event_prominence
```

`data.events` now has an `event_prominence` column, and you can bin by it in the notebooks and in the all-runs viewer. Section 11 of the walkthrough notebook runs this example live.

---

## 9. Map of `po2_transients.py`

| Section | Contents |
|---|---|
| Module docstring | Data-flow diagram, table schemas, how to extend, **known issues A–C** (left unfixed on purpose) |
| 0 | `REPO_ROOT`, `DATA_DIR` (reads `.env`) |
| 1 | Key column lists (`RUN_KEYS`, `EVENT_KEYS`, …, `TABLE_KEYS`) |
| 2 | `EventParams` |
| 3 | `TransientData` |
| 4 | Stage 1: `.mat` readers, `RAW_SIGNALS`, `load_synced_run` |
| 5 | Stage 2: filters, `EVENT_METHODS`, `preprocess_event_source` |
| 6 | Stage 3: `DetectionContext`, `find_event_peaks`, `EVENT_METRICS`, `detect_events` |
| 7 | Stage 4: `extract_event_segments` |
| 8 | Stage 5: `POINT_METRICS`, `PO2_VARIANTS`, `PO2_EVENT_METRICS` |
| 9 | Drivers: `LoadedRun`, `load_run`, `analyze_run`, `process_run`, `process_all_runs` |
| 10 | Event traces and averages: `interpolate_event`, `event_responses`, `group_average` |
| 11 | Binning: `BIN_METRICS`, `classify_events`, `bin_labels`, `bar_values`, `correlation_values` |
| 12 | Statistics: `point_qc`, `event_point_table`, `fit_mixed_model`, `extract_*` |

In VS Code, **Cmd/Ctrl-click** any function name to jump to it. The **Outline** view (Explorer sidebar → Outline) lists every function in the open file.
