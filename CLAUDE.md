# CLAUDE.md

Guidance for Claude Code (and other AI coding assistants) working in this repository.

## What this repo is

An event-triggered pO₂ transient analysis. Synced MATLAB recordings (`data/<YY-MM-DD_MOUSE>/pO2_all_synced_runXX.mat`, plus `sleepVikram_runXX.mat` for sleep state) are turned into pandas tables. Behavioral events (pupil dilations, locomotion bouts, sleep-state transitions) are detected, pO₂ windows are cut around each event, and the windows are averaged, binned and modeled with a linear mixed model. Users work in the Jupyter notebooks in `notebooks/`.

**Primary users are lab members who know MATLAB but are new to Python, the terminal and git.** When you explain or suggest terminal commands, say what each one does. Prefer MATLAB analogies, and don't assume knowledge of conda, git or Jupyter internals.

## Layout and layer rules

```
po2_transients.py  (pt)  DATA layer: .mat → tables, and every computation. Never imports matplotlib or plots.
po2_plots.py       (pp)  FIGURE layer: draws pt tables. Does no analysis; it asks pt for every number.
po2_gui.py         (gui) WIDGET layer: ipywidgets/ipympl state and callbacks; calls pt for data, pp for drawing.
presentation_figures.py (pf)  styled slide figures; uses pt (and its own drawing code).
notebooks/0X_*.ipynb     set parameters, call pt/pp/gui, show results. No analysis logic in notebooks.
```

Dependencies point one way only: `pt ← pp ← gui`. Keep it that way: a change to *what is computed* goes in `pt`, and a change to *how it looks* goes in `pp`.

- **Pipeline stages** (numbered the same in the `po2_transients.py` docstring and section headers): 1 load → 2 preprocess → 3 detect → 4 extract → 5 pO₂ metrics → `TransientData` with the five tables `dataset, po2_signals, extract_signals, events, po2_events`. `load_run` = stages 1–2, `analyze_run` = 3–5, and `process_all_runs` stacks every run.
- **Parameters** live in one `pt.EventParams` dataclass, built in each notebook's configuration cell. The notebooks are the only place parameter values are set.
- **Extension points are registries** (`RAW_SIGNALS`, `EVENT_METHODS`, `EVENT_METRICS`, `POINT_METRICS`, `PO2_VARIANTS`, `PO2_EVENT_METRICS`). Add a new metric or event type by registering a function. Don't special-case it in plotting or GUI code. See "How to extend" at the top of `po2_transients.py`.
- **Known issues A–C** in the `po2_transients.py` docstring (and the `NOTE` comments) are deliberately left unfixed so results stay reproducible. Don't "fix" them unless explicitly asked, and if asked, say that results will change.

## Data and paths

- `pt.DATA_DIR` = `PO2_DATA_DIR` from `.env` / the environment, defaulting to `<repo>/data`. Every file read goes through it (`synced_run_path`, and the sleep-state file in `load_synced_run`). `session` keys in the tables stay the bare folder name.
- `pf.FIG_DIR` = `PO2_FIG_DIR`, defaulting to `<repo>/figures`.
- **Never commit data** (`data/*` except `data/README.md`), `figures/`, or `.env`. The repo is public and the data are unpublished.

## Environment

```bash
conda env create -f environment.yml      # env "po2"; pip-installs exact versions + this repo in editable mode (-e .)
conda activate po2
python -m ipykernel install --user --name po2 --display-name "Python (po2)"
```

The scientific packages are pinned **via pip** on purpose. The mixed-model fits are near a boundary, and conda-forge (OpenBLAS) builds shift them slightly compared with the pip wheels the results were produced with. Don't move them to conda packages or loosen the pins without re-checking the model output.

Because of `pip install -e .`, notebooks import the modules from any folder, and `%autoreload 2` picks up edits to the `.py` files without a restart.

## Conventions

- **Notebooks are committed with outputs cleared**: `jupyter nbconvert --clear-output --inplace notebooks/*.ipynb`. Outputs contain figures of unpublished data and bloat the diffs.
- Notebook kernelspec is `po2` ("Python (po2)").
- Every public function has a docstring whose first line says what it returns. `01_pipeline_walkthrough.ipynb` builds its function index from these, so keep them.
- Keep the notebooks' beginner annotations ("Before you run this", "Behind the scenes", "What to look for") accurate when you change the pipeline.
- Comments explain *why*. Keep the existing style of the module docstrings (tables of stages, tables, registries).

## Checking a change

There is no test suite. To verify, execute the notebooks headlessly and check for errors and unchanged event counts:

```bash
cd notebooks
jupyter nbconvert --to notebook --execute --allow-errors --ExecutePreprocessor.kernel_name=po2 --output-dir /tmp/po2_check 0*.ipynb
```

Reference counts with the committed parameters (all 11 runs): pupil 42 events / 790 event×point rows, locomotion 86 / 1550, sleep state 48 / 870. Locomotion LMM log-likelihood: −698.6356.

**Expected, pre-existing failure:** the mixed-model cell raises `LinAlgError: Singular matrix` in `05_analyze_sleep_po2_transients.ipynb` and in §10 of `01_pipeline_walkthrough.ipynb` (which is configured for `sleep_state`). Sleep-state `event_auc` has only two values, so the lowest percentile class is empty. It is documented in both notebooks. Don't treat it as a regression, and don't change the analysis to "fix" it unless asked.
