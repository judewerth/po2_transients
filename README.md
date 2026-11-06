# po2_transients

**Event-triggered pO₂ transient analysis.** Does the oxygen partial pressure (pO₂) measured along a vessel change around spontaneous **pupil dilations**, **locomotion bouts** or **sleep-state transitions**?

This repository holds the full analysis: Python code that loads the synced `pO2_all_synced_runXX.mat` files, detects behavioral events, cuts pO₂ windows around each one, averages and bins them, and fits a linear mixed model. It also has Jupyter notebooks that run each analysis step by step, plus interactive viewers for exploring single runs and the pooled dataset.

> **Never used Python or a terminal?** You are the intended reader. Follow [Full setup from zero](#full-setup-from-zero) from top to bottom. Every command is explained, and nothing is assumed beyond some MATLAB experience. Setup takes about 30–45 minutes once, and after that you only need [Everyday use](#everyday-use).

---

## Contents

1. [What's in this repository](#whats-in-this-repository)
2. [Quick start (if you already know conda + git)](#quick-start-if-you-already-know-conda--git)
3. [Full setup from zero](#full-setup-from-zero)
   - [Step 0: The terminal in 5 minutes](#step-0-the-terminal-in-5-minutes)
   - [Step 1: Install Git](#step-1-install-git)
   - [Step 2: Install Miniconda (Python)](#step-2-install-miniconda-python)
   - [Step 3: Install VS Code](#step-3-install-vs-code)
   - [Step 4: Download this repository](#step-4-download-this-repository-git-clone)
   - [Step 5: Create the Python environment](#step-5-create-the-python-environment)
   - [Step 6: Add the data](#step-6-add-the-data)
   - [Step 7: Open the project and run your first notebook](#step-7-open-the-project-and-run-your-first-notebook)
4. [Which notebook do I use?](#which-notebook-do-i-use)
5. [How the code is organized: back end to front end](#how-the-code-is-organized-back-end-to-front-end)
6. [Changing the analysis](#changing-the-analysis) (and [Reproducibility](#reproducibility))
7. [Everyday use](#everyday-use)
8. [Saving your own changes with Git](#saving-your-own-changes-with-git)
9. [Troubleshooting](#troubleshooting)
10. [Glossary](#glossary)
11. [Using Claude Code / AI assistants with this repo](#using-claude-code--ai-assistants-with-this-repo)
12. [License](#license)

---

## What's in this repository

```
po2_transients/
├── README.md                  ← you are here
├── po2_transients.py          ← DATA layer: reads .mat files, builds tables, does every computation
├── po2_plots.py               ← FIGURE layer: draws those tables with matplotlib
├── po2_gui.py                 ← WIDGET layer: the two interactive viewers
├── presentation_figures.py    ← publication-style figures in one consistent style
├── notebooks/                 ← where you work: each notebook runs one analysis
│   ├── 01_pipeline_walkthrough.ipynb          what every function does, step by step, on one run
│   ├── 02_visualize_po2_gui.ipynb             interactive viewers (single run + all runs)
│   ├── 03_analyze_pupil_po2_transients.ipynb  pupil-dilation events → pO2
│   ├── 04_analyze_locomotion_po2_transients.ipynb  locomotion events → pO2
│   ├── 05_analyze_sleep_po2_transients.ipynb  sleep-state transitions → pO2
│   └── 06_presentation_figures.ipynb          builds the slide figures into figures/
├── data/                      ← put the session folders here (NOT uploaded to GitHub; see data/README.md)
├── figures/                   ← created when you save figures (NOT uploaded)
├── docs/
│   ├── architecture.md        ← how the code fits together, in depth
│   ├── terminal_and_git.md    ← every terminal / git command used here, explained
│   └── matlab_to_python.md    ← cheat sheet for MATLAB users
├── environment.yml            ← the "recipe" listing every Python package (and version) this code needs
├── pyproject.toml             ← makes the .py files importable from the notebooks
├── .env.example               ← template for your local settings (data location)
├── .gitignore                 ← tells git which files never to upload (data, outputs, .env)
├── .vscode/                   ← shared VS Code settings and recommended extensions
├── CLAUDE.md                  ← instructions for AI coding assistants
└── LICENSE
```

**MATLAB analogy.** The `.py` files are like a folder of `.m` function files. The notebooks are like MATLAB **Live Scripts** (`.mlx`): code, text and figures in one document, run section by section.

---

## Quick start (if you already know conda + git)

```bash
git clone https://github.com/judewerth/po2_transients.git
cd po2_transients
conda env create -f environment.yml        # creates env "po2" and pip-installs this repo in editable mode
conda activate po2
python -m ipykernel install --user --name po2 --display-name "Python (po2)"
# copy session folders into data/  (or: cp .env.example .env and set PO2_DATA_DIR)
code .                                     # open in VS Code, pick the "po2" kernel, run notebooks/01_...
```

---

## Full setup from zero

Instructions are written for **macOS**. Windows differences are marked with 🪟.

### Step 0: The terminal in 5 minutes

The **terminal** (also called the *command line* or *shell*) is a text window where you type commands instead of clicking. It is a lot like MATLAB's **Command Window**: you type a line, press Enter, and it runs. The terminal runs programs and moves around folders instead of MATLAB functions.

**Open it:** press `Cmd + Space`, type `Terminal`, press Enter.
🪟 *Windows:* after Step 2, use **Anaconda Prompt** from the Start menu. It is a terminal that already knows where conda is.

You'll see a **prompt** that ends in `%` or `$`, something like:

```
yourname@MacBook ~ %
```

The `~` is the **current folder**, and `~` is shorthand for your home folder (`/Users/yourname`). This is just like MATLAB's "Current Folder". Commands act on the current folder unless you give a full path.

Try these three. Type each one and press Enter:

| Command | What it does | MATLAB equivalent |
|---|---|---|
| `pwd` | **p**rint **w**orking **d**irectory: shows which folder you're in | `pwd` |
| `ls` | **l**i**s**t the files and folders here (🪟 `dir`) | `ls` / `dir` |
| `cd Documents` | **c**hange **d**irectory into `Documents` | `cd Documents` |
| `cd ..` | go up one folder (`..` means "the parent folder") | `cd ..` |

Useful habits:
- **Tab completes names.** Type `cd Doc` and press Tab, and it fills in `Documents/`. This avoids typos.
- **↑ (up arrow)** brings back the previous command, as in MATLAB.
- **Ctrl + C** stops a running command (MATLAB: Ctrl + C too).
- When this README shows a block like the one below, type or paste the line *without* any leading `%` or `$`. Lines starting with `#` are comments (explanations) and you don't need to type them.

```bash
# this is a comment, the line below is the command
pwd
```

More detail on every command: [docs/terminal_and_git.md](docs/terminal_and_git.md).

### Step 1: Install Git

**Git** is a version-control tool. It records snapshots ("commits") of a folder of code, so you can see what changed, undo mistakes, and sync with a copy on **GitHub** (a website that hosts git repositories). You need it to download this repository and get future updates.

**macOS:** in the terminal, run

```bash
xcode-select --install
```

A dialog pops up asking to install the "command line developer tools". Click **Install** and wait. (This installs Apple's developer tools, which include `git`. If it says they are already installed, you're done.)

🪟 *Windows:* download and run the installer from <https://git-scm.com/download/win>. The default options are fine.

**Check it worked:**

```bash
git --version
```

It should print something like `git version 2.39.0`. Any version is fine. *What happened:* you asked the `git` program to report its version. If the terminal answers, it knows where git is installed.

### Step 2: Install Miniconda (Python)

Your Mac may already have a Python, but **don't use it**. It belongs to the operating system. Instead we install **Miniconda**, which provides:

- **Python** itself, the language interpreter (like the MATLAB engine), and
- **conda**, a package and environment manager.

**Why environments?** An *environment* is a self-contained folder with one Python version and a specific set of packages (libraries). Each project gets its own, so updating a package for one project can't break another. Imagine being able to keep a separate MATLAB install with exactly the toolboxes (and toolbox versions) each project needs. This repo's environment is called **`po2`**, and its recipe is the file [`environment.yml`](environment.yml).

1. Go to <https://www.anaconda.com/download/success> and scroll to **Miniconda Installers**.
2. **macOS:** download the **Apple Silicon** `.pkg` if your Mac has an M1/M2/M3/M4 chip ( → About This Mac → "Chip: Apple M…"). Otherwise download the **Intel** one. Open it and click through the installer.
   🪟 *Windows:* download the Windows 64-bit `.exe`, run it, and choose "Just Me". Leave the other options at their defaults.
3. **Close the terminal window and open a new one.** Programs installed while a terminal is open aren't visible to it until it restarts.
4. Check it worked:

   ```bash
   conda --version
   ```

   It should print `conda 25.x.x` or similar. Your prompt should now start with **`(base)`**, which is the name of the environment currently *active*. `base` is conda's default environment, and we won't install anything into it.

   > If you get `conda: command not found`, run `~/miniconda3/bin/conda init zsh`, then close and reopen the terminal. (`conda init` adds a few lines to your shell's startup file, `~/.zshrc`, so every new terminal can find conda.)

### Step 3: Install VS Code

**Visual Studio Code** (VS Code) is a free code editor. It plays the role of the MATLAB desktop: file browser on the left, editor in the middle, terminal at the bottom, and it runs Jupyter notebooks.

1. Download from <https://code.visualstudio.com/> and install it. On macOS, drag it into Applications.
2. Open VS Code. Click the **Extensions** icon in the left bar (four squares), then search for and install:
   - **Python** (by Microsoft)
   - **Jupyter** (by Microsoft)
   (When you open this repository later, VS Code will also offer to install these, because they're listed in `.vscode/extensions.json`.)
3. *Optional but handy:* press `Cmd + Shift + P`, type `shell command`, and choose **Shell Command: Install 'code' command in PATH**. After this, typing `code .` in a terminal opens the current folder in VS Code. (🪟 Windows does this automatically.)

### Step 4: Download this repository (`git clone`)

Pick where the project should live, for example your `Documents/GitHub` folder:

```bash
mkdir -p ~/Documents/GitHub      # make the folder (-p: no error if it already exists; also makes parent folders)
cd ~/Documents/GitHub            # go into it
git clone https://github.com/judewerth/po2_transients.git
cd po2_transients                # go into the downloaded project
ls                               # you should see README.md, notebooks, po2_transients.py, ...
```

*What happened:* `git clone <url>` downloads the whole repository (every file plus its full history) into a new folder named `po2_transients`. The folder is a **clone**: a linked copy that can later fetch updates with `git pull`.

> 🪟 On Windows, in Anaconda Prompt, use `mkdir %USERPROFILE%\Documents\GitHub` and `cd %USERPROFILE%\Documents\GitHub` for the first two lines.

### Step 5: Create the Python environment

Still inside the `po2_transients` folder (check with `pwd`), run:

```bash
conda env create -f environment.yml
```

*What happened:* conda read the recipe file (`-f environment.yml`), created a new environment named `po2`, and downloaded Python 3.10 plus every package listed there (numpy, scipy, pandas, matplotlib, statsmodels, Jupyter support, …). Its last step ran `pip install -e .`, which **installs this repository's own `.py` files into the environment in "editable" mode**. That's why a notebook inside `notebooks/` can `import po2_transients` even though the file lives one folder up. "Editable" means Python reads the files from this folder directly, so if you edit `po2_transients.py`, the change is live and there's nothing to reinstall. (The MATLAB analogy is `addpath` that you never have to repeat.) This step takes a few minutes and prints a lot. A final line like `done` or no red `Error` means success.

Now **activate** the environment and register it with Jupyter:

```bash
conda activate po2
python -m ipykernel install --user --name po2 --display-name "Python (po2)"
```

*What happened:*
- `conda activate po2` switches this terminal to the `po2` environment, and your prompt now starts with **`(po2)`**. From now on, `python` in this terminal means the environment's Python with all the packages. You need to activate it again in every new terminal window.
- `python -m ipykernel install …` registers the environment as a Jupyter **kernel** named "Python (po2)". A kernel is the background Python process that runs notebook cells. Registering it makes it show up in VS Code's and Jupyter's kernel pickers.

**Check it worked:**

```bash
python -c "import po2_transients as pt; print('OK, data folder:', pt.DATA_DIR)"
```

(`python -c "…"` runs one line of Python without opening a file.) It should print `OK, data folder: …/po2_transients/data`.

### Step 6: Add the data

The recordings are **not on GitHub**. They are unpublished and too large, so ask Jude or check the lab drive for the synced session folders. Each session folder is named `YY-MM-DD_MOUSE` and contains `pO2_all_synced_runXX.mat` files. See [data/README.md](data/README.md) for the exact layout.

**Option A (simplest): copy them into `data/`.** Drag the session folders into the repo's `data/` folder in Finder (or VS Code's file panel). You should end up with:

```
po2_transients/data/26-10-05_B6_318/pO2_all_synced_run05.mat
po2_transients/data/26-10-05_B6_318/pO2_all_synced_run06.mat
...
```

**Option B: leave them where they are** (e.g. on a lab drive) and point the code at them with a `.env` file:

```bash
cp .env.example .env        # copy the template (cp = copy). On Windows: copy .env.example .env
```

Open `.env` in VS Code and fill in the folder that *contains* the session folders:

```
PO2_DATA_DIR=/Volumes/LabDrive/pO2/synced
```

*What is `.env`?* A small text file of `NAME=value` settings that are specific to **your** computer. `po2_transients.py` reads it when it is imported (via the `python-dotenv` package) to find the data. It is listed in `.gitignore`, so your personal paths are never uploaded. `.env.example` is the shared template that *is* uploaded.

### Step 7: Open the project and run your first notebook

```bash
code .          # "." means "this folder". Or: open VS Code → File → Open Folder… → po2_transients
```

1. If VS Code asks *"Do you trust the authors?"*, click **Yes**. If it offers to install recommended extensions, click **Install**.
2. In the left file panel, open **`notebooks/01_pipeline_walkthrough.ipynb`**.
3. Click **Select Kernel** (top-right) → **Jupyter Kernel…** (or *Python Environments…*) → **Python (po2)**.
4. Click into the first code cell and press **Shift + Enter**. The cell runs (a ✓ and a number appear next to it) and the cursor moves to the next cell. Keep pressing Shift + Enter and read along, or click **Run All** at the top.

You should see printed tables and figures of the example run's pupil, locomotion, sleep state and pO₂. If you do, **setup is complete** 🎉. If not, see [Troubleshooting](#troubleshooting).

> **Notebook basics.** A notebook is a list of **cells**. *Markdown* cells hold formatted text. *Code* cells hold Python, and their output (prints, tables, figures) appears right below them. Variables are shared between cells like the MATLAB workspace, so run cells top to bottom. **Restart** (top bar) clears every variable, like `clear all`. Notebooks in this repo are saved **without** outputs, so you see figures only after running them.

> *Prefer a browser?* With the environment active, run `jupyter lab` in the terminal. It opens Jupyter Lab in your web browser, where you pick the "Python (po2)" kernel. Press Ctrl + C in the terminal to stop it when you're done.

---

## Which notebook do I use?

| Notebook | Use it to… | Typical run time |
|---|---|---|
| [01_pipeline_walkthrough](notebooks/01_pipeline_walkthrough.ipynb) | **Learn the code.** Calls every function on one example run and shows what it returns and why. Start here. | ~1 min |
| [02_visualize_po2_gui](notebooks/02_visualize_po2_gui.ipynb) | **Explore interactively.** Pick runs, change thresholds, click events/points to exclude them, and choose bins for the pooled data. | interactive |
| [03_analyze_pupil_po2_transients](notebooks/03_analyze_pupil_po2_transients.ipynb) | Full analysis with **pupil-dilation** events: detection QC, pooled responses, classes, mixed model | a few min |
| [04_analyze_locomotion_po2_transients](notebooks/04_analyze_locomotion_po2_transients.ipynb) | Same analysis with **locomotion** events | a few min |
| [05_analyze_sleep_po2_transients](notebooks/05_analyze_sleep_po2_transients.ipynb) | Same analysis with **sleep-state transitions** | a few min |
| [06_presentation_figures](notebooks/06_presentation_figures.ipynb) | Rebuild the slide figures (PNG + SVG) into `figures/` | ~1 min |

The three analysis notebooks (03–05) have **the same structure**. Only the event source and its parameters differ:

1. **Configuration**: every number you might change, in one cell.
2. **Steps 1–5**: `pt.process_all_runs` loads, filters, detects events and builds the tables for every run.
3. **Per-run detection QC**: one figure per run showing which events were detected.
4. **Signal quality**: pO₂ point noise (CV) and the point QC cutoff.
5. **Pooled responses**: every event window and their average.
6. **Event classes**: events binned (e.g. into tertiles of size) and averaged per class, per point, plus an interpolation check.
7. **Linear mixed model**: does the pO₂ peak depend on event class, accounting for mouse/date/run/event/point?

---

## How the code is organized: back end to front end

The code is split into **layers**, and each layer only uses the ones below it:

```
             ┌────────────────────────────────────────────────────────────┐
FRONT END    │ notebooks/*.ipynb    set parameters, call functions, show   │  ← you work here
             └───────────────┬───────────────────────────┬────────────────┘
                             │                           │
             ┌───────────────▼──────────┐   ┌────────────▼───────────────┐
             │ po2_gui.py   (gui)       │   │ presentation_figures (pf)  │
             │ widgets + click handlers │   │ styled slide figures       │
             └───────┬─────────┬────────┘   └────────────┬───────────────┘
                     │         │                         │
             ┌───────│─────────▼─────────────────────────│───────────────┐
             │       │   po2_plots.py   (pp)              │               │
             │       │   draws tables onto matplotlib axes│               │
             └───────│─────────┬─────────────────────────│───────────────┘
                     │         │                         │
             ┌───────▼─────────▼─────────────────────────▼───────────────┐
BACK END     │ po2_transients.py  (pt)                                    │
             │ .mat → tables; filtering, detection, metrics, binning, LMM │
             │ never plots                                                │
             └───────────────────────────┬────────────────────────────────┘
                                         │ reads
                               data/<session>/pO2_all_synced_runXX.mat
```

**Why split it this way?** Each kind of change has one place to go. Changing *what is computed* (a filter, a metric) means editing `po2_transients.py`, and every notebook, viewer and figure picks it up. Changing *how a figure looks* means editing `po2_plots.py`. The notebooks stay short: they read like a methods section, with parameters followed by calls followed by results.

### What happens to one run (the 5 stages in `po2_transients.py`)

```
pO2_all_synced_run06.mat
   │  1 LOAD        load_synced_run          → raw behavior signal(s) + pO2 matrix (time × points)
   │  2 PREPROCESS  preprocess_event_source  → filtered signal; events are detected on one of them
   │                  pupil: moving median → low-pass → derivative
   │                  locom: |accel − mean| → Gaussian smoothing
   │                  sleep: clean short bouts → derivative (nonzero at transitions)
   │  3 DETECT      detect_events            → event times + per-event metrics (height, half-width, AUC)
   │  4 EXTRACT     extract_event_segments   → a fixed window (e.g. −30…+30 s) around each event, from every signal and pO2 point
   │  5 pO2 METRICS POINT_METRICS, PO2_EVENT_METRICS → per point (rms, cv, …), per event × point (event_peak, …)
   ▼
TransientData: five pandas tables
```

`process_all_runs(runs, params)` does this for every run and stacks the results into one object, `data`, holding **five tables** (pandas DataFrames, which behave much like MATLAB tables):

| Table | One row per | Holds |
|---|---|---|
| `data.dataset` | run | raw time vectors and signals |
| `data.po2_signals` | run × pO₂ point | whole-run pO₂ trace and its QC numbers (rms, cv, baseline) |
| `data.extract_signals` | run | filtered behavior signals and the detection threshold |
| `data.events` | event | event time, size metrics, and the behavior trace around it |
| `data.po2_events` | event × pO₂ point | pO₂ window around the event, normalized, and its peak |

Every table starts with the same **key columns** (`session, mouse_id, date, run`, plus `event` and/or `point`), so tables can be joined the way you would `join` MATLAB tables.

### Tracing a figure back to the code

Take the **"Event responses by event class"** figure in notebook 03:

1. The notebook cell calls `pt.classify_events(data.events, "event_auc", class_bins)`. This adds an `event_class` column to `data.events` (back end, a pure computation).
2. It then calls `pp.plot_class_averages(data, "event_class", params)` (figure layer).
3. Inside, `pp.plot_class_averages` asks the back end for numbers: `pt.event_responses` interpolates each event window onto a 1 s grid, and `pt.group_average` gives the mean ± SD per class.
4. `pp` draws those numbers with matplotlib (`ax.plot`, `ax.fill_between`), the same idea as MATLAB's `plot` and `fill`.

The interactive viewers work the same way, plus **callbacks**. When you type a new threshold, `po2_gui.py` runs a function that calls `pt.analyze_run(...)` with the new value and asks `pp` to redraw the affected panels. Nothing is computed inside the GUI itself.

**To see this yourself:** in VS Code, Cmd-click (🪟 Ctrl-click) any function name in a notebook, such as `pp.plot_class_averages`, to jump to its source. Every function has a docstring (the text in `"""…"""` under its name) that says what goes in and what comes out.

The full tour, including the "registries" that let you add a new event type or metric with a few lines, is in **[docs/architecture.md](docs/architecture.md)**.

---

## Changing the analysis

| I want to… | Do this |
|---|---|
| Analyze different runs or add a session | Put the folder in `data/` and edit the `runs = {...}` dict in the notebook's configuration cell |
| Detect more / fewer events | Lower / raise `height_thresh` in `params` and check the per-run QC figures |
| Change the window around events | `extract_window_s=(-30, 30)` (seconds before, after) |
| Change how events are binned | `class_column`, `class_bins`, `class_percentile`, or pick bins interactively in notebook 02 (`all_viewer.bins()`) |
| Exclude noisy pO₂ points | `max_cv` |
| Change the statistical model | `lmm_outcome`, `lmm_fixed`, `lmm_random`, `lmm_groups` |
| Change a filter, or add a metric or event type | Edit `po2_transients.py`. See "How to extend" at the top of that file and [docs/architecture.md](docs/architecture.md#8-extending-the-pipeline-registries) |
| Change a figure's look | Edit the function in `po2_plots.py` (or `presentation_figures.py` for slides) |

After editing a `.py` file you don't need to restart. The notebooks use `%autoreload 2`, which re-imports changed files before each cell. Re-run the cells that use the changed function.

> **Known issues** are documented at the top of `po2_transients.py` (e.g. a ≤ 1-sample alignment offset in event windows). They are deliberately left unchanged so results stay reproducible. Read them before interpreting fine timing.

### Reproducibility

- `environment.yml` pins **exact** package versions, and installs the scientific packages with pip rather than conda on purpose. Several mixed-model fits sit near a parameter boundary (statsmodels warns *"MLE may be on the boundary"*), and there the optimizer's answer shifts slightly with the linear-algebra library underneath numpy. The pip packages reproduce the original results on macOS. On Windows or Linux the last digits of a mixed-model fit may differ slightly. Event detection, counts and averages are unaffected.
- With the committed settings and all 11 runs, you should get **42** pupil events (790 event × point rows), **86** locomotion events (1550), and **48** sleep-state transitions (870). If you get these numbers, your setup matches.

---

## Everyday use

Each time you sit down to work:

```bash
cd ~/Documents/GitHub/po2_transients   # go to the project
conda activate po2                     # use the project's Python (prompt shows "(po2)")
git pull                               # download any updates from GitHub (see below)
code .                                 # open VS Code
```

In VS Code you can also open a terminal *inside* the editor: **Terminal → New Terminal** (or Ctrl + \`). It starts in the project folder.

**If `environment.yml` changed in an update** (git pull lists it), update your environment with:

```bash
conda env update -f environment.yml --prune    # install new packages, remove ones no longer listed
```

---

## Saving your own changes with Git

Git tracks changes to the *code and notebooks*. Data, figures and `.env` are ignored on purpose (see `.gitignore`).

```bash
git status              # what changed since the last commit? (red = changed, not yet staged)
git diff                # show the exact line changes
```

**Before committing a notebook, clear its outputs** (VS Code: the **…** menu at the top of the notebook → *Clear All Outputs*, then save). Saved outputs contain figures of unpublished data and make the files huge.

To keep your own version of the analysis without affecting the shared one, work on a **branch**:

```bash
git switch -c my-analysis          # create a branch called "my-analysis" and move onto it
# ...edit notebooks / code...
git add notebooks/03_analyze_pupil_po2_transients.ipynb   # "stage" a file: mark it for the next snapshot
git commit -m "Pupil analysis with 20 s window"           # take the snapshot, with a message saying why
```

A commit is saved **only on your computer**. To share it, either ask Jude to add you as a collaborator (then `git push -u origin my-analysis` and open a Pull Request on GitHub), or **fork** the repo on GitHub (the *Fork* button makes your own online copy). Every command is explained in [docs/terminal_and_git.md](docs/terminal_and_git.md).

---

## Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `conda: command not found` | Terminal opened before Miniconda was installed. Close and reopen it. Still failing: run `~/miniconda3/bin/conda init zsh` and reopen. |
| `FileNotFoundError: File …/data/26-10-05_B6_318/pO2_all_synced_run05.mat does not exist` | The data isn't where the code looks. The printed path shows where that is. Copy the session folders into `data/`, or fix `PO2_DATA_DIR` in `.env`, then **restart the kernel** (`.env` is read at import). |
| `ModuleNotFoundError: No module named 'po2_transients'` (or `numpy`, …) | The notebook is using the wrong kernel. Select **Python (po2)** (top-right). If it's not listed, run Step 5's `ipykernel install` line again. |
| The `po2` kernel isn't in VS Code's list | Command palette (`Cmd + Shift + P`) → *Developer: Reload Window*. Or choose *Select Another Kernel → Python Environments → po2*. |
| GUI figures show as text (`Canvas(…)`) or don't respond to clicks | Run the GUI notebook's first cell (`%matplotlib widget`) **before** any figure, and use the `po2` kernel (it includes `ipympl`). Restart the kernel and run from the top. |
| `CondaValueError: prefix already exists` when creating the env | The `po2` environment already exists. Update it with `conda env update -f environment.yml --prune`, or start over with `conda env remove -n po2` and then create it again. |
| `LinAlgError: Singular matrix` in the mixed-model cell of notebook 05 (or §10 of notebook 01) | Known and expected with the committed settings: sleep-state events have only two `event_auc` values, so the lowest percentile class is empty and the model can't be fitted. The note above that cell explains how to bin by transition type instead. |
| Lots of pink `RuntimeWarning` / `ConvergenceWarning` text | Warnings only. The cell still finished. Mixed models often warn when a variance component is near zero. |
| A cell runs forever | `process_all_runs` reads every `.mat` file, so give it a few minutes (longer from a network drive). To stop it: the ■ button (*Interrupt*). |
| `git pull` says "Your local changes would be overwritten" | You edited a file that was also updated online. Save your version on a branch first (see above), or discard your edits to that file with `git restore <file>`. |

Still stuck? Copy the **last lines** of the error (the "Traceback") and send them to Jude, or paste them into an AI assistant along with the cell that failed.

---

## Glossary

| Term | Meaning |
|---|---|
| **Terminal / shell** | Text interface to your computer. You type commands, it runs them (macOS default shell: `zsh`). |
| **Path** | A file's address, e.g. `/Users/you/Documents/GitHub/po2_transients/data`. `~` = home folder, `.` = this folder, `..` = parent folder. |
| **Python** | The programming language (and the program, the *interpreter*, that runs it). |
| **Package / library** | Reusable code you import, like a MATLAB toolbox (numpy, pandas, matplotlib, …). |
| **Module** | One `.py` file of Python code. `import po2_transients as pt` loads it under the short name `pt`. |
| **conda** | Installs Python and packages, and manages **environments**. |
| **Environment** | An isolated set of Python + packages for one project (here: `po2`). *Activate* it to use it. |
| **pip** | Python's own package installer. conda uses it for this repo (`pip install -e .`). |
| **Jupyter notebook (`.ipynb`)** | Document of text + code cells + outputs, like a MATLAB Live Script. |
| **Kernel** | The Python process that runs a notebook's cells and holds its variables. |
| **numpy array** | N-dimensional numeric array, like a MATLAB matrix. **Indexing starts at 0**, not 1. |
| **pandas DataFrame** | A table with named columns, like a MATLAB `table`. |
| **dict** | Key → value lookup (`runs["26-10-05_B6_318"]`), like a MATLAB `struct` or `containers.Map`. |
| **dataclass** | A struct-like object with named fields and defaults (`pt.EventParams`). |
| **Git** | Version control: records snapshots (*commits*) of the code. |
| **Repository (repo)** | A project folder tracked by git, here `po2_transients`. |
| **Clone** | Your local, linked copy of a repo on GitHub. |
| **Commit** | A saved snapshot of changes, with a message. |
| **Branch** | An independent line of commits, so you can experiment without touching `main`. |
| **Push / pull** | Upload your commits to GitHub / download others' commits from GitHub. |
| **Fork** | Your own copy of someone else's repo on GitHub. |
| **`.gitignore`** | List of files git should never track (data, outputs, `.env`). |
| **`.env`** | Local, private `NAME=value` settings (here: where the data lives). Never committed. |
| **Editable install (`pip install -e .`)** | Makes this folder's `.py` files importable from anywhere, using the files in place. |

The MATLAB ↔ Python syntax cheat sheet is in [docs/matlab_to_python.md](docs/matlab_to_python.md).

---

## Using Claude Code / AI assistants with this repo

[CLAUDE.md](CLAUDE.md) gives AI coding assistants (Claude Code, and others that read such files) the context they need: the layer rules, where parameters live, which known issues not to "fix", and that data must never be committed. If you use one, open it in the repo folder and ask questions like *"explain what `pt.detect_events` does"* or *"add a per-event metric for peak prominence"*.

---

## License

[MIT](LICENSE): free to use and adapt, with attribution. The **data** are not part of this repository or its license.
