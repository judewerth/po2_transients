# Terminal and Git: every command used in this repo, explained

A reference for the commands in the [README](../README.md): what each one does, what its parts mean, and what happens behind the scenes. Commands are for macOS (zsh). Windows differences are marked 🪟.

**Reading a command.** Commands follow the pattern `program  subcommand  --options  arguments`:

```
conda   env create   -f environment.yml
  │        │          │      └── argument: the file to use
  │        │          └── option (flag): "-f" = "use this file"
  │        └── subcommand: what conda should do
  └── the program being run
```

Short options start with one dash (`-f`), long ones with two (`--file`). They usually mean the same thing. Every program has built-in help: `conda --help`, `git help clone`.

---

## 1. Moving around and looking at files

| Command | What it does | Notes |
|---|---|---|
| `pwd` | Print the folder you're in | MATLAB: `pwd` |
| `ls` | List files in the current folder | `ls -la` also shows hidden files (names starting with `.`, like `.env`, `.gitignore`) and details. 🪟 `dir` |
| `cd <folder>` | Change into a folder | `cd ..` = up one, `cd ~` = home, `cd -` = back to the previous folder. Paths with spaces need quotes: `cd "My Folder"` |
| `mkdir -p <folder>` | Make a folder | `-p`: also make any missing parent folders, and don't complain if it exists |
| `cp <from> <to>` | Copy a file | `cp -r` copies a folder. 🪟 `copy` |
| `mv <from> <to>` | Move or rename | 🪟 `move` |
| `rm <file>` | Delete a file, **permanently**. There is no Trash | Be careful. `rm -r` deletes folders |
| `open .` | Open the current folder in Finder | 🪟 `explorer .` |
| `code .` | Open the current folder in VS Code | needs the `code` command, see README Step 3 |
| `cat <file>` | Print a text file's contents | |
| `clear` | Clear the screen | Ctrl + L also works |

**Paths.** An *absolute* path starts at the top of the disk (`/Users/you/Documents/...`). A *relative* path starts from the current folder (`notebooks/01_pipeline_walkthrough.ipynb`). `~` is your home folder, `.` is "here", `..` is "one up".

**Keyboard.** Tab completes names (press it twice to list the options). ↑/↓ scroll through previous commands. Ctrl + C stops the running command. Ctrl + A / Ctrl + E jump to the start / end of the line.

---

## 2. Conda: Python and environments

| Command | What it does |
|---|---|
| `conda --version` | Check that conda is installed and on the PATH |
| `conda init zsh` | Adds a few lines to `~/.zshrc` (the script every new terminal runs at start-up) so that `conda` and `conda activate` work. Needed once. 🪟 the installer does this for Anaconda Prompt. |
| `conda env create -f environment.yml` | Reads the recipe, creates the environment it names (`po2`), and installs every listed package. Behind the scenes it makes a folder like `~/miniconda3/envs/po2/` with its own `python` and `site-packages/` (where libraries live). |
| `conda activate po2` | Use that environment in **this terminal**. It changes `PATH` (the list of folders the shell searches for programs) so that `python`, `pip` and `jupyter` refer to `~/miniconda3/envs/po2/bin/...`. The prompt shows `(po2)`. |
| `conda deactivate` | Go back to the previous environment |
| `conda env list` | List environments. `*` marks the active one |
| `conda list` | List the packages installed in the active environment, with versions |
| `conda env update -f environment.yml --prune` | Bring the environment in line with an updated recipe. `--prune` removes packages no longer listed |
| `conda env remove -n po2` | Delete the environment (your code and data are untouched). Re-create it with `conda env create` |

**Why not just `pip install` into my system Python?** Packages would be shared by every project, and one project's upgrade could break another's. Environments isolate them, and `environment.yml` lets anyone rebuild the exact same set.

**What the `pip:` section of `environment.yml` does.** `pip` is Python's own installer. conda installs Python and pip, then hands the `pip:` list to pip: the scientific packages at exact versions, and `-e .`:

- `.` means "the current folder", which pip installs as a package using `pyproject.toml`.
- `-e` (*editable*) means: don't copy the files, point at them. Your edits to `po2_transients.py` are used immediately.

---

## 3. Running Python and Jupyter from the terminal

| Command | What it does |
|---|---|
| `python` | Start an interactive Python prompt (`>>>`), like the MATLAB command window. `exit()` or Ctrl + D to leave |
| `python -c "import po2_transients as pt; print(pt.DATA_DIR)"` | Run one line of Python and exit (`-c` = "command") |
| `python presentation_figures.py` | Run a Python file top to bottom, like `run script.m`. This one rebuilds every slide figure into `figures/` |
| `python -m ipykernel install --user --name po2 --display-name "Python (po2)"` | `-m` runs a module as a program. This registers the active environment as a Jupyter kernel called "Python (po2)" (a small JSON file in `~/Library/Jupyter/kernels/po2/` saying which `python` to start) |
| `jupyter lab` | Start Jupyter Lab, a notebook editor in your web browser. The terminal must stay open while you use it, and Ctrl + C stops it |
| `jupyter nbconvert --clear-output --inplace notebooks/*.ipynb` | Remove all outputs from every notebook (do this before committing). `*` is a wildcard: every file ending in `.ipynb` |
| `jupyter nbconvert --to notebook --execute <nb> --output-dir <dir>` | Run a notebook top to bottom without opening it, saving the executed copy to `<dir>` |

---

## 4. Git: version control

**Mental model.** A repository is a folder plus a hidden `.git/` subfolder that stores every saved snapshot (**commit**). Your files live in three places:

```
working folder  ──git add──▶  staging area  ──git commit──▶  history (.git/)  ──git push──▶  GitHub
   (your edits)               (next snapshot)                 (local commits)               (shared)
                                                                     ◀──git pull──  GitHub
```

| Command | What it does |
|---|---|
| `git --version` | Check that git is installed |
| `git clone <url>` | Download a repository from GitHub into a new folder, with its full history. The online copy becomes the **remote** named `origin` |
| `git status` | **The most useful command.** It shows your branch, which files changed (red = not staged, green = staged) and whether you're behind or ahead of GitHub |
| `git diff` | Show line-by-line what changed in files you haven't staged (`git diff --staged` for staged ones). `q` quits the viewer |
| `git pull` | Fetch new commits from GitHub and merge them into your current branch. Do this before starting work |
| `git add <file>` | Stage a file: include its current state in the next commit. `git add .` stages everything changed (`.gitignore` still applies) |
| `git restore <file>` | **Discard** your uncommitted edits to a file (back to the last commit). This can't be undone |
| `git commit -m "message"` | Save the staged changes as a snapshot with a message. Write *why* you made the change |
| `git log --oneline` | List commits, newest first. `q` quits |
| `git switch -c <name>` | Create a new **branch** and move onto it. A branch is an independent line of commits, so experiments don't touch `main` |
| `git switch <name>` | Move to an existing branch (your files change to match it) |
| `git branch` | List branches. `*` marks the current one |
| `git push -u origin <branch>` | Upload your branch's commits to GitHub. `-u` remembers the link, so later a plain `git push` is enough. Needs write access to the repo |

**First-time git setup** (once per computer, before your first commit). Git labels every commit with a name and email:

```bash
git config --global user.name "Your Name"
git config --global user.email "you@example.com"     # use your GitHub email (or its no-reply address)
```

**What `.gitignore` does.** It lists patterns git should never track. Here that's `data/*` (recordings), `figures/`, `.env`, and Python/Jupyter clutter. `git status` won't even show them, and `git add .` skips them. This is what keeps unpublished data off a public GitHub.

**Typical day**

```bash
git pull                                    # get updates
# ... work ...
git status                                  # what did I change?
git add notebooks/04_analyze_locomotion_po2_transients.ipynb
git commit -m "Locomotion: try 2 s smoothing"
git push                                    # if you have access and want to share
```

**Merge conflicts.** If you and someone else changed the same lines, `git pull` stops and marks the file with `<<<<<<<`, `=======` and `>>>>>>>` sections. Edit the file to keep what you want, delete the markers, then run `git add <file>` and `git commit`. Notebooks are JSON files and are hard to merge by hand. Clearing outputs before committing (and working on your own branch) avoids most conflicts.

---

## 5. Publishing this repository to GitHub (for the repo owner)

Done once, by Jude:

```bash
# 1. On github.com: "+" → New repository → name "po2_transients", Public, do NOT add a README/license
# 2. In this folder:
git remote add origin https://github.com/judewerth/po2_transients.git   # link the local repo to the online one
git push -u origin main                                                  # upload the main branch
```

Then add lab members as collaborators: on GitHub, open the repo → *Settings → Collaborators*.
