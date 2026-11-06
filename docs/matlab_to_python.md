# MATLAB → Python cheat sheet

Python syntax you'll meet in this repo, next to its MATLAB equivalent. Python libraries used here: **numpy** (`np`, arrays), **pandas** (`pd`, tables), **matplotlib** (`plt`, plots), **scipy** (signal processing, `.mat` reading), **statsmodels** (mixed models).

## The five biggest differences

1. **Indexing starts at 0.** `x[0]` is the first element, and `x[-1]` is the last (MATLAB: `x(1)`, `x(end)`).
2. **Square brackets index, parentheses call.** `x[3]` indexes, `f(3)` calls a function. MATLAB uses `()` for both.
3. **Slices exclude the end.** `x[0:3]` gives elements 0, 1, 2 (three elements). MATLAB's `x(1:3)` includes 3.
4. **Indentation is syntax.** Blocks (`if`, `for`, `def`) end when the indentation ends. There is no `end`.
5. **Tools are imported explicitly.** `import numpy as np`, then `np.mean(x)`. Nothing is on the "path" by default.

## Basics

| Python | MATLAB | Notes |
|---|---|---|
| `# comment` | `% comment` | |
| `x = 5` | `x = 5;` | Python never echoes assignments. Use `print(x)` or put `x` alone on the last line of a cell |
| `print(f"n = {n}")` | `fprintf('n = %d\n', n)` | f-string: `{...}` is filled in from variables |
| `True, False, None` | `true, false, []` | `None` = "no value" |
| `and, or, not` | `&&, \|\|, ~` | For arrays, use `&, \|, ~` (elementwise) |
| `x == 1`, `x != 1` | `x == 1`, `x ~= 1` | |
| `x ** 2` | `x.^2` | `*` on numpy arrays is elementwise. Matrix multiply is `@` |
| `"text"` or `'text'` | `"text"` / `'text'` | |
| `len(x)` | `numel(x)` / `length(x)` | |
| `help(fn)` or `fn?` (in Jupyter) | `help fn` | Shift + Tab in Jupyter, or hover in VS Code |

## Containers

| Python | MATLAB | Notes |
|---|---|---|
| `[1, 2, 3]` (list) | `{1, 2, 3}` / `[1 2 3]` | Lists can mix types and grow: `lst.append(4)` |
| `(-30, 30)` (tuple) | `[-30 30]` | Like a list, but can't be changed |
| `{"a": 1, "b": 2}` (dict) | `struct('a',1,'b',2)` / `containers.Map` | Access with `d["a"]`. Keys can be strings with any characters |
| `d.keys()`, `d.items()` | `fieldnames(s)` | `for k, v in d.items():` loops over key/value pairs |
| `params.height_thresh` (dataclass field) | `params.height_thresh` (struct field) | `pt.EventParams` is a dataclass: a struct with fixed fields and defaults |

## numpy arrays

| Python | MATLAB |
|---|---|
| `np.array([[1, 2], [3, 4]])` | `[1 2; 3 4]` |
| `A.shape` | `size(A)` |
| `A[0, 1]` | `A(1, 2)` |
| `A[:, 0]` (first column) | `A(:, 1)` |
| `A[1:, :]` (all but first row) | `A(2:end, :)` |
| `np.arange(0, 10, 2)` | `0:2:8` |
| `np.linspace(0, 1, 5)` | `linspace(0, 1, 5)` |
| `np.zeros((3, 4))`, `np.full((3, 4), np.nan)` | `zeros(3, 4)`, `nan(3, 4)` |
| `np.nanmean(A, axis=0)` | `mean(A, 1, 'omitnan')` |
| `A[A > 0]` | `A(A > 0)` |
| `np.where(A > 0)` | `find(A > 0)` |
| `np.concatenate([a, b])` | `[a, b]` / `cat` |

`axis=0` runs **down the rows** (one result per column), like MATLAB's dimension 1.

Arrays here are often **(n_time, n_points)**: `po2[:, 3]` is the whole trace of point 3, and `po2[100, :]` is every point at sample 100.

## pandas DataFrames (≈ MATLAB tables)

| Python | MATLAB |
|---|---|
| `df.head()` | `head(T)` |
| `df.columns` | `T.Properties.VariableNames` |
| `df["event_auc"]` | `T.event_auc` |
| `df[df["run"] == 6]` | `T(T.run == 6, :)` |
| `df.loc[3, "event_auc"]` (by label), `df.iloc[3]` (by position) | `T(4, :)` |
| `df["new"] = values` | `T.new = values` |
| `a.merge(b, on=keys)` | `innerjoin(a, b, 'Keys', keys)` |
| `df.groupby("event_class")["event_peak"].mean()` | `groupsummary(T, 'event_class', 'mean', 'event_peak')` |
| `df["event_class"].value_counts()` | `groupcounts(T, 'event_class')` |
| `len(df)` | `height(T)` |

In this repo some cells hold whole arrays (traces), shown as `array(600,)`. Get one with `df.loc[0, "event_trace"]`.

## Control flow and functions

```python
# Python                                   % MATLAB
for run in [5, 6, 8]:                      for run = [5 6 8]
    print(run)                                 disp(run)
                                           end
if n > 10:                                 if n > 10
    msg = "many"                               msg = "many";
elif n > 0:                                elseif n > 0
    msg = "some"                               msg = "some";
else:                                      else
    msg = "none"                               msg = "none";
                                           end
def rms(x):                                function y = rms(x)
    """Root mean square."""                    % Root mean square.
    return np.sqrt(np.mean(x ** 2))            y = sqrt(mean(x.^2));
                                           end
```

**Keyword arguments.** `pt.EventParams(event_method="pupil", min_distance_s=30)` names each argument, so order doesn't matter and any left out keep their defaults. MATLAB's closest equivalent is name-value pairs: `fn('event_method', 'pupil')`.

**List comprehension.** `[x ** 2 for x in values if x > 0]` builds a list in one line, like `arrayfun` plus logical indexing.

## Plotting

| Python (matplotlib) | MATLAB |
|---|---|
| `fig, ax = plt.subplots()` | `figure; ax = gca;` |
| `fig, axs = plt.subplots(2, 1, sharex=True)` | `tiledlayout(2, 1)` |
| `ax.plot(t, x, color="k", lw=2)` | `plot(t, x, 'k', 'LineWidth', 2)` |
| `ax.fill_between(t, lo, hi, alpha=0.3)` | `fill([t fliplr(t)], [lo fliplr(hi)], ...)` |
| `ax.set_xlabel("Time (s)")` | `xlabel('Time (s)')` |
| `ax.axvline(0)` | `xline(0)` |
| `plt.show()` | `drawnow` |
| `fig.savefig("f.png", dpi=300)` | `exportgraphics(gcf, 'f.png', 'Resolution', 300)` |

## Reading .mat files

`po2_transients.loadmat_nested_dict(path)` reads both MATLAB formats and returns nested dicts: `raw["VS"]["classification"]["frame_state"]` is MATLAB's `VS.classification.frame_state`. Remember that MATLAB's 1-based indices stored *in* the data (e.g. frame numbers) are still 1-based after loading.

## Jupyter-specific

| Thing | Meaning |
|---|---|
| `%something` | A notebook "magic" command, not Python (`%autoreload 2`, `%matplotlib widget`) |
| Last line of a cell is a bare variable | Its value is displayed nicely (tables render as HTML) |
| `display(df)` | Show a value nicely from anywhere in a cell |
| Kernel → Restart | Clear every variable (`clear all`) |
