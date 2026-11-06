# data/

This folder is where the code looks for recordings by default. **Its contents are never uploaded to GitHub.** `.gitignore` excludes everything here except this README, because the data are unpublished and large (~250 MB for three sessions).

## Getting the data

Ask Jude Werth for the synced session folders, or copy them from the lab drive. Put each **session folder** directly inside `data/`.

## Expected layout

```
data/
├── 26-10-05_B6_318/                 ← one folder per session, named YY-MM-DD_MOUSEID
│   ├── pO2_all_synced_run05.mat     ← required: pO2 + pupil + locomotion for run 5 (two-digit run number)
│   ├── sleepVikram_run05.mat        ← needed only for sleep-state events (event_method="sleep_state")
│   ├── pO2_all_synced_run06.mat
│   ├── sleepVikram_run06.mat
│   └── ...                          ← other files (protocol notes, slides) are fine and ignored by the code
├── 26-10-07_B6_322/
└── 26-10-21_B6_322/
```

- **Folder name**: must match `YY-MM-DD_MOUSEID` (e.g. `26-10-05_B6_318`). The code splits it into `date` and `mouse_id` (`pt.parse_session`), so a folder with another name raises an error.
- **Run files**: `pO2_all_synced_runXX.mat` must contain `tImaging`, `pO2OverTime`, `tPupilShort`, `pupilShort`, `tAccelShort` and `AccelShort`. `sleepVikram_runXX.mat` must contain `VS.classification.t_img` and `VS.classification.frame_state`. Both MATLAB formats work: v7.3 (HDF5) files are read with `h5py`, older ones with `scipy.io.loadmat`. The variable mapping is `RAW_SIGNALS` in `po2_transients.py`.
- The notebooks choose which runs to analyze in their `runs = {...}` configuration dict. Files that aren't listed there are ignored.

## Keeping the data somewhere else

You don't have to copy the data here. Point the code at another folder, one that *contains* the session folders, by creating a `.env` file at the repo root:

```bash
cp .env.example .env
```

and setting, for example:

```
PO2_DATA_DIR=/Volumes/LabDrive/pO2/synced
```

Restart the notebook kernel afterwards. Every notebook prints `Reading session folders from: ...` in its first code cell so you can check the path.
