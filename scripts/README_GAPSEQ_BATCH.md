# gapseq Batch Pipeline

NCBI GCF accession → gap-filled SBML GEM → PeptoMatch gem_cache

## Files

| File | Where to run | Purpose |
|------|--------------|---------|
| `data/strain_registry.csv` | (metadata) | 10 priority strains + Gram + medium |
| `scripts/download_genomes.sh` | **WSL/Ubuntu** | Bulk FASTA download from NCBI |
| `scripts/gapseq_batch.sh` | **WSL/Ubuntu** | find → find-transport → draft → fill |
| `scripts/copy_gems.sh` | **WSL/Ubuntu** | Copy xml → `/mnt/d/folder1/peptomatch/outputs/gem_cache/` |

## Prerequisites (one-time, already done)

- gapseq installed at `~/gapseq` (Ubuntu)
- cobrar R package installed
- `gapseq find/draft/fill -h` all work
- LP manually verified (pipeline sanity check)

## Workflow

### 1) Download genomes (2 min)

```bash
bash /mnt/d/folder1/peptomatch/scripts/download_genomes.sh
```

Fetches all 10 NCBI accessions into `~/genomes/{abbrev}.fna`.
Already-downloaded strains are skipped (idempotent).

Verify:
```bash
ls -lh ~/genomes/*.fna
```

### 2) Batch gapseq (~3 hours for 9 strains with 2 parallel)

LP is already done, so skip it:

```bash
bash /mnt/d/folder1/peptomatch/scripts/gapseq_batch.sh --skip LP --parallel 2
```

Options:
- `--parallel N` — run N strains concurrently (needs `parallel` installed: `sudo apt install parallel`)
- `--priority 1` — only priority-1 strains (default)
- `--only LR,LA` — run only these
- `--skip LP` — exclude these
- `--force` — re-run even if `xml` already exists

Resume: re-run the same command — it auto-skips strains with existing `{abbrev}/{abbrev}.xml`.

Logs: `~/gapseq_work/logs/{abbrev}.log`

### 3) Copy to PeptoMatch (10 sec)

```bash
bash /mnt/d/folder1/peptomatch/scripts/copy_gems.sh
```

Copies `~/gapseq_work/{abbrev}/{abbrev}.xml` → `/mnt/d/folder1/peptomatch/outputs/gem_cache/{abbrev}.xml`

### 4) Verify in Windows

```cmd
cd D:\folder1\peptomatch
python -c "from peptomatch.fba_simulator import FBASimulator; import os; [print(f'{s}: mu =', FBASimulator(f'outputs/gem_cache/{s}').predict_growth()) for s in sorted(os.listdir('outputs/gem_cache')) if s.endswith('.xml')]"
```

Expected μ ranges:
- Lactobacillus (LP, LR, LA, LC, LPC, LS): **0.2 – 0.6 h⁻¹**
- Enterococcus / Streptococcus (EF, STT): **0.3 – 0.7**
- Bacillus (BS, BC): **0.5 – 1.0**

## Troubleshooting

| Symptom | Fix |
|---------|-----|
| `download: cannot resolve FTP path` | NCBI changed ASM suffix — open the URL in a browser, copy new FTP path |
| `fill: no medium file found` | gapseq missing `ALLmed.csv` — try `ls ~/gapseq/dat/media/`, adjust `FALLBACK_MEDIA` |
| `μ = 0` or very small | Strain is fastidious; try `--force` with a richer medium (manually edit gapseq_batch.sh `FALLBACK_MEDIA`) |
| `parallel: command not found` | `sudo apt install parallel` or use `--parallel 1` |
| Some strains fail with exit 10/11 | `find`/`find-transport` issue — check `~/gapseq_work/logs/{abbrev}.log` |

## After completion

With `outputs/gem_cache/` populated (10 xml files), PeptoMatch is FBA-ready:

- `gateway.py` `/api/fba/*` endpoints will find models
- `fba_simulator.py` can predict μ for any strain
- Next: wire `gem_manager.py` to auto-discover via `strain_registry.csv`

ML models (`peptomatch_ridge_max_od.joblib` etc.) are separate — train those once experimental data accumulates (`python -m peptomatch.ml_train`).
