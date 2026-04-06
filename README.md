# 384-Well Plate Sample Data — HTS Analysis Pipeline

Sample high-throughput screening (HTS) dataset with 5 simulated 384-well plates,
including QC analysis, normalization, and hit identification.

What I am looking to accomplish here? 
- Plot the QC metrics for all plates
- Short-list hits from single-point data for dose-responses to eventually get at IC50s

---

## Plate Layout

Each plate is a 16-row (A–P) × 24-column (1–24) grid with the following control layout:

```
     01  02  03  04 ... 20  21  22  23  24
A  [  L][ M][ E][ E]...[ E][ E][ E][ L][ H]
B  [  L][ M][ E][ E]...[ E][ E][ E][ L][ H]
C  [  L][ M][ E][ E]...[ E][ E][ E][ L][ H]
D  [  L][ M][ E][ E]...[ E][ E][ E][ L][ H]
E  [  L][ M][ E][ E]...[ E][ E][ E][ L][ H]
F  [  L][ M][ E][ E]...[ E][ E][ E][ L][ H]
G  [  L][ M][ E][ E]...[ E][ E][ E][ L][ H]
H  [  L][ M][ E][ E]...[ E][ E][ E][ L][ H]
I  [  L][ H][ E][ E]...[ E][ E][ E][ M][ H]
J  [  L][ H][ E][ E]...[ E][ E][ E][ M][ H]
K  [  L][ H][ E][ E]...[ E][ E][ E][ M][ H]
L  [  L][ H][ E][ E]...[ E][ E][ E][ M][ H]
M  [  L][ H][ E][ E]...[ E][ E][ E][ M][ H]
N  [  L][ H][ E][ E]...[ E][ E][ E][ M][ H]
O  [  L][ H][ E][ E]...[ E][ E][ E][ M][ H]
P  [  L][ H][ E][ E]...[ E][ E][ E][ M][ H]

L = Low control     (~100–500 RFU)      — Col 1 (all rows), Col 23 rows A–H
M = Medium control  (~5,000–15,000 RFU) — Col 2 rows A–H, Col 23 rows I–P
H = High control    (~25,000–35,000 RFU)— Col 2 rows I–P, Col 24 (all rows)
E = Experimental    (~500–30,000 RFU)   — Columns 3–22 (all rows)
```

---

## Steps

### Step 1: Generate Plate Data
Raw RFU values were simulated using `numpy.random.normal()` with the following parameters:

| Well Type | Mean (RFU) | Std Dev | Clipped Range |
|---|---|---|---|
| Low control | 300 | 80 | 50–600 |
| Medium control | 10,000 | 2,000 | 4,000–18,000 |
| High control | 30,000 | 2,000 | 22,000–38,000 |
| Experimental | 15,000 | 8,000 | 200–38,000 |

5 plates were generated with different random seeds. Raw data saved as:
- `sample_384_plate.csv` — Plate 1
- `sample_384_plate2.csv` through `sample_384_plate5.csv` — Plates 2–5

### Step 2: Calculate QC Metrics
Calculated per plate using all low control wells (Col 1 + Col 23 rows A–H, n=24)
and all high control wells (Col 24 + Col 2 rows I–P, n=24).

| Metric | Formula | Pass Threshold |
|---|---|---|
| Z' | `1 - (3σ_high + 3σ_low) / abs(μ_high - μ_low)` | > 0.5 |
| %CV Low | `(σ_low / μ_low) × 100` | < 10% |
| %CV High | `(σ_high / μ_high) × 100` | < 10% |
| Signal Window | `(μ_high - 3σ_high - μ_low - 3σ_low) / σ_low` | > 2 |

Results saved in `plate_QC_metrics.csv` (Plate 1) and `all_plates_QC_metrics.csv` (all 5 plates).

QC summary across plates:

| Plate | Z' | %CV Low | %CV High | Signal Window |
|---|---|---|---|---|
| Plate 1 | 0.795 ✓ | 17.07% ✗ | 6.60% ✓ | 465.08 ✓ |
| Plate 2 | 0.768 ✓ | 24.79% ✗ | 7.41% ✓ | 292.64 ✓ |
| Plate 3 | 0.792 ✓ | 23.83% ✗ | 6.63% ✓ | 333.85 ✓ |
| Plate 4 | 0.820 ✓ | 26.09% ✗ | 5.70% ✓ | 322.64 ✓ |
| Plate 5 | 0.776 ✓ | 28.04% ✗ | 7.12% ✓ | 279.40 ✓ |

Note: %CV Low consistently fails because low control RFU values (~300 RFU) are near
background noise, making them inherently more variable. This is typical in real HTS assays.

### Step 3: Normalize Experimental Data
Each well normalized to the plate's mean low and high control values:

```
% Activity = (Signal - Mean_Low) / (Mean_High - Mean_Low) × 100
```

- 0% = low control level (no activity)
- 100% = high control level (full activity)
- Values outside 0–100% reflect assay noise (expected in real screens)

Normalized plates saved as `sample_384_plate_normalized.csv` through `sample_384_plate5_normalized.csv`.

### Step 4: Trellis QC Plot
A trellis plot (`trellis_QC_plot.png`) was generated showing:
- **Top two rows:** Heatmaps of normalized % activity for all 5 plates (red = low, green = high)
- **Bottom row:** Bar charts comparing Z', %CV Low, %CV High, and Signal Window across all 5 plates,
  with red dashed lines indicating pass/fail thresholds

### Step 5: Hit Identification — >50% Activity
Experimental wells (columns 3–22) with normalized % activity > 50% were flagged as hits.
Control wells were excluded.

- Total hits: **794 across 5 plates** (~50% of experimental wells, as expected from random data)
- Saved in `hits_above_50pct.csv`

### Step 6: Hit Refinement — >75% Activity
Hits were further filtered to >75% activity, representing stronger candidates.

- Total hits: **253 across 5 plates**
- Saved in `hits_above_75pct.csv`

| Plate | Hits >75% |
|---|---|
| Plate 1 | 44 |
| Plate 2 | 54 |
| Plate 3 | 44 |
| Plate 4 | 58 |
| Plate 5 | 53 |

Top 10 strongest hits:

| Plate | Well | % Activity |
|---|---|---|
| Plate 5 | D14 | 129.24% |
| Plate 5 | F03 | 129.24% |
| Plate 5 | K19 | 129.24% |
| Plate 3 | I14 | 127.68% |
| Plate 2 | N05 | 127.58% |
| Plate 4 | C10 | 127.05% |
| Plate 1 | I18 | 125.97% |
| Plate 1 | H12 | 121.83% |
| Plate 2 | O15 | 121.54% |
| Plate 4 | L13 | 120.53% |

---

## Files

| File | Description |
|---|---|
| `sample_384_plate.csv` | Raw RFU data, Plate 1 |
| `sample_384_plate2.csv` – `sample_384_plate5.csv` | Raw RFU data, Plates 2–5 |
| `sample_384_plate_normalized.csv` | Normalized % activity, Plate 1 |
| `sample_384_plate2_normalized.csv` – `sample_384_plate5_normalized.csv` | Normalized % activity, Plates 2–5 |
| `plate_QC_metrics.csv` | QC metrics for Plate 1 with pass/fail |
| `all_plates_QC_metrics.csv` | QC metrics summary across all 5 plates |
| `trellis_QC_plot.png` | Trellis heatmap and QC bar charts |
| `hits_above_50pct.csv` | All experimental hits >50% activity |
| `hits_above_75pct.csv` | Refined hits >75% activity |
