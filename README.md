# 384-Well Plate HTS Analysis Pipeline

A fully automated high-throughput screening (HTS) analysis pipeline using 5 simulated 384-well plates.
The pipeline runs as a sequence of Jupyter notebooks, taking raw plate RFU data through to IC50 calculation.

What do I want to accomplish in this project? 
- Use assay-aware computational analysis to identify promising hits
- Demonstrate use of automated pipeline to anaysis any set of raw data generated


---

## Plate Layout

Each plate is a 16-row (A–P) × 24-column (1–24) grid:

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

L = Low control     (~100–500 RFU)       Col 1 all rows, Col 23 rows A–H
M = Medium control  (~5,000–15,000 RFU)  Col 2 rows A–H, Col 23 rows I–P
H = High control    (~25,000–35,000 RFU) Col 2 rows I–P, Col 24 all rows
E = Experimental    (~500–30,000 RFU)    Columns 3–22 all rows
```

---

## Pipeline

Run the notebooks in order. Each notebook reads from the previous one's output.

| Notebook | Description | Output |
|---|---|---|
| `01_plate_metadata.ipynb` | Define plate layout and control positions | `data/plate_metadata.json` |
| `02_load_raw_data.ipynb` | Load all raw plate CSVs into a unified long-format dataset | `data/processed/all_plates_raw.csv` |
| `03_qc_and_normalization.ipynb` | Calculate Z', %CV, Signal Window; normalize to controls | `data/processed/qc_metrics.csv`, `all_plates_normalized.csv` |
| `04_plot_qc_metrics.ipynb` | Trellis heatmaps and QC metric bar charts across all plates | `data/results/trellis_QC_plot.png` |
| `05_hit_calling.ipynb` | Flag experimental wells >50% and >75% activity | `data/results/hits_50pct.csv`, `hits_75pct.csv` |
| `06_dose_response_setup.ipynb` | Cherry-pick top hits and generate 8-point dose-response data | `data/results/dose_response_data.csv` |
| `07_ic50_calculation.ipynb` | Fit 4PL curve per compound and calculate IC50 | `data/results/ic50_results.csv`, `ic50_curves.png` |

---

## QC Metrics

| Metric | Formula | Pass Threshold |
|---|---|---|
| Z' | `1 - (3σ_high + 3σ_low) / abs(μ_high - μ_low)` | > 0.5 |
| %CV Low | `(σ_low / μ_low) × 100` | < 10% |
| %CV High | `(σ_high / μ_high) × 100` | < 10% |
| Signal Window | `(μ_high - 3σ_high - μ_low - 3σ_low) / σ_low` | > 2 |

| Plate | Z' | %CV Low | %CV High | Signal Window |
|---|---|---|---|---|
| Plate 1 | 0.795 ✓ | 17.07% ✗ | 6.60% ✓ | 465.08 ✓ |
| Plate 2 | 0.768 ✓ | 24.79% ✗ | 7.41% ✓ | 292.64 ✓ |
| Plate 3 | 0.792 ✓ | 23.83% ✗ | 6.63% ✓ | 333.85 ✓ |
| Plate 4 | 0.820 ✓ | 26.09% ✗ | 5.70% ✓ | 322.64 ✓ |
| Plate 5 | 0.776 ✓ | 28.04% ✗ | 7.12% ✓ | 279.40 ✓ |

Note: %CV Low consistently fails because low control RFU values (~300 RFU) sit close to
background noise, making them inherently variable. This is typical in real HTS assays.

---

## Normalization

```
% Activity = (Signal - Mean_Low) / (Mean_High - Mean_Low) × 100
```

- 0% = low control (baseline)
- 100% = high control (full effect)
- Values outside 0–100% reflect assay noise

---

## Hit Calling

- **>50% activity:** 794 hits across 5 plates
- **>75% activity:** 253 hits across 5 plates (used for dose-response follow-up)

Top 10 hits by % activity:

| Plate | Well | % Activity |
|---|---|---|
| Plate 5 | D14, F03, K19 | 129.2% |
| Plate 3 | I14 | 127.7% |
| Plate 2 | N05 | 127.6% |
| Plate 4 | C10 | 127.1% |
| Plate 1 | I18 | 126.0% |
| Plate 1 | H12 | 121.8% |
| Plate 2 | O15 | 121.5% |
| Plate 4 | L13 | 120.5% |

---

## IC50 Calculation

Top 20 hits by % activity were taken forward for dose-response (8 concentrations: 0.03–100 µM).
A 4-parameter logistic (4PL) curve was fitted per compound using `scipy.optimize.curve_fit`:

```
f(x) = A + (D - A) / (1 + (C/x)^B)
```

Results saved in `data/results/ic50_results.csv`. See also the companion repo:
[4-PL-Curve-fit-in-Python](https://github.com/sraashi/4-PL-Curve-fit-in-Python---numpy-and-scipy)

---

## Repo Structure

```
├── 01_plate_metadata.ipynb
├── 02_load_raw_data.ipynb
├── 03_qc_and_normalization.ipynb
├── 04_plot_qc_metrics.ipynb
├── 05_hit_calling.ipynb
├── 06_dose_response_setup.ipynb
├── 07_ic50_calculation.ipynb
├── data/
│   ├── plate_metadata.json
│   ├── raw/                  ← raw plate CSVs (5 plates)
│   ├── processed/            ← normalized data and QC metrics
│   └── results/              ← plots, hit lists, IC50 results
└── README.md
```
