# DEFUSE-Implement

Implementation package for [**DEFUSE: Double-Robust Efficiency-Fused Semi-Supervised Estimation**](https://arxiv.org/abs/2405.18722)

---

## Project Structure

This repository is organized into **three main components**:

1. **MCAR Scenario**  
2. **MAR Scenario**  
3. **Real Data**  

---

## Simulation Studies

The simulation section covers two types of missingness mechanisms:

### **1. MCAR and MAR Scenarios**

This section includes **six simulation settings** as described in our main paper:

| Setting | Description |
|----------|-------------|
| **MCAR I** | Single BM, continuous Y, correctly specified control variate. |
| **MCAR II** | Single BM, binary Y, misspecified control variate. |
| **MCAR III** | Multiple BM, continuous Y, correctly specified control variate. |
| **MCAR IV** | Single and multiple BM with continuous outcome. |
| **MCAR V** | Multiple BM with continuous outcome. |
| **MCAR VI** | Auxiliary covariates with continuous outcome. |

All simulations can be executed using the main script:  
`Implement.R`

---

### **File Organization**

Each simulation folder (`Source/`) contains **four core R scripts**:

| File | Description |
|------|--------------|
| **`computing_function.R`** | Auxiliary computational utilities including Hessian matrix evaluation, score function computation, DML learning, asymptotic variance estimation, and polynomial basis construction. |
| **`data_generation.R`** | Data-generation routines for all simulation designs. |
| **`DEFUSE.R`** | Implements the **preliminary estimator**. |
| **`main_function_xxx.R`** | Implements the **enhanced two-step estimator**, including `DEFUSE_LM` and `DEFUSE`. |

---

## Real Data Analysis

The **Real Data** component consists of two application studies:

1. **MIMIC 3**  
2. **NACC**

Each dataset folder shares a similar internal structure:

- `Implement.R` — main execution script (analogous to simulation implementation).  
- `Source/` — contains main estimation functions.  
- `data/` — includes pseudo (de-identified) data, column-wise shuffled while keeping the same variable names and missing structure.
- `Check_misalign.R` — verifies source misalignment.  
- `misalign/` — includes supporting functions for misalignment checks.

---

## Computation time

Below summarizes the **average computation time** for one replication under different settings, based on execution on a **MacBook Pro (M3 Max, 36 GB RAM, macOS 15.6)**:

| Component | Setting / Sample Size | Approx. Time per Replication |
|------------|----------------------|------------------------------|
| **Simulation (MCAR)** | *n = 500* | ~0.5 minute |
| **Simulation (MAR)** | *n = 500* | ~0.4 minute |
|  | *n = 1000* | ~1-2 minutes |
|  | *n = 2000* | ~5-8 minutes |
| **Real Data: MIMIC-III** | check paper | ~2–3 minutes |
| **Real Data: NACC** | check paper | ~15–20 minutes |


## Mapping Between Paper and Repository

The following table summarizes the correspondence between the **tables and figures in the paper** and their **source scripts / settings** in this repository:

| Paper Component | Source of Result | Description |
|------------------|------------------|--------------|
| **Table 2** | `Scenario/MCAR II` and `Scenario/MCAR III` | Generated from single- and multiple-BM simulations with binary and continuous outcomes. |
| **Table 3** | `Scenario/MAR IV` | Based on single BM simulations under MAR with continuous outcome. |
| **Table 4** | `Scenario/MAR V` | Based on the Auxiliary single BM simulation under MAR with continuous outcome. |
| **Figure 1** | `Scenario/MAR V` | Visualization of the selection performance under multiple-BM MAR setting. |
| **Table 5** | `Real Data/NACC` | Results derived from the NACC dataset analysis. |
| **Table 6** | `Real Data/MIMIC 3` | Results derived from the MIMIC-III dataset analysis. |
| **Table S2** | `Scenario/MCAR I` | Supplementary results fo the single BM continuous outcome simulation with correctly specified control variate. |
| **Figure S1 & Figure S2** | `Scenario/MCAR I` | Supplementary figures visualizing performance under the single BM continuous Y setting. |
| **Table S3** | `Scenario/MCAR II` | Supplementary results corresponding to single-BM binary Y scenario. |
| **Table S4** | `Scenario/MAR IV` | Based on multiple BM simulations under MAR with continuous outcome. |
| **Table S5** | `Real Data/NACC/Check Misalign` and `Real Data/MIMIC 3/Check Misalign` | Results from the misalignment diagnostic analysis in real data section. |
| **Table S6 & Table S7** | `Real Data/NACC` and `Real Data/MIMIC 3` | Coefficients estimated from real data, with SE's and P values |



---
