<p align="center">
  <h1 align="center">Semiconductor AI Portfolio</h1>
  <p align="center">
    <strong>AI-driven spectral analysis for semiconductor metrology — from PhD research to production-ready tools</strong>
  </p>
  <p align="center">
    <a href="https://github.com/kyb8801/semiconductor-ai-portfolio/actions/workflows/ci.yml">
      <img src="https://github.com/kyb8801/semiconductor-ai-portfolio/actions/workflows/ci.yml/badge.svg" alt="CI">
    </a>
    <img src="https://img.shields.io/badge/Python-3.9%20%7C%203.10%20%7C%203.11-blue?logo=python&logoColor=white" alt="Python">
    <img src="https://img.shields.io/badge/License-MIT-green.svg" alt="License">
  </p>
</p>

---

## About

This portfolio demonstrates AI/ML applications in semiconductor optical metrology, built on **8+ years of hands-on experience** with AFM, SEM, Raman, NSOM, and TERS. Each project uses real research data from doctoral studies in optics and nanoscience.

**Author:** Yongbeom Kim (Ph.D. in Optics/Nanoscience)
**Focus:** AI for semiconductor metrology — Optics PhD × AI × uncertainty quantification

---

## Projects

### 01 — MoSe2 PL Spectrum Analysis
Automated peak detection and FWHM calculation for photoluminescence spectra of MoSe2 monolayers.
- Savitzky-Golay smoothing + `scipy.find_peaks`
- Data: WITec confocal Raman/PL measurements (personal PhD data)

### 02 — MoSe2 NSOM Defect Mapping
Hyperspectral near-field optical microscopy data analysis with unsupervised clustering.
- K-Means clustering on 30x30 NSOM spatial maps
- Spectral feature extraction + defect region identification

### 03 — TMD NSOM Comparison (WSe2 / WS2)
Cross-material comparison of transition metal dichalcogenide optical properties.
- WSe2 bilayer CVD vs WS2 bilayer analysis
- Peak position mapping and strain analysis

### 04 — Metrology Toolkit
Simulation tools for semiconductor optical metrology:
- **Fresnel OCD Study** — Thin-film reflectance simulation for optical critical dimension metrology
- **Optical Models** — Refractive index models (Cauchy, Sellmeier, Drude-Lorentz)
- **SPC Control Charts** — Statistical process control for semiconductor manufacturing

---

## Tech Stack

![Python](https://img.shields.io/badge/Python-3.10-3776AB?logo=python&logoColor=white)
![NumPy](https://img.shields.io/badge/NumPy-1.24-013243?logo=numpy&logoColor=white)
![SciPy](https://img.shields.io/badge/SciPy-1.11-8CAAE6?logo=scipy&logoColor=white)
![scikit-learn](https://img.shields.io/badge/scikit--learn-1.3-F7931E?logo=scikit-learn&logoColor=white)
![Matplotlib](https://img.shields.io/badge/Matplotlib-3.7-11557C)
![pandas](https://img.shields.io/badge/pandas-2.0-150458?logo=pandas&logoColor=white)
![Jupyter](https://img.shields.io/badge/Jupyter-Notebook-F37626?logo=jupyter&logoColor=white)

---

## Project Structure

```
semiconductor-ai-portfolio/
├── 01_MoSe2_PL_analysis/
│   └── MoSe2_PL_analysis.ipynb        # PL peak detection & FWHM
├── 02_MoSe2_NSOM_defect_mapping/
│   ├── MoSe2_NSOM_analysis.ipynb       # Hyperspectral NSOM analysis
│   └── MoSe2_NSOM_30x30_kmeans.ipynb   # K-Means defect clustering
├── 03_TMD_NSOM_comparison/
│   ├── WSe2_bilayer_CVD_analysis.ipynb  # WSe2 bilayer analysis
│   └── WSe2_CVD_vs_WS2_bilayer_analysis.ipynb
├── metrology/
│   ├── fresnel/                        # OCD thin-film simulation
│   ├── optical_models/                 # Refractive index models
│   └── spc/                            # SPC control charts
├── demo.ipynb                          # Quick-start demo notebook
├── requirements.txt
└── README.md
```

---

## Getting Started

```bash
# Clone the repository
git clone https://github.com/kyb8801/semiconductor-ai-portfolio.git
cd semiconductor-ai-portfolio

# Install dependencies
pip install -r requirements.txt

# Launch Jupyter and explore
jupyter notebook demo.ipynb
```

---

## Demo

See [`demo.ipynb`](demo.ipynb) for a guided walkthrough of all projects, including:
- PL spectrum peak detection with interactive visualization
- NSOM defect clustering results
- Thin-film reflectance simulation

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

## Contact

- **GitHub:** [kyb8801](https://github.com/kyb8801)
- **Kaggle:** [aioptic](https://www.kaggle.com/aioptic)

---

*Part of my metrology × AI × UQ portfolio.*
