# Target-Focused Virtual Screening Pipeline for Bacterial Peptide Deformylase (PDF) Inhibitors

This repository features a computational biochemistry pipeline designed to discover hit candidates for **Peptide Deformylase (PDF)**—a critical bacterial enzyme target for novel antibiotics—by screening the **Tox21 chemical library** (10,000+ compounds).

The core framework leverages ligand-based **Molecular Similarity Search** combined with multi-task experimental filtering to isolate structurally viable, non-cytotoxic hits.

## 🧬 Biological & Computational Design

### 1. Template Target Anchoring
Standard toxicity repositories lack specific biological labels for bacterial proteins like PDF. To circumvent this, this pipeline utilizes **Actinonin**—a known, potent natural peptide deformylase inhibitor—as a template structural query. 

### 2. Molecular Featurization
Chemical profiles (SMILES representations) are mapped into machine-readable mathematical formats by generating **2048-bit Morgan Fingerprints** (radius=2) via RDKit. Structural overlaps are computed using **Tanimoto Similarity Coefficients**.

### 3. Toxicological Host Filtering
To ensure identified hits can function as safe therapeutic options, the pipeline integrates real wet-lab experimental endpoints from the Tox21 panel. Compounds exhibiting high structural similarity are automatically **flagged and dropped** if they express:
*   `SR-MMP`: Mitochondrial Membrane Potential Disruption (Cellular metabolic failure).
*   `SR-p53`: Genotoxicity and DNA Damage responses.

---

## 🛠️ Technology Stack
*   **Language:** Python 3.10+
*   **Cheminformatics Framework:** RDKit
*   **AI/Data Engineering:** DeepChem, NumPy
*   **Machine Learning Ecosystem:** Scikit-Learn

---

## 🚀 Getting Started

### Installation
Clone this repository and install the dependencies:
```bash
git clone https://github.com
cd tox21-pdf-virtual-screening
pip install -r requirements.txt
```

### Run the Screening Pipeline
Execute the virtual screening engine directly:
```bash
python src/pipeline.py
```

## 📈 Architecture Overview
1. **Ingest Template:** Parse and featurize the reference PDF inhibitor ligand (Actinonin).
2. **Stream Library:** Query DeepChem's MoleculeNet client to isolate structural geometries from the Tox21 repository.
3. **Compute Vector Distances:** Quantify Tanimoto similarity metrics across all fingerprints.
4. **Enforce Biological Logic:** Discard any structural analogs displaying mammalian cytotoxicity matrices.
5. **Rank Hits:** Output top target candidates prioritized by pure ligand-binding safety profiles.
