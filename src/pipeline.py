"""
Target-Focused Virtual Screening on Tox21 Library
Target: Bacterial Peptide Deformylase (PDF)
Template Ligand: Actinonin
"""

import os
import deepchem as dc
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

def run_virtual_screen():
    # 1. Define Target Anchor (Actinonin)
    actinonin_smiles = "CCCCCC(CC(=O)NO)C(=O)N1CCCC1C(=O)NC(C(C)C)CO"
    actinonin_mol = Chem.MolFromSmiles(actinonin_smiles)
    if not actinonin_mol:
        raise ValueError("Invalid SMILES format for reference ligand.")
    
    actinonin_fp = AllChem.GetMorganFingerprintAsBitVect(actinonin_mol, radius=2, nBits=2048)
    print("[INFO] Featurized reference ligand (Actinonin) successfully.")

    # 2. Ingest Tox21 Library
    print("[INFO] Fetching and streaming Tox21 dataset from MoleculeNet...")
    tasks, datasets, transformers = dc.molnet.load_tox21(featurizer='Raw')
    train_dataset, _, _ = datasets

    # Map toxicological endpoints
    mmp_idx = tasks.index('SR-MMP')  # Mitochondrial Membrane Potential Disruption
    p53_idx = tasks.index('SR-p53')  # Genotoxicity Indicator

    safe_candidates = []
    skipped_count = 0

    print("[INFO] Initiating ligand similarity metrics and toxicity profiling...")
    
    # 3. Stream and evaluate candidates
    for X, y, w, ids in train_dataset.itersamples():
        if X is None:
            continue
            
        # Calculate chemical features
        fp = AllChem.GetMorganFingerprintAsBitVect(X, radius=2, nBits=2048)
        similarity = DataStructs.TanimotoSimilarity(actinonin_fp, fp)
        
        # Extract true wet-lab assay markers
        is_mito_toxic = (y[mmp_idx] == 1) and (w[mmp_idx] > 0)
        is_dna_toxic = (y[p53_idx] == 1) and (w[p53_idx] > 0)
        
        # Domain Filter: Discard structural matches causing host toxicity
        if not is_mito_toxic and not is_dna_toxic:
            safe_candidates.append({
                'id': ids,
                'similarity': similarity
            })
        else:
            if similarity > 0.35:
                skipped_count += 1

    print(f"[INFO] Excluded {skipped_count} structural analogs due to hostile toxicological profiles.")

    # 4. Sort and return results
    safe_candidates = sorted(safe_candidates, key=lambda x: x['similarity'], reverse=True)
    
    print("\n================ TOP 5 SAFE PDF INHIBITOR CANDIDATES ================")
    for idx, hit in enumerate(safe_candidates[:5], 1):
        print(f"{idx}. ID: {hit['id']:<15} | Tanimoto Score: {hit['similarity']:.4f} | Status: PASS")
    print("=====================================================================\n")

if __name__ == "__main__":
    run_virtual_screen()
