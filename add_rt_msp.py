import pandas as pd
import argparse
import os
from rdkit import Chem
import numpy as np

from typing import Optional, Tuple
import math
import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi
from tqdm import tqdm
def add_predicted_rt_from_df(input_file: str,
                             output_file: str,
                             df: pd.DataFrame,
                             rewrite_smiles: bool = False,
                             decimals: int = 4) -> Tuple[int, int]:
    """
    Insert 'Predicted RT: <value>' before each 'Num Peaks' line in an MSP file.
    RT is looked up by InChIKey computed from the SMILES line.
    
    Returns:
        (total_spectra, total_predicted_added)
    """
    # Build lookup: InChIKey (upper, stripped) -> predicted_rt
    df_lk = df.copy()
    df_lk['inchikey'] = df_lk['inchikey'].astype(str).str.strip().str.upper()
    rt_map = dict(zip(df_lk['inchikey'], df_lk['predicted_rt']))

    def _format_rt(val: Optional[float]) -> Optional[str]:
        if val is None:
            return None
        try:
            fv = float(val)
            if math.isnan(fv):
                return None
            return f"{fv:.{decimals}f}"
        except Exception:
            return None

    # per-entry state
    current_inchikey = None
    smiles_valid = False
    pending_rt_val: Optional[float] = None

    # counters
    total_spectra = 0
    total_predicted_added = 0

    with open(input_file, "r", encoding="utf-8", errors="ignore") as fin, \
         open(output_file, "w", encoding="utf-8") as fout:

        for raw in tqdm(fin, desc="Processing MSP"):
            line = raw
            s = line.strip()

            # Blank line -> end of entry, reset state
            if s == "":
                current_inchikey = None
                smiles_valid = False
                pending_rt_val = None
                fout.write(line)
                continue

            # Parse SMILES line
            if s.lower().startswith("smiles"):
                parts = line.split(":", 1)
                smiles_str = parts[1].strip() if len(parts) == 2 else s.split(None, 1)[-1].strip()
                mol = Chem.MolFromSmiles(smiles_str)
                if mol is None:
                    # invalid SMILES
                    smiles_valid = False
                    current_inchikey = None
                    pending_rt_val = None
                    fout.write(line)  # keep original SMILES text
                else:
                    smiles_valid = True
                    can_smi = Chem.MolToSmiles(mol, canonical=True)
                    try:
                        ikey = inchi.MolToInchiKey(mol)
                        current_inchikey = ikey.strip().upper()
                    except Exception:
                        current_inchikey = None
                    pending_rt_val = rt_map.get(current_inchikey) if current_inchikey else None

                    if rewrite_smiles:
                        fout.write(f"SMILES: {can_smi}\n")
                    else:
                        fout.write(line)
                continue

            # Right before spectrum starts
            if s.lower().startswith("num peaks"):
                total_spectra += 1

                if not smiles_valid:
                    fout.write("Predicted RT: -1\n")
                else:
                    rt_formatted = _format_rt(pending_rt_val)
                    if rt_formatted is None:
                        fout.write("Predicted RT: -1\n")
                    else:
                        fout.write(f"Predicted RT: {rt_formatted}\n")
                        total_predicted_added += 1

                fout.write(line)
                continue

            # passthrough
            fout.write(line)

    print(f"Done. Spectra processed: {total_spectra}, Predicted RT added: {total_predicted_added}")
    return total_spectra, total_predicted_added
def add_predicted_rt_msp(data_name, col, msp_dir, prediction_dir):
    all_predicted_rt = pd.read_csv(os.path.join(prediction_dir, f'all_predicted_rt_{col}.csv'))
    add_predicted_rt_from_df(os.path.join(msp_dir, f'{data_name}-raw.msp'), os.path.join(msp_dir, f'{data_name}-{col}-predicted.msp'), all_predicted_rt)
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--library_name', type=str)
    parser.add_argument('--column_name', type=str)
    parser.add_argument('--msp_dir', type=str, default = '/quobyte/metabolomicsgrp/fanzhou/raw_msp')
    parser.add_argument('--prediction_dir', type=str, default = '/quobyte/metabolomicsgrp/fanzhou/predictions')
    args = parser.parse_args()
    library_name = args.library_name
    column_name = args.column_name
    msp_dir = args.msp_dir
    prediction_dir = args.prediction_dir
    add_predicted_rt_msp(library_name, column_name, msp_dir, prediction_dir)
    # will read the {library_name}-raw.msp from msp_dir, look up for rt in all_predicted_rt_{column_name}.csv in prediction dir, then add the rt to the msp file
    # then save the msp file to {library_name}-{column_name}-predicted.msp in msp_dir

    print('done')
    # from utils.chemutils import calc_descriptors_df
