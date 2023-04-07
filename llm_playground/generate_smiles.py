import json
from typing import List

from rdkit import Chem


def generate_smiles_from_mol(mol: Chem.Mol, times: int) -> List[str]:
    smiles_list = []
    for _ in range(times):
        # Mol オブジェクトを SMILES に変換（ランダムな原子を始点とする）
        smiles = Chem.MolToSmiles(mol, doRandom=True)
        smiles_list.append(smiles)
    return smiles_list


# # 使用例
# mol = Chem.MolFromSmiles("C1CC2=C3C(=CC=C2)C(=C1)C=C3")
# times = 5
# smiles_list = generate_smiles_from_mol(mol, times)

# for idx, smiles in enumerate(smiles_list):
#     print(f"SMILES {idx + 1}: {smiles}")

sdf_file = "resources/CBLB_inhibitors_vsF.sdf"  # ディレクトリをマニュアルで追加
sdf_supplier = Chem.SDMolSupplier(sdf_file)

output_data = []

for mol in sdf_supplier:
    if mol is None:
        continue

    compound_name = mol.GetProp("COMPOUND_NAMES")
    reference = mol.GetProp("Reference")
    ic50_range = mol.GetProp("IC50_range_nM")

    smiles_list = generate_smiles_from_mol(mol, times=10)

    mol_data = {
        "COMPOUND_NAMES": compound_name,
        "SMILES": smiles_list,
        "Reference": reference,
        "IC50_range_nM": ic50_range,
    }
    output_data.append(mol_data)

# output_data を JSON 形式で保存
with open("resources/CBLB_inhibitors_vsF.json", "w") as f:  # ディレクトリをマニュアルで追加
    json.dump(output_data, f, indent=4)  # indent=4 をマニュアルで追加
