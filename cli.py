import sys
import os
import requests
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
from memory_profiler import profile

def download_pdb(protein_name):
    pdb_url = f'https://files.rcsb.org/download/{protein_name.upper()}.pdb'
    pdb_file_path = f'static/pdb_files/{protein_name}.pdb'

    os.makedirs(os.path.dirname(pdb_file_path), exist_ok=True)

    if not os.path.exists(pdb_file_path):
        response = requests.get(pdb_url)
        with open(pdb_file_path, 'wb') as pdb_file:
            pdb_file.write(response.content)

    return pdb_file_path

def superimpose_proteins(protein1_pdb, protein2_pdb):
    traj1 = md.load(protein1_pdb)
    traj2 = md.load(protein2_pdb)

    common_atoms = set(traj1.topology.select('protein')) & set(traj2.topology.select('protein'))

    if not common_atoms:
        print("No common atoms found, cannot perform interaction simulation.")
        return None, None

    mean_structure1 = np.mean(traj1.xyz[:, list(common_atoms), :], axis=1)
    mean_structure2 = np.mean(traj2.xyz[:, list(common_atoms), :], axis=1)

    com_structure1 = np.mean(mean_structure1, axis=0)
    com_structure2 = np.mean(mean_structure2, axis=0)

    translation_vector = com_structure1 - com_structure2

    traj2_aligned = md.Trajectory(traj2.xyz + translation_vector, traj2.top)

    return traj1.xyz, traj2_aligned.xyz

def visualize_interaction_result(protein1_coords, protein2_coords, interaction_strength=None):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    if interaction_strength is not None:
        color_map = plt.cm.get_cmap('viridis')
        colors = color_map(interaction_strength)
        scatter1 = ax.scatter(*protein1_coords.T, c=colors[0], s=50, label='Protein 1')
        scatter2 = ax.scatter(*protein2_coords.T, c=colors[1], s=50, label='Protein 2')
    else:
        scatter1 = ax.scatter(*protein1_coords.T, s=50, label='Protein 1')
        scatter2 = ax.scatter(*protein2_coords.T, s=50, label='Protein 2')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Protein Interaction')

    plt.legend()
    plt.show()

@profile
def docking_score(protein_pair):
    protein_name1, protein_name2 = protein_pair

    pdb_file1 = download_pdb(protein_name1)
    pdb_file2 = download_pdb(protein_name2)

    mol1 = Chem.MolFromPDBFile(pdb_file1)
    mol2 = Chem.MolFromPDBFile(pdb_file2)

    mol1_atoms_subset = mol1.GetSubstructMatch(Chem.MolFromSmiles('C'))
    mol2_atoms_subset = mol2.GetSubstructMatch(Chem.MolFromSmiles('C'))

    mol1_subset = Chem.MolFromSmiles(''.join([mol1.GetAtomWithIdx(atom_idx).GetSmarts() for atom_idx in mol1_atoms_subset]))
    mol2_subset = Chem.MolFromSmiles(''.join([mol2.GetAtomWithIdx(atom_idx).GetSmarts() for atom_idx in mol2_atoms_subset]))

    mol1_subset = Chem.AddHs(mol1_subset)
    mol2_subset = Chem.AddHs(mol2_subset)

    AllChem.EmbedMultipleConfs(mol1_subset, numConfs=10, randomSeed=42)
    AllChem.EmbedMultipleConfs(mol2_subset, numConfs=10, randomSeed=42)

    for confId in range(mol1_subset.GetNumConformers()):
        AllChem.UFFOptimizeMolecule(mol1_subset, confId=confId, maxIters=500)
    for confId in range(mol2_subset.GetNumConformers()):
        AllChem.UFFOptimizeMolecule(mol2_subset, confId=confId, maxIters=500)

    docking_scores = []
    for confId in range(mol1_subset.GetNumConformers()):
        ff1 = AllChem.UFFGetMoleculeForceField(mol1_subset, confId=confId)
        ff2 = AllChem.UFFGetMoleculeForceField(mol2_subset, confId=confId)

        if ff1 is not None and ff2 is not None:
            energy_diff = ff1.CalcEnergy() - ff2.CalcEnergy()
            docking_scores.append(energy_diff)
        else:
            print(f"Warning: Energy calculation failed for conformation {confId}")
            if ff1 is None:
                print(f"FF1 is None for conformation {confId}")
            if ff2 is None:
                print(f"FF2 is None for conformation {confId}")

    print(f"Number of conformers: {mol1_subset.GetNumConformers()}")
    print(f"Number of docking scores: {len(docking_scores)}")

    mol1_subset = None
    mol2_subset = None

    valid_scores = [score for score in docking_scores if score is not None]

    final_score = np.mean(valid_scores) if valid_scores else None

    return protein_pair, final_score

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python 2-in-1.py <protein_name1> <protein_name2>")
        sys.exit(1)

    protein_name1 = sys.argv[1]
    protein_name2 = sys.argv[2]

    # Example: Batch processing for docking scores
    protein_names_batch = [
        (protein_name1, protein_name2),
        # Add more pairs as needed
    ]

    results = []

    for protein_pair in protein_names_batch:
        result = docking_score(protein_pair)
        results.append(result)

    for protein_pair, score in results:
        print(f'Docking Score for {protein_pair}: {score}')

    # Interaction simulation using superimposition
    protein1_pdb = download_pdb(protein_name1)
    protein2_pdb = download_pdb(protein_name2)

    protein1_coords, protein2_coords = superimpose_proteins(protein1_pdb, protein2_pdb)

    if not np.any(protein1_coords) or not np.any(protein2_coords):
        print("No common atoms found, cannot perform interaction simulation.")
    else:
        visualize_interaction_result(protein1_coords, protein2_coords)
