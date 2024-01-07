from flask import Flask, render_template, request
import os
import requests
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
from memory_profiler import profile

import matplotlib
matplotlib.use('Agg')

app = Flask(__name__)

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

    # Save the plot to a temporary file
    tmp_file_path = 'static/interaction_plot.png'
    plt.savefig(tmp_file_path)
    plt.close()

    return tmp_file_path  # Return the file path to be embedded in the HTML

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

    docking_scores = []  # Fix: declare docking_scores list here

    calculation_steps_data = []

    for confId in range(mol1_subset.GetNumConformers()):
        ff1 = AllChem.UFFGetMoleculeForceField(mol1_subset, confId=confId)
        ff2 = AllChem.UFFGetMoleculeForceField(mol2_subset, confId=confId)

        step_data = {
            'step': f'Conformation {confId + 1} Energy Calculation',
            'protein_pair': protein_pair,
            'energy_diff': None
        }

        if ff1 is not None and ff2 is not None:
            energy_diff = ff1.CalcEnergy() - ff2.CalcEnergy()
            docking_scores.append(energy_diff)
            step_data['energy_diff'] = energy_diff
        else:
            print(f"Warning: Energy calculation failed for conformation {confId}")
            if ff1 is None:
                print(f"FF1 is None for conformation {confId}")
            if ff2 is None:
                print(f"FF2 is None for conformation {confId}")

        calculation_steps_data.append(step_data)

    mol1_subset = None
    mol2_subset = None

    valid_scores = [score for score in docking_scores if score is not None]

    final_score = np.mean(valid_scores) if valid_scores else None

    return protein_pair, final_score, calculation_steps_data

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        protein_name1 = request.form['protein_name1']
        protein_name2 = request.form['protein_name2']

        protein_names_batch = [(protein_name1, protein_name2)]

        results = []

        for protein_pair in protein_names_batch:
            result = docking_score(protein_pair)
            results.append(result)

        for protein_pair, score, calculation_steps_data in results:
            docking_result = f'Docking Score for {protein_pair}: {score}'

        protein1_pdb = download_pdb(protein_name1)
        protein2_pdb = download_pdb(protein_name2)

        protein1_coords, protein2_coords = superimpose_proteins(protein1_pdb, protein2_pdb)

        if not np.any(protein1_coords) or not np.any(protein2_coords):
            interaction_result = "No common atoms found, cannot perform interaction simulation."
            interaction_plot_path = None
        else:
            interaction_result = "Interaction simulation successful."
            interaction_plot_path = visualize_interaction_result(protein1_coords, protein2_coords)

        return render_template('index.html', docking_result=docking_result, interaction_result=interaction_result, interaction_plot_path=interaction_plot_path, calculation_steps_data=calculation_steps_data)

    return render_template('index.html', docking_result=None, interaction_result=None, interaction_plot_path=None, calculation_steps_data=None)

if __name__ == '__main__':
    app.run(debug=False)
