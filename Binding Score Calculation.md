### Binding Score Calculation

The CLI script provides a command-line interface for calculating the binding score between two proteins based on their molecular conformations. The process involves several steps:

1. **Download PDB Files:**
   
   - The script starts by downloading the PDB (Protein Data Bank) files for the specified proteins (`Entry ID 1` and `Entry ID 2`) using the RCSB PDB web service.
   - The PDB files are saved locally in the `static/pdb_files/` directory.

2. **Molecular Conformation Superimposition:**
   
   - The script uses the MDTraj library to load the molecular trajectories (`md.Trajectory`) from the downloaded PDB files for both proteins.
   - It identifies the common atoms in both proteins' structures using the `select` method with the 'protein' keyword.
   - If no common atoms are found, the script prints a message and exits, as the interaction simulation cannot proceed without common atoms.

3. **Calculate Center of Mass (COM) and Translation:**
   
   - The script calculates the mean structure and center of mass (COM) for each protein based on the identified common atoms.
   - It computes the translation vector needed to superimpose the two proteins, aligning their common atoms.

4. **Apply Superimposition:**
   
   - The MDTraj library is used to apply the calculated translation vector to the second protein's trajectory, aligning it with the first protein.

5. **Binding Score Calculation:**
   
   - The binding score is calculated using the Universal Force Field (UFF) implemented in the RDKit library.
   - The script converts the proteins' molecular structures to RDKit Molecule objects and selects a subset of atoms for docking score calculation based on the substructure matching with a simple organic group (e.g., methane, 'C').
   - The selected substructures are then embedded in 3D space and optimized using the UFF method.
   - For each conformation, the script calculates the energy difference between the two proteins using the UFF force fields.
   - The calculated energy differences for each conformation are stored in a list (`docking_scores`).

6. **Result Output:**
   
   - The script prints the number of conformers and the number of valid docking scores.
   - If any energy calculations fail for a specific conformation, a warning message is printed, and the failed conformations are tracked.
   - The final binding score is calculated as the mean of valid docking scores.

7. **Interaction Simulation Visualization:**
   
   - After the binding score calculation, the script attempts to visualize the interaction by superimposing the two proteins' structures.
   - It uses a 3D scatter plot to display the protein structures in different colors, and the plot is shown using Matplotlib.
   - The visualization is provided to help understand the spatial arrangement of the proteins during the interaction.
