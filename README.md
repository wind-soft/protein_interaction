# Protein Docking and Interaction Simulation

This repository contains Python scripts for performing protein docking and interaction simulation. The project consists of two main components: a command-line interface (CLI) script for docking and a web application for interactive use.

## CLI Script (cli.py)

### Dependencies

- Python 3.7 or later
- Required Python packages (install using `pip install -r requirements.txt`):
  - requests
  - numpy
  - mdtraj
  - matplotlib
  - rdkit
  - memory-profiler

### Usage

To use the CLI script, run the following command in the terminal:

```bash
python cli.py <Entry_ID_1> <Entry_ID_2>
```

Replace `<Entry_ID_1>` and `<Entry_ID_2>` with the RCSB Entry ID of the proteins you want to analyze.

#### Example

```bash
python cli.py 6M0J proteinB
```

### Functionality

- Downloads PDB files for the specified proteins.
- Superimposes the protein structures for interaction simulation.
- Calculates docking scores based on molecular conformations.
- Visualizes the interaction results in a 3D plot.

## Web Application (web.py)

### Dependencies

- Python 3.7 or later
- Required Python packages (install using `pip install -r requirements.txt`):
  - Flask
  - requests
  - numpy
  - mdtraj
  - matplotlib
  - rdkit
  - memory-profiler

### Usage

1. Run the web application:

```bash
`python web.py`
```

2. Open a web browser and navigate to `http://localhost:5000`.

3. Enter the names of the proteins you want to analyze in the provided form.

4. Submit the form to see docking scores and interaction simulation results.

### Functionality

- Downloads PDB files for the specified proteins.
- Superimposes the protein structures for interaction simulation.
- Calculates docking scores based on molecular conformations.
- Displays results on a web interface, including docking scores and interaction simulation plots.

## Notes

- The web application saves interaction plots to the `static` folder.
- For detailed profiling information, you can use the `memory_profiler` module by decorating functions with `@profile`.

## Acknowledgments

- The project uses various Python libraries, including mdtraj, matplotlib, rdkit, and Flask.

Feel free to contribute, report issues, or suggest improvements!
