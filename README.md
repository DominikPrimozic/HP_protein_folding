# HP_protein_folding

This repository contains a simple implementation of the HP model for protein folding in 2D and 3D. The model is based on the lattice representation of proteins and uses Monte Carlo methods to simulate the folding process.

## The Model

The HP (Hydrophobic-Polar) model is a lattice-based approach for simulating protein folding. In this model:
- **Each amino acid residue** is represented as a point (2D) or sphere (3D) on a rectangular lattice.
- **Residue Classification**: 
  - Hydrophobic residues are labeled as H.
  - Polar residues are labeled as P.
- **Folding Energy**:
  - Favorable H-H spatial neighbors are rewarded with a negative energy contribution.
  - Hydrophobic collapse drives the protein's folding process.

### Simulation
The folding process is simulated using:
1. **Monte Carlo Sampling**: Randomly explores conformational space.
2. **Energy Evaluation**:
   - Any structure that improves energy is immediately accepted.
   - Unfavorable conformations are accepted probabilistically using the Metropolis algorithm:
     - The probability of accepting an unfavorable move decreases with worsening energy and lower temperatures.
     - This ensures the simulation can escape local minima.

## Implementation
The core implementation of the HP model is written in **C++**, leveraging its computational efficiency for the simulation. Additionally:
- A simple file reader is included to load input sequences of amino acids and lattice parameters.
- Input files are preconfigured for testing the simulation.

## Visualization
Since C++ is not well-suited for creating visuals, visualization scripts are written in **Python** using the **Plotly** library. These scripts:
- Visualize the conformations during simulation in both 2D and 3D.
- Provide an interactive interface to analyze the folding pathways.

## How to Use
1. **Compile the C++ Code**:
   - Use a compiler like `g++` to compile the main simulation program.

2. **Run the Python Visualization**:
   - Use the output from the C++ simulation as input for the Python script.
   - Example:
     ```bash
     python 3D_animator.py 
     ```
