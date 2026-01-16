#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Topology Generator for Polyurethane Simulations

Generates LAMMPS-compatible topology files from SMILES strings using
the MoSDeF (Molecular Simulation Design Framework) toolkit.

The workflow:
    1. Parse SMILES string to create a molecular structure
    2. Energy minimize the monomer geometry using OPLS-AA
    3. Fill a simulation box with multiple monomers
    4. Assign OPLS-AA force field parameters
    5. Write LAMMPS data file

USAGE:
    python generate_topology.py "SMILES_STRING" N_MONOMERS OUTPUT_FILE

EXAMPLE:
    python generate_topology.py "OCCOCCO" 1000 "Ediol_1000monomer.data"
    
    Common SMILES strings for diols:
    - Ethylene glycol diol: "OCCOCCO"
    - Hexamethylene diol:   "OCCCCCCO"  
    - Fluorinated diol:     "OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)CO"

REQUIREMENTS:
    conda env create -f mosdef_env.yml
    conda activate mosdef_env

Author: Tristan Myers
Created: April 2025
"""

import sys
import mbuild as mb
from foyer import Forcefield
import gmso.external
import gmso.formats


def generate_topology(smiles_string, n_monomers, output_filename, density=20):
    """
    Generate a LAMMPS data file from a SMILES string.
    
    Parameters
    ----------
    smiles_string : str
        SMILES representation of the monomer (e.g., "OCCOCCO")
    n_monomers : int
        Number of monomers to place in the simulation box
    output_filename : str
        Output LAMMPS data file name
    density : float, optional
        Initial packing density in kg/m³ (default: 20, deliberately low 
        for initial structure to avoid overlaps)
    """
    print(f"Creating monomer from SMILES: {smiles_string}")
    monomer = mb.load(smiles_string, smiles=True)
    
    print("Energy minimizing monomer structure...")
    monomer.energy_minimize(forcefield='oplsaa', steps=10**4)
    
    print(f"Filling box with {n_monomers} monomers...")
    full_box = mb.fill_box(compound=monomer, n_compounds=n_monomers, density=density)
    
    print("Applying OPLS-AA force field parameters...")
    opls = Forcefield(name='oplsaa')
    full_box = opls.apply(full_box)
    
    print(f"Writing LAMMPS data file: {output_filename}")
    top_box = gmso.external.from_parmed(full_box)
    gmso.formats.write_lammpsdata(top_box, output_filename)
    
    print("Done!")
    return full_box


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(__doc__)
        print("\nERROR: Expected 3 arguments:")
        print("  python generate_topology.py SMILES N_MONOMERS OUTPUT_FILE")
        sys.exit(1)
    
    input_smiles = sys.argv[1]
    n_monomers = int(sys.argv[2])
    output_filename = sys.argv[3]
    
    generate_topology(input_smiles, n_monomers, output_filename)
