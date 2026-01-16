#!/bin/bash -l
# =============================================================================
# SLURM Job Submission: Topology Generation
# =============================================================================
#
# Generates polymer topology from SMILES string and automatically submits
# the equilibration job.
#
# USAGE:
#   sbatch submit_topology.sh <N_MONOMERS> <Z_LENGTH>
#
# EXAMPLE:
#   sbatch submit_topology.sh 1000 120
#
# =============================================================================

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --job-name="topgen"
#SBATCH --partition=idle
#SBATCH --time=7-0
#SBATCH --mail-user='YOUR_EMAIL@example.com'
#SBATCH --mail-type=ALL
#SBATCH --export=NONE
#SBATCH -o stdout_topgen_%j.txt
#SBATCH -e stderr_topgen_%j.txt
#SBATCH --requeue
#SBATCH --open-mode=append

# -----------------------------------------------------------------------------
# USER CONFIGURATION
# -----------------------------------------------------------------------------

# SMILES string for the monomer
# Options:
#   Ethylene glycol: "OCCOCCO"
#   Hexamethylene:   "OCCCCCCO"
#   Fluorinated:     "OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)CO"
SMILES="OCCOCCO"
PREFIX="Ediol"

# -----------------------------------------------------------------------------
# Environment Setup for DARWIN at U. Delaware (modify for your cluster)
# -----------------------------------------------------------------------------

vpkg_require anaconda/5.3.1:python2 openmpi/4.0.5
conda activate mosdef_env


. /opt/shared/slurm/templates/libexec/openmpi.sh

# -----------------------------------------------------------------------------
# Run Topology Generation
# -----------------------------------------------------------------------------

N_MONOMERS=$1
Z_LENGTH=$2

python3 generate_topology.py "${SMILES}" ${N_MONOMERS} "${PREFIX}_${N_MONOMERS}monomer.data"

# -----------------------------------------------------------------------------
# Submit Next Stage (Equilibration)
# -----------------------------------------------------------------------------

sbatch submit_equilibration.sh lammps_equilibration.in "${PREFIX}_${N_MONOMERS}monomer.data" "${Z_LENGTH}"
