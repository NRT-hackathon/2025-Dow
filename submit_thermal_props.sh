#!/bin/bash -l
# =============================================================================
# SLURM Job Submission: Thermal Property Measurements
# =============================================================================
#
# Runs NEMD simulation for thermal conductivity using Müller-Plathe method.
#
# USAGE:
#   sbatch submit_thermal_props.sh <LAMMPS_INPUT> <DATA_FILE>
#
# EXAMPLE:
#   sbatch submit_thermal_props.sh lammps_thermal_props.in equilibrated.data
#
# =============================================================================


#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --tmp=1M
#SBATCH --job-name="thermal_props"
#SBATCH --partition=idle
#SBATCH --time=7-00:00:00
#SBATCH --output=stdout_%j.txt
#SBATCH --error=stderr_%j.txt
#SBATCH --mail-user='YOUR_EMAIL@example.com'
#SBATCH --mail-type=ALL
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --requeue

# -----------------------------------------------------------------------------
# USER CONFIGURATION
# -----------------------------------------------------------------------------

LAMMPS_EXE="lmp"

# -----------------------------------------------------------------------------
# Environment Setup
# -----------------------------------------------------------------------------

vpkg_devrequire fftw/3.3.9:gcc-10.1.0,openmpi-4.0.5
vpkg_require anaconda/5.3.1:python3 openmpi/4.0.5

. /opt/shared/slurm/templates/libexec/openmpi.sh

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

file_loc=$PWD
mkdir -p "$file_loc/sim_output"
mkdir -p "$file_loc/restart_files"

# -----------------------------------------------------------------------------
# Run LAMMPS
# -----------------------------------------------------------------------------

LAMMPS_INPUT="$1"
DATA_FILE="$2"

echo "Running thermal property calculation"
echo "Input: $LAMMPS_INPUT"
echo "Data: $DATA_FILE"

${UD_MPIRUN:-mpirun} ${LAMMPS_EXE} \
    -in "$LAMMPS_INPUT" \
    -var datafile "$DATA_FILE" \
    -var vel_seed $RANDOM

exit $?
