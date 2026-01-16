#!/bin/bash -l
# =============================================================================
# SLURM Job Submission: System Equilibration
# =============================================================================
#
# Runs NPT equilibration of polymer system using LAMMPS.
# Supports automatic restart from checkpoint files.
#
# USAGE:
#   sbatch submit_equilibration.sh <LAMMPS_INPUT> <DATA_FILE> <Z_LENGTH>
#
# EXAMPLE:
#   sbatch submit_equilibration.sh lammps_equilibration.in system.data 120
#
# =============================================================================


#SBATCH --ntasks=30
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --tmp=1M
#SBATCH --job-name="equilibration"
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

# LAMMPS executable name (modify for your installation)
LAMMPS_EXE="lmp"

# For custom builds:
# LAMMPS_EXE="lmp_mpi"

# -----------------------------------------------------------------------------
# Environment Setup (modify for your cluster)
# -----------------------------------------------------------------------------

# Example for DARWIN cluster:
vpkg_devrequire fftw/3.3.9:gcc-10.1.0,openmpi-4.0.5
vpkg_require anaconda/5.3.1:python3 openmpi/4.0.5

# For other clusters:
# module load lammps
# module load openmpi

. /opt/shared/slurm/templates/libexec/openmpi.sh

# -----------------------------------------------------------------------------
# Setup Working Directory
# -----------------------------------------------------------------------------

file_loc=$PWD
mkdir -p "$file_loc/sim_output"
mkdir -p "$file_loc/restart_files"

# -----------------------------------------------------------------------------
# Determine Start Mode (fresh or restart)
# -----------------------------------------------------------------------------

restart_file=$(ls -t1 "$file_loc/restart_files/" 2>/dev/null | head -n 1)
if [ -n "$restart_file" ]; then
    datafile="restart_files/$restart_file"
    restart_bool=1
    echo "Restarting from: $datafile"
else
    datafile="$2"
    restart_bool=0
    echo "Starting fresh from: $datafile"
fi

# Track restart count for log naming
if [ -z "${SLURM_RESTART_COUNT}" ]; then
    restart_count=0
else
    restart_count=${SLURM_RESTART_COUNT}
fi

# -----------------------------------------------------------------------------
# Run LAMMPS
# -----------------------------------------------------------------------------

LAMMPS_INPUT="$1"
Z_LENGTH="$3"

echo "Using MPI launcher: ${UD_MPIRUN:-mpirun}"
echo "LAMMPS input: $LAMMPS_INPUT"
echo "Z-axis length: $Z_LENGTH Angstroms"

${UD_MPIRUN:-mpirun} ${LAMMPS_EXE} \
    -in "$LAMMPS_INPUT" \
    -var datafile "$file_loc/$datafile" \
    -var vel_seed $RANDOM \
    -var restart_bool "$restart_bool" \
    -var zlen "$Z_LENGTH" \
    -log "log_${restart_count}.lammps"

exit $?
