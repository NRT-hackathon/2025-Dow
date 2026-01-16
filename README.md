# Polyurethane Thermal Conductivity Simulation Workflow

A complete workflow for simulating dense polyurethane and calculating thermal conductivity using LAMMPS molecular dynamics. The code, including this README, is written for a user with working knowledge of Python, LAMMPS, and the SLURM job management system.

## NOTE ##
This workflow is only currently set up to simulate & analyze collections of individual polyol monomers, rather than polyurethane chains.

#####

## Overview

This workflow enables:
1. **Topology Generation** — Create polymer structures from SMILES strings using MoSDeF
2. **System Equilibration** — Relax packed monomers to equilibrium density via NPT
3. **Thermal Property Calculation** — Measure thermal conductivity (Müller-Plathe NEMD) and heat capacity (DOS/VACF)
4. **Post-Processing** — Analyze simulation output with Python

## Quick Start

```bash
# 1. Set up environment
conda env create -f mosdef_env.yml
conda activate mosdef_env

# 2. Generate topology (1000 ethylene glycol monomers)
python generate_topology.py "OCCOCCO" 1000 "system.data"

# 3. Run equilibration (on HPC cluster)
sbatch submit_equilibration.sh lammps_equilibration.in system.data 120

# 4. Run thermal property measurement
sbatch submit_thermal_props.sh lammps_thermal_props.in stage3.data

# 5. Analyze results
python run_analysis.py
```

## File Descriptions

### Python Scripts
| File | Description |
|------|-------------|
| `generate_topology.py` | Creates LAMMPS data files from SMILES strings |
| `analysis_functions.py` | Reusable library for post-processing |
| `run_analysis.py` | Example analysis script with usage patterns |

### LAMMPS Input Scripts
| File | Description |
|------|-------------|
| `lammps_equilibration.in` | NVT minimization → box deformation → NPT relaxation |
| `lammps_thermal_props.in` | VACF sampling + Müller-Plathe NEMD |

### SLURM Submission Scripts
| File | Description |
|------|-------------|
| `submit_topology.sh` | Job submission for topology generation |
| `submit_equilibration.sh` | Job submission for equilibration |
| `submit_thermal_props.sh` | Job submission for thermal measurements |

## Requirements

### Python Environment

Create the conda environment:
```bash
conda env create -f mosdef_env.yml
conda activate mosdef_env
```

Required packages (included in `mosdef_env.yml`):
- Python 3.9
- mBuild
- Foyer
- GMSO
- OpenBabel

For analysis, also install:
```bash
conda install numpy scipy matplotlib
pip install MDAnalysis  # optional, for trajectory analysis
```

### LAMMPS

LAMMPS must be compiled with these packages:
- `MOLECULE` — molecular systems
- `KSPACE` — long-range electrostatics (PPPM/Ewald)
- `EXTRA-MOLECULE` - additional molecular features
- `MISC` - additional features

**Installation options:**


1. **Build from source (recommended for production):**
   ```bash
   git clone https://github.com/lammps/lammps.git
   cd lammps/src
   make yes-molecule yes-kspace yes-extra-molecule yes-misc
   make mpi
   ```

2. **Use cluster module, if suitable:**
   ```bash
   module load lammps
   ```

## Workflow Details

### Step 1: Topology Generation

The `generate_topology.py` script creates LAMMPS-compatible data files:

```bash
python generate_topology.py "SMILES" N_MONOMERS "output.data"
```

**Common SMILES strings:**
| Material | SMILES |
|----------|--------|
| Ethylene glycol diol | `OCCOCCO` |
| Hexamethylene diol | `OCCCCCCO` |
| Fluorinated diol | `OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)CO` |

### Step 2: System Equilibration

The equilibration runs ~34 ns with these stages:
1. **NVT minimization** (0–1 ns) — Energy minimize and thermalize
2. **Box deformation** (1–3 ns) — Compress to target density
3. **NVT equilibration** (3–4 ns) — Relax at target density
4. **NPT equilibration** (4–34 ns) — Full equilibration at 350 K, 1 atm

### Step 3: Thermal Conductivity (NEMD)

Uses the Müller-Plathe method:
- Kinetic energy is exchanged between "hot" and "cold" slabs
- Temperature gradient develops along z-axis
- Thermal conductivity: κ = Q/(A × dT/dz)

**Stages:**
1. **VACF sampling** (optional) — Collect velocity autocorrelation for DOS
2. **NEMD relaxation** (0–2 ns) — Reach steady-state temperature gradient
3. **NEMD production** (2–3 ns) — Sample for averaging

### Step 4: Post-Processing

Edit `run_analysis.py` to set your file paths, then:

```python
from analysis_functions import *

# Load data
profiles, timesteps = readdata('profile.mp')
thermo = readfile_log('log.lammps', np.arange(0, 6))

# Calculate thermal conductivity
kappa = get_thermcond_NEMD(
    thermo[-1000:, 5],
    (profiles[-1000:], timesteps[-1000:]),
    delta_t=1000,
    z_axis_len=120,
    Xsectional_len=45
)

print(f"κ = {np.mean(kappa):.3f} ± {np.std(kappa):.3f} W/(m·K)")
```

## Output Files

| File | Contents |
|------|----------|
| `sim_output/traju.txt` | Unwrapped atomic trajectory |
| `sim_output/profile.mp` | Temperature profiles (NEMD) |
| `vacf.txt` | Velocity autocorrelation function |
| `sys.thermo` | Thermodynamic properties vs time |
| `restart_files/` | Checkpoint files for restart |
| `stage2.data`, `stage3.data` | Intermediate configurations |

## Customization

### Changing the Polymer

Edit `submit_topology.sh`:
```bash
SMILES="YOUR_SMILES_HERE"
PREFIX="your_prefix"
```

### Changing Temperature

Edit the LAMMPS input files:
```
variable temp equal YOUR_TEMP
```

### Cluster Configuration

Edit the `submit_*.sh` scripts to match your cluster:
- Replace `vpkg_require` with `module load`
- Update `LAMMPS_EXE` to your executable name
- Adjust `--partition`, `--mail-user`, etc.

## Troubleshooting

**"LAMMPS data file not found"**
- Check the path in your submission script
- Ensure topology generation completed successfully

**"restart file not found" on first run**
- This is normal — simulation will start from the data file

**Thermal conductivity values seem wrong**
- Check that `delta_t` matches your thermo output frequency
- Verify `z_axis_len` matches your simulation box
- Ensure you're using production data (after steady-state is reached)

## Citation

If you use this workflow, please cite:
- MoSDeF: [Cummings et al., Mol. Simul. 2021](https://doi.org/10.1080/08927022.2021.1957877)
- LAMMPS: [Thompson et al., Comp. Phys. Comm. 2022](https://doi.org/10.1016/j.cpc.2021.108171)
- OPLS-AA: [Jorgensen et al., JACS 1996](https://doi.org/10.1021/ja9621760)

## License

[Add your license here]

## Author

Tristan Myers  
University of Delaware  
May 2025
