#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis Script for Polyurethane Thermal Properties

This script demonstrates how to use the analysis_functions module to
calculate thermal conductivity and related properties from LAMMPS
simulation output.

USAGE:
    1. Update the file paths in the "USER CONFIGURATION" section
    2. Uncomment the analyses you want to run
    3. Run: python run_analysis.py

Author: Tristan Myers
"""

import numpy as np
import matplotlib.pyplot as plt
from analysis_functions import (
    readdata,
    readfile_log,
    read_vac_data,
    compute_fft,
    calc_quantum_correction,
    get_thermcond_NEMD,
    EMA_series
)

# =============================================================================
# USER CONFIGURATION
# =============================================================================

# Directory containing simulation output
SIM_DIR = "/path/to/your/simulation/"

# Simulation parameters
TEMPERATURE = 350       # K
Z_AXIS_LENGTH = 120     # Angstroms (simulation box z-length)
CROSS_SECTION = 45      # Angstroms (x and y dimensions)
SAMPLING_INTERVAL = 1000  # fs (timestep * thermo frequency)

# File paths (update these to match your simulation output)
PROFILE_FILE = SIM_DIR + "profile.mp"        # Temperature profile from ave/chunk
LOG_FILE = SIM_DIR + "log.lammps"            # LAMMPS log file
VACF_FILE = SIM_DIR + "vacf.txt"             # Velocity autocorrelation output


# =============================================================================
# EXAMPLE 1: Calculate Thermal Conductivity from NEMD Simulation
# =============================================================================

def example_thermal_conductivity():
    """
    Calculate thermal conductivity from Müller-Plathe NEMD simulation.
    
    Required files:
        - Temperature profiles from 'fix ave/chunk' (profile.mp)
        - LAMMPS log file with energy transfer from 'fix thermal/conductivity'
    """
    print("=" * 60)
    print("THERMAL CONDUCTIVITY CALCULATION")
    print("=" * 60)
    
    # Load temperature profile data
    # Columns: [bin_coord, temperature]
    profiles, timesteps = readdata(PROFILE_FILE)
    print(f"Loaded {len(timesteps)} temperature profiles")
    
    # Load thermodynamic data from log file
    # Column indices: [step, temp, epair, etotal, f_3 (energy), v_tdiff, ...]
    # Adjust indices based on your thermo_style output
    thermo = readfile_log(LOG_FILE, np.arange(0, 6))
    print(f"Loaded {len(thermo)} thermo entries")
    
    # Use last 1000 samples for production average
    n_production = 1000
    
    # Calculate thermal conductivity
    kappa = get_thermcond_NEMD(
        deltaE_MP=thermo[-n_production:, 5],  # Energy transfer column
        MP_profile=(profiles[-n_production:], timesteps[-n_production:]),
        delta_t=SAMPLING_INTERVAL,
        z_axis_len=Z_AXIS_LENGTH,
        Xsectional_len=CROSS_SECTION,
        avg_bool=False
    )
    
    # Apply exponential smoothing for cleaner visualization
    kappa_smoothed = EMA_series(kappa, alpha=0.05)
    
    # Calculate statistics
    kappa_avg = np.mean(kappa)
    kappa_std = np.std(kappa)
    
    print(f"\nThermal Conductivity: {kappa_avg:.4f} ± {kappa_std:.4f} W/(m·K)")
    
    # Plot results
    plt.figure(figsize=(10, 6))
    plt.plot(kappa, 'lightgray', alpha=0.5, label='Raw')
    plt.plot(kappa_smoothed, 'k', linewidth=2, label='Smoothed (EMA)')
    plt.axhline(kappa_avg, color='r', linestyle='--', label=f'Mean: {kappa_avg:.3f}')
    plt.xlabel('Timestep')
    plt.ylabel('Thermal Conductivity (W/m·K)')
    plt.title('NEMD Thermal Conductivity')
    plt.legend()
    plt.tight_layout()
    plt.savefig('thermal_conductivity.png', dpi=300)
    plt.show()
    
    return kappa_avg, kappa_std


# =============================================================================
# EXAMPLE 2: Plot Temperature Profiles
# =============================================================================

def example_temperature_profiles():
    """
    Visualize temperature gradients from NEMD simulation.
    """
    print("=" * 60)
    print("TEMPERATURE PROFILE VISUALIZATION")
    print("=" * 60)
    
    profiles, timesteps = readdata(PROFILE_FILE)
    
    # Get spatial coordinates (in nm for plotting)
    z_coords = Z_AXIS_LENGTH * profiles[0, :, 0] / 10  # Convert Å to nm
    
    # Plot several profiles from production period
    plt.figure(figsize=(10, 6))
    
    n_production = min(1000, len(timesteps))
    indices = np.linspace(len(timesteps) - n_production, len(timesteps) - 1, 5, dtype=int)
    
    for idx in indices:
        t_ns = timesteps[idx] / 1e6  # Convert to ns
        plt.plot(z_coords, profiles[idx, :, 1], label=f'{t_ns:.1f} ns')
    
    # Plot average
    avg_profile = np.mean(profiles[-n_production:, :, 1], axis=0)
    plt.plot(z_coords, avg_profile, 'k', linewidth=2, label='Average')
    
    plt.axhline(TEMPERATURE, color='gray', linestyle='--', alpha=0.5)
    plt.xlabel('Z-axis Distance (nm)')
    plt.ylabel('Temperature (K)')
    plt.title('Temperature Profiles from NEMD')
    plt.legend()
    plt.tight_layout()
    plt.savefig('temperature_profiles.png', dpi=300)
    plt.show()


# =============================================================================
# EXAMPLE 3: DOS Analysis and Quantum Correction
# =============================================================================

def example_dos_analysis():
    """
    Calculate density of states from velocity autocorrelation
    and apply quantum correction for heat capacity.
    
    Required files:
        - VACF output from 'compute vacf' + 'fix print'
    """
    print("=" * 60)
    print("DOS ANALYSIS AND QUANTUM CORRECTION")
    print("=" * 60)
    
    # Load VACF data
    timestep, vac_dict = read_vac_data(VACF_FILE)
    timestep *= 1000  # Convert ps to fs if needed
    
    print(f"Loaded VACF with {len(timestep)} points")
    
    # Compute FFT to get DOS
    freq, fft_vals = compute_fft(timestep, vac_dict['vac_total'])
    
    # Use only positive frequencies (first half of FFT output)
    cutoff = len(freq) // 2
    freq_pos = freq[:cutoff]
    dos = fft_vals.real[:cutoff]
    
    # Skip first point (zero frequency)
    freq_pos = freq_pos[1:]
    dos = dos[1:]
    
    # Calculate quantum correction factor
    correction = calc_quantum_correction(freq_pos, dos, TEMPERATURE)
    
    # Classical heat capacity: 3R per mole of atoms
    R = 8.314  # J/(mol·K)
    Cv_classical = 3 * R
    Cv_quantum = Cv_classical * correction
    
    print(f"\nQuantum Correction Factor: {correction:.4f}")
    print(f"Classical C_v: {Cv_classical:.2f} J/(mol·K)")
    print(f"Quantum-corrected C_v: {Cv_quantum:.2f} J/(mol·K)")
    
    # Plot VACF
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    axes[0].plot(timestep / 1000, vac_dict['vac_total'], 'k')
    axes[0].set_xlabel('Time (ps)')
    axes[0].set_ylabel('VACF')
    axes[0].set_title('Velocity Autocorrelation Function')
    
    # Plot DOS (limit frequency range for visibility)
    max_freq_plot = 200  # THz
    mask = freq_pos < max_freq_plot
    axes[1].plot(freq_pos[mask], dos[mask], 'k')
    axes[1].set_xlabel('Frequency (THz)')
    axes[1].set_ylabel('DOS (arb. units)')
    axes[1].set_title('Density of States')
    
    plt.tight_layout()
    plt.savefig('dos_analysis.png', dpi=300)
    plt.show()
    
    return correction


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    # Uncomment the analyses you want to run:
    
    # example_thermal_conductivity()
    # example_temperature_profiles()
    # example_dos_analysis()
    
    print("Update the USER CONFIGURATION section and uncomment desired analyses.")
    print("See docstrings for required files and expected data formats.")
