#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis Functions for Polyurethane Thermal Properties Simulation

This module contains reusable functions for post-processing LAMMPS simulation
output, including:
- NEMD thermal conductivity calculation (Müller-Plathe method)
- Velocity autocorrelation and DOS analysis
- Quantum correction factors for heat capacity
- LAMMPS log file parsing

Author: Tristan Myers
Created: May 2025
"""

import numpy as np
import scipy.constants
import re


# =============================================================================
# LAMMPS Output Parsing Functions
# =============================================================================

def readdata(filename):
    """
    Parse LAMMPS ave/chunk output files (e.g., temperature profiles from fix ave/chunk).
    
    Parameters
    ----------
    filename : str
        Path to LAMMPS ave/chunk output file (e.g., profile.mp)
    
    Returns
    -------
    datapoint : np.ndarray
        3D array of shape (n_timesteps, n_bins, 2) containing [coord, value] pairs
    timesteps : list
        List of timestep values
    """
    with open(filename, 'r') as inputdata:
        lines = inputdata.readlines()

    all_data = []
    timesteps = []
    block_line = 0
    npoint = 0

    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue

        block_line += 1

        if block_line == 1:
            val = line.split()
            timestep = int(val[0])
            npoint = int(val[1])
            timesteps.append(timestep)
            current_timestep_data = []
        elif 1 < block_line <= 1 + npoint:
            val = re.findall(r'\S+', line)
            coord = float(val[1])
            value = float(val[3])
            current_timestep_data.append([coord, value])
        if block_line == npoint + 1:
            all_data.append(current_timestep_data)
            block_line = 0

    datapoint = np.array(all_data, dtype=float)
    return datapoint, timesteps


def readfile_log(path, property_indices, Nevery=1):
    """
    Parse LAMMPS log files to extract thermodynamic properties.
    
    Searches for simulation output blocks (starting with "Step" header and 
    ending with "Loop time") and extracts specified properties.
    
    Parameters
    ----------
    path : str
        Path to LAMMPS log file
    property_indices : array-like
        Column indices to extract (0-indexed, excluding Step column)
    Nevery : int, optional
        Only include every Nth line (default: 1)
    
    Returns
    -------
    prop_list : np.ndarray
        2D array with columns [timestep, property1, property2, ...]
    """
    with open(path, "r") as file:
        lines = file.readlines()
    
    prop_list = np.array([])
    iterator = 0
    sim_printout_switch = False
    
    for i in range(len(lines)):
        if len(lines[i].split()) == 0:
            continue
        if lines[i].split()[0] == "Step":
            sim_printout_switch = True
            continue
        if sim_printout_switch and lines[i].split()[0] == "Loop" and lines[i].split()[1] == "time":
            sim_printout_switch = False
        if i == len(lines) - 1:
            break
        if sim_printout_switch:
            iterator += 1
            if iterator % Nevery == 0:
                split_line = lines[i].split()
                newline = np.array([int(split_line[0])])
                for j in range(len(property_indices)):
                    newline = np.append(newline, float(split_line[property_indices[j]]))
                prop_list = np.append(prop_list, newline, axis=0)
    
    prop_list = np.resize(prop_list, (int(len(prop_list)/(len(property_indices)+1)), len(property_indices)+1))
    return prop_list


def read_vac_data(filename):
    """
    Parse LAMMPS velocity autocorrelation output file.
    
    Expected format: time vac_x vac_y vac_z vac_total
    
    Parameters
    ----------
    filename : str
        Path to VACF output file (e.g., vacf.txt)
    
    Returns
    -------
    timestep : np.ndarray
        Time values
    vac_dict : dict
        Dictionary with keys 'vac_x', 'vac_y', 'vac_z', 'vac_total'
    """
    data = []
    with open(filename, 'r') as file:
        for line in file:
            if line.strip() == "" or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) == 5:
                data.append([float(p) for p in parts])
    
    data = np.array(data)
    timestep = data[:, 0]
    vac_dict = {
        'vac_x': data[:, 1],
        'vac_y': data[:, 2],
        'vac_z': data[:, 3],
        'vac_total': data[:, 4]
    }
    
    return timestep, vac_dict


# =============================================================================
# NEMD Thermal Conductivity Calculation
# =============================================================================

def get_thermcond_NEMD(deltaE_MP, MP_profile, delta_t, z_axis_len, Xsectional_len, avg_bool=False):
    """
    Calculate thermal conductivity from Müller-Plathe NEMD simulation.
    
    Uses the cumulative energy transferred by fix thermal/conductivity and
    the resulting temperature gradient to compute thermal conductivity.
    
    Parameters
    ----------
    deltaE_MP : array-like
        Cumulative energy transferred (kcal/mol) from fix thermal/conductivity
    MP_profile : tuple
        (temperature_profiles, timesteps) from ave/chunk output
        - temperature_profiles: 3D array (n_timesteps, n_bins, 2)
        - timesteps: list of timestep values
    delta_t : float
        Sampling interval in fs (femtoseconds)
    z_axis_len : float
        Length of simulation box in z-direction (Angstroms)
    Xsectional_len : float
        Cross-sectional length in x and y (Angstroms), assuming square cross-section
    avg_bool : bool, optional
        If True, return (average, std_dev); if False, return full profile
    
    Returns
    -------
    kappa : float or np.ndarray
        Thermal conductivity in W/(m·K)
    """
    tprofiles_data, timesteps = MP_profile[0], MP_profile[1]
    
    deltaE_MP = np.array(deltaE_MP)
    delta_KE = np.gradient(deltaE_MP) / delta_t  # Energy transfer rate, (kcal/mol)/fs
    
    temp_profiles = tprofiles_data[:, :, 1]
    axis_spacing = z_axis_len * tprofiles_data[0, :, 0]
    
    n_timesteps = len(timesteps)
    dtdz_arr = np.zeros(n_timesteps)  # K/Å
    
    for i in range(n_timesteps):
        # Fit linear slopes to hot and cold regions
        slope1 = np.polyfit(axis_spacing[:11], temp_profiles[i, :11], 1)[0]
        slope2 = np.polyfit(axis_spacing[11:], temp_profiles[i, 11:], 1)[0]
        
        # Weight by bin count (11 bins on one side, 9 on other for 20 total)
        dtdz_arr[i] = slope1 * (11 / len(axis_spacing)) - slope2 * (9 / len(axis_spacing))
    
    # Heat flux: 0.5 factor accounts for bidirectional heat flow
    flux_profile = 0.5 * delta_KE / Xsectional_len**2  # kcal/(mol·Å²·fs)
    
    # Thermal conductivity
    kappa_profile = flux_profile / dtdz_arr  # kcal/(mol·Å·fs·K)
    
    # Convert to SI units: W/(m·K)
    # 1 kcal/mol = 4184 J / N_A
    # 1 Å = 1e-10 m
    # 1 fs = 1e-15 s
    kappa_conv = (4184) / scipy.constants.N_A * (1e15) * (1e10)
    kappa_profile_SI = kappa_profile * kappa_conv
    
    if avg_bool:
        return np.average(kappa_profile_SI), np.std(kappa_profile_SI)
    else:
        return kappa_profile_SI


# =============================================================================
# DOS and Quantum Correction Functions
# =============================================================================

def compute_fft(timestep, signal):
    """
    Compute FFT of a time-domain signal.
    
    Parameters
    ----------
    timestep : np.ndarray
        Time values (determines sampling frequency)
    signal : np.ndarray
        Signal values (e.g., velocity autocorrelation)
    
    Returns
    -------
    freq : np.ndarray
        Frequency values in inverse time units (e.g., THz if timestep in ps)
    fft_vals : np.ndarray
        Complex FFT values
    """
    dt = timestep[1] - timestep[0]
    freq = np.fft.fftfreq(len(signal), d=dt)
    fft_vals = np.fft.fft(signal)
    return freq, fft_vals


def calc_quantum_correction(freq, dos, T):
    """
    Calculate quantum correction factor for heat capacity from DOS.
    
    Applies the Bose-Einstein quantum correction to convert classical MD
    heat capacity to quantum-corrected values.
    
    Parameters
    ----------
    freq : np.ndarray
        Frequency array in THz (ps⁻¹), as returned by np.fft.fftfreq
        Note: Only positive frequencies should be used
    dos : np.ndarray
        Density of states (real-valued, from FFT of VACF)
    T : float
        Temperature in Kelvin
    
    Returns
    -------
    float
        Quantum correction factor. Multiply by 3R (24.94 J/(mol·K)) 
        to get molar heat capacity
    
    Notes
    -----
    The dimensionless parameter u = ℏω/(kB·T) where ω = 2πf.
    The weighting function is u²exp(u)/(exp(u)-1)².
    
    For classical limit (u→0), correction factor → 1.
    For quantum limit (u→∞), correction factor → 0.
    """
    # Convert ordinary frequency (THz) to angular frequency (rad/s)
    omega = 2 * np.pi * freq * 1e12  # rad/s
    
    # Bose-Einstein parameter: u = ℏω/(kB·T)
    u = scipy.constants.hbar * omega / (scipy.constants.k * T)
    
    # Handle u=0 to avoid division by zero
    with np.errstate(divide='ignore', invalid='ignore'):
        weight = np.where(
            np.abs(u) > 1e-10,
            u**2 * np.exp(u) / (np.exp(u) - 1)**2,
            1.0  # Classical limit
        )
    
    # Integrate: correction = ∫ DOS(ω)·weight(ω) dω / ∫ DOS(ω) dω
    numerator = np.trapz(dos * weight, freq)
    denominator = np.trapz(dos, freq)
    
    return numerator / denominator


# =============================================================================
# Utility Functions
# =============================================================================

def EMA_series(prop_list, alpha):
    """
    Calculate exponential moving average of a time series.
    
    Parameters
    ----------
    prop_list : array-like
        Input time series data
    alpha : float
        Smoothing parameter (0 < alpha < 1)
        - Small alpha (e.g., 0.05): more smoothing, slow response
        - Large alpha (e.g., 0.3): less smoothing, fast response
    
    Returns
    -------
    EMA_profile : np.ndarray
        Exponentially smoothed time series
    """
    series_len = len(prop_list)
    EMA_profile = np.zeros(series_len)
    
    for i in range(series_len):
        if i == 0:
            EMA_profile[0] = prop_list[0]
        else:
            EMA_profile[i] = alpha * prop_list[i] + (1 - alpha) * EMA_profile[i-1]
    
    return EMA_profile
