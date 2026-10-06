# DNA Hybridization Thermodynamics Calculator

A biophysical tool designed to model the coil-globule (melting) transition of short DNA oligonucleotides, with a specific focus on highly degraded ancient DNA (aDNA) fragments.

Building upon previous work in modeling stochastic polymer phase transitions (Markov Chains and Monte Carlo simulations), this project applies rigorous statistical mechanics to biological systems by implementing the **Nearest-Neighbor thermodynamic model** (SantaLucia, 1998).

## Biophysical Features

* **Thermodynamic Profiling:** Calculates the total Enthalpy ($\Delta H$) and Entropy ($\Delta S$) of a DNA duplex by evaluating base-stacking interactions between adjacent nucleotides and incorporating sequence-dependent end effects.
* **Phase Transition (Melting) Estimation:** Computes the melting temperature ($T_m$) at which the DNA duplex separates into single strands.
* **Electrostatic Corrections:** Models the counterion condensation effect (Owczarzy et al., 2004) by applying salt-concentration corrections. This simulates how the ionic strength of different buffers screens the electrostatic repulsion of the phosphate backbone, directly dictating duplex thermodynamic stability.

## Research Context

This calculator serves as a foundational theoretical framework for understanding and optimizing the physical parameters of **hybridization capture** protocols. By fine-tuning the thermodynamic and electrostatic environment, it is possible to maximize physical affinity and target enrichment efficiency when recovering highly degraded Pleistocene DNA from complex archaeological matrices.

## Usage

Run the `dna_thermo_kinetics.py` script to simulate how changing the sodium ion concentration ($[Na^+]$) in the buffer alters the thermal stability of a sample DNA fragment.
