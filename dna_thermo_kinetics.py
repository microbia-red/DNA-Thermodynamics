import math

class DNAThermodynamics:
    """
    Thermodynamic calculator for DNA oligonucleotides (e.g., highly degraded ancient DNA)
    based on the Nearest-Neighbor model (SantaLucia, 1998).
    """
    def __init__(self):
        # Nearest-Neighbor thermodynamic parameters (dH in kcal/mol, dS in cal/K*mol)
        # Unified parameters based on SantaLucia 1998, Table 2
        self.nn_params = {
            'AA': (-7.9, -22.2), 'TT': (-7.9, -22.2),
            'AT': (-7.2, -20.4), 'TA': (-7.2, -21.3),
            'CA': (-8.5, -22.7), 'TG': (-8.5, -22.7),
            'GT': (-8.4, -22.4), 'AC': (-8.4, -22.4),
            'CT': (-7.8, -21.0), 'AG': (-7.8, -21.0),
            'GA': (-8.2, -22.2), 'TC': (-8.2, -22.2),
            'CG': (-10.6, -27.2), 'GC': (-9.8, -24.4),
            'GG': (-8.0, -19.9), 'CC': (-8.0, -19.9)
        }
        self.R = 1.987 # Ideal gas constant in cal/(K*mol)

    def calculate_thermo(self, sequence):
        """Calculates the total enthalpy (dH) and entropy (dS) of the sequence."""
        sequence = sequence.upper()
        total_dH = 0.0
        total_dS = 0.0
        
        # 1. Initiation parameters depending on terminal base pairs (SantaLucia 1998, Table 2)
        # We check the first and the last nucleotide of the sequence to account for end effects
        terminals = [sequence[0], sequence[-1]]
        for base in terminals:
            if base in ['G', 'C']:
                total_dH += 0.1   # Init w/term G-C
                total_dS += -2.8
            elif base in ['A', 'T']:
                total_dH += 2.3   # Init w/term A-T
                total_dS += 4.1
                
        # 2. Nearest-Neighbor propagation parameters
        for i in range(len(sequence) - 1):
            pair = sequence[i:i+2]
            if pair in self.nn_params:
                dH, dS = self.nn_params[pair]
                total_dH += dH
                total_dS += dS
                
        # Convert dH to cal/mol to match dS units
        return total_dH * 1000, total_dS

    def melting_temperature(self, sequence, dna_conc=1e-6, na_conc=0.05):
        """
        Calculates the melting temperature (Tm) in degrees Celsius.
        dna_conc: DNA concentration in Molar (default 1 uM)
        na_conc: Sodium (Na+) concentration in Molar (for salt correction)
        """
        dH, dS = self.calculate_thermo(sequence)
        
        # Base Tm calculation in Kelvin
        # Formula: Tm = dH / (dS + R * ln(C/4))
        tm_kelvin = dH / (dS + self.R * math.log(dna_conc / 4.0))
        
        # Salt concentration correction (Schildkraut-Lifson / SantaLucia / Owczarzy)
        # Buffer ionic strength affects the electrostatic repulsion of the polymer backbone
        tm_celsius = (tm_kelvin - 273.15) + 16.6 * math.log10(na_conc)
        
        return tm_celsius

# USAGE EXAMPLE
if __name__ == "__main__":
    thermo = DNAThermodynamics()
    
    # Example of DNA fragment (e.g., 35 base pairs)
    dna_fragment = "GATTACAGATCGATCGATCGATCGATCGATCGATC" 
    
    # Buffers with different ionic strengths (Sodium concentrations)
    buffers = [0.01, 0.05, 0.1, 0.5] # Molar
    
    print(f"aDNA Sequence ({len(dna_fragment)} bp): {dna_fragment}")
    print("-" * 50)
    print("Influence of buffer ionic strength (Na+) on thermal stability:")
    
    for na in buffers:
        tm = thermo.melting_temperature(dna_fragment, na_conc=na)
        print(f"[Na+] = {na:4.2f} M  --> Tm = {tm:.2f} °C")