"""
constants.py

Contains amino acid class dictionaries and physiochemical property descriptor lists
used for encoding protein sequences in DeepPPI-style neural network models.

Author: Kiana Seraj
"""

from utils import *

# ============================== #
# Amino Acid Groupings by Property
# ============================== #

# Hydrophobicity classification (3 classes)
hydrophobicity = {
    "1": "RKEDQN",    # Hydrophilic
    "2": "GASTPHY",   # Neutral
    "3": "CLVIMFW"    # Hydrophobic
}

# Normalized Van der Waals volume
Normalized_van_der_waals_vol = {
    "1": "GASTPCD",
    "2": "NVEQIL",
    "3": "MHKFRYW"
}

# Polarity classification
polarity = {
    "1": "LIFWCMVY",   # Non-polar
    "2": "PATGS",      # Moderately polar
    "3": "HQRKNED"     # Highly polar
}

# Polarizability classification
polarizability = {
    "1": "GASDT",
    "2": "CPNVEQIL",
    "3": "KMHFRYW"
}

# Charge classification
charge = {
    "1": "KR",                      # Positive
    "2": "ANCQGHILMFPSTWYV",       # Neutral
    "3": "DE"                      # Negative
}

# Secondary structure preference
secondary_structure = {
    "1": "EALMQKRH",   # Alpha-helix
    "2": "VIYCWFT",    # Beta-sheet
    "3": "GNPSD"       # Coil
}

# Solvent accessibility
solvent_accessible = {
    "1": "ALFGCIVW",   # Buried
    "2": "PKQEND",     # Intermediate
    "3": "MPSTHY"      # Exposed
}

# Surface tension
surface_tension = {
    "1": "GQDNAHR",
    "2": "KTSEC",
    "3": "ILMFPWYV"
}

# Protein-protein interaction hot spot classification
prot_prot_hotspot = {
    "1": "DHIKNPRWY",
    "2": "EQSTGAMF",
    "3": "CLV"
}

# Protein-protein interaction propensity
prot_prot_propensity = {
    "1": "CDFMPQRWY",
    "2": "AGHVLNST",
    "3": "EIK"
}

# Protein-DNA interaction propensity (Schneider et al.)
prot_dna_propensity_Schneider = {
    "1": "GKNQRSTY",
    "2": "ADEFHILVW",
    "3": "CMP"
}

# Protein-DNA interaction propensity (Ahmad et al.)
prot_dna_propensity_Ahmad = {
    "1": "GHKNQRSTY",
    "2": "ADEFIPVW",
    "3": "CLM"
}

# Protein-RNA interaction propensity (Kim et al.)
prot_RNA_propensity_Kim = {
    "1": "HKMRY",
    "2": "FGILNPQSVW",
    "3": "CDEAT"
}

# Protein-RNA interaction propensity (Ellis et al.)
prot_RNA_propensity_Ellis = {
    "1": "HGKMRSYW",
    "2": "AFINPQT",
    "3": "CDELV"
}

# Protein-RNA interaction propensity (Phipps et al.)
prot_RNA_propensity_Phipps = {
    "1": "HKMQRS",
    "2": "ADEFGLNPVY",
    "3": "CITW"
}

# Protein-ligand interaction propensity (general)
prot_ligand_propensity = {
    "1": "CFHWY",
    "2": "GILNMSTR",
    "3": "AEDKPQV"
}

# Valid protein-ligand binding residues (experimentally confirmed)
prot_ligand_valid_propensity = {
    "1": "CFHWYM",
    "2": "DGILNSTV",
    "3": "AEKPQR"
}

# Polar protein-ligand residues
prot_ligand_polar_propensity = {
    "1": "DEHRY",
    "2": "CFKMNQSTW",
    "3": "AGILPV"
}

# Molecular weight groupings
molecular_weight = {
    "1": "AGS", 
    "2": "CDEHIKLMNQPTV", 
    "3": "FRWY"
}

# Calculated LogP (hydrophobicity)
cLogP = {
    "1": "RKDNEQH", 
    "2": "PYSTGACV", 
    "3": "WMFLI"
}

# Hydrogen bond donor classification
hydrogen_bond_donor = {
    "1": "HKNQR", 
    "2": "DESTWY", 
    "3": "ACGFILMPV"
}

# Hydrogen bond acceptor classification
hydrogen_bond_acceptor = {
    "1": "DEHNQR", 
    "2": "KSTWY", 
    "3": "ACGFILMPV"
}

# Water solubility groupings
water_Solubility = {
    "1": "ACGKRT", 
    "2": "EFHILMNPQSVW", 
    "3": "DY"
}

# Amino acid flexibility grouping
Amino_acid_flexibility = {
    "1": "EGKNQS", 
    "2": "ADHIPRTV", 
    "3": "CFLMWY"
}

# ============================== #
# List of Physiochemical Descriptors (functions from utils.py)
# ============================== #

physiochemical_properties = [
    hydrophobicity_descriptor,
    Normalized_van_der_waals_vol_descriptor,
    polarity_descriptor,
    polarizability_descriptor,
    charge_descriptor,
    secondary_structure_descriptor,
    solvent_accessible_descriptor,
    surface_tension_descriptor,
    prot_prot_hotspot_descriptor,
    prot_prot_propensity_descriptor,
    prot_dna_propensity_Schneider_descriptor,
    prot_dna_propensity_Ahmad_descriptor,
    prot_RNA_propensity_Kim_descriptor,
    prot_RNA_propensity_Ellis_descriptor,
    prot_RNA_propensity_Phipps_descriptor,
    prot_ligand_propensity_descriptor,
    prot_ligand_valid_propensity_descriptor,
    prot_ligand_polar_propensity_descriptor,
    molecular_weight_descriptor,
    cLogP_descriptor,
    hydrogen_bond_donor_descriptor,
    hydrogen_bond_acceptor_descriptor,
    water_Solubility_descriptor,
    Amino_acid_flexibility_descriptor  
]

# ============================== #
# Grantham Distance Matrix
# ============================== #

# Chemical distance between amino acid pairs based on physicochemical properties
# Source: Grantham R. Science (1974)
# Only upper triangular entries stored to avoid redundancy
grantham_distances = {
    ('A', 'C'): 195, ('A', 'D'): 126, ('A', 'E'): 107, ('A', 'F'): 113, ('A', 'G'): 60,
    ...
    ('V', 'W'): 88, ('V', 'Y'): 55,
    ('W', 'Y'): 37
}
