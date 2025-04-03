"""
utils.py

Utility functions for computing protein sequence descriptors including:
- Amino acid composition
- Physicochemical property encodings
- Sequence order and transition metrics
- Quasi sequence order, pseudo amino acid composition
- UniProt sequence fetching

Author: Kiana Seraj
"""

from constants import *
import math
import numpy as np
import protpy
import requests
from protFeat.feature_extracter import extract_protein_feature

# =========================
# Feature 1: Basic Composition Descriptors
# =========================

def AAC(sequence):
    """
    Amino Acid Composition (AAC), 20 features.

    Args:
        sequence (str): Amino acid sequence.

    Returns:
        list: Normalized AAC features.
    """
    amino_acid_composition = protpy.amino_acid_composition(sequence)
    return [i / 100 for i in amino_acid_composition.values[0]]


def DPC(sequence):
    """
    Dipeptide Composition (DPC), 400 features.

    Args:
        sequence (str): Amino acid sequence.

    Returns:
        list: Normalized DPC features.
    """
    dipeptide_composition = protpy.dipeptide_composition(sequence)
    return [i / 100 for i in dipeptide_composition.values[0]]


def APAAC(sequence):
    """
    Amphiphilic Pseudo Amino Acid Composition, 80 features.

    Args:
        sequence (str): Amino acid sequence.

    Returns:
        list: APAAC features.
    """
    return protpy.amphiphilic_pseudo_amino_acid_composition(sequence, lamda=30, weight=0.5)


# =========================
# Feature 2: Physicochemical Property Encodings (Each returns a 1/2/3-encoded string)
# =========================

def _property_descriptor(sequence, prop_dict):
    return "".join(key for aa in sequence for key, val in prop_dict.items() if aa in val)

# Automatically generated wrappers
hydrophobicity_descriptor = lambda seq: _property_descriptor(seq, hydrophobicity)
Normalized_van_der_waals_vol_descriptor = lambda seq: _property_descriptor(seq, Normalized_van_der_waals_vol)
polarity_descriptor = lambda seq: _property_descriptor(seq, polarity)
polarizability_descriptor = lambda seq: _property_descriptor(seq, polarizability)
charge_descriptor = lambda seq: _property_descriptor(seq, charge)
secondary_structure_descriptor = lambda seq: _property_descriptor(seq, secondary_structure)
solvent_accessible_descriptor = lambda seq: _property_descriptor(seq, solvent_accessible)
surface_tension_descriptor = lambda seq: _property_descriptor(seq, surface_tension)
prot_prot_hotspot_descriptor = lambda seq: _property_descriptor(seq, prot_prot_hotspot)
prot_prot_propensity_descriptor = lambda seq: _property_descriptor(seq, prot_prot_propensity)
prot_dna_propensity_Schneider_descriptor = lambda seq: _property_descriptor(seq, prot_dna_propensity_Schneider)
prot_dna_propensity_Ahmad_descriptor = lambda seq: _property_descriptor(seq, prot_dna_propensity_Ahmad)
prot_RNA_propensity_Kim_descriptor = lambda seq: _property_descriptor(seq, prot_RNA_propensity_Kim)
prot_RNA_propensity_Ellis_descriptor = lambda seq: _property_descriptor(seq, prot_RNA_propensity_Ellis)
prot_RNA_propensity_Phipps_descriptor = lambda seq: _property_descriptor(seq, prot_RNA_propensity_Phipps)
prot_ligand_propensity_descriptor = lambda seq: _property_descriptor(seq, prot_ligand_propensity)
prot_ligand_valid_propensity_descriptor = lambda seq: _property_descriptor(seq, prot_ligand_valid_propensity)
prot_ligand_polar_propensity_descriptor = lambda seq: _property_descriptor(seq, prot_ligand_polar_propensity)
molecular_weight_descriptor = lambda seq: _property_descriptor(seq, molecular_weight)
cLogP_descriptor = lambda seq: _property_descriptor(seq, cLogP)
hydrogen_bond_donor_descriptor = lambda seq: _property_descriptor(seq, hydrogen_bond_donor)
hydrogen_bond_acceptor_descriptor = lambda seq: _property_descriptor(seq, hydrogen_bond_acceptor)
water_Solubility_descriptor = lambda seq: _property_descriptor(seq, water_Solubility)
Amino_acid_flexibility_descriptor = lambda seq: _property_descriptor(seq, Amino_acid_flexibility)


# =========================
# Feature 3: CTD Encoding (Composition, Transition, Distribution)
# =========================

def composition_descriptor(encoding):
    """
    Composition descriptor: frequency of classes 1, 2, 3 in encoded string.

    Args:
        encoding (str): Encoded property string (e.g., '1213231...').

    Returns:
        list: Composition feature vector (length 3).
    """
    return [encoding.count("1")/len(encoding), encoding.count("2")/len(encoding), encoding.count("3")/len(encoding)]


def transition_descriptor(encoding):
    """
    Transition descriptor: frequency of transitions between classes (1-2, 1-3, 2-3).

    Args:
        encoding (str): Encoded property string.

    Returns:
        list: Transition frequencies (length 3).
    """
    transitions_ = {}
    transitions = ["12", "13", "21", "23", "31", "32"]
    for i in range(len(encoding) - 1):
        transition = encoding[i:i+2]
        if transition in transitions:
            transitions_[transition] = transitions_.get(transition, 0) + 1

    combined = {"12": 0, "13": 0, "23": 0}
    for t, count in transitions_.items():
        if t in ["12", "21"]: combined["12"] += count
        elif t in ["13", "31"]: combined["13"] += count
        else: combined["23"] += count

    denom = len(encoding) - 1 if len(encoding) > 1 else 1
    return [combined["12"]/denom, combined["13"]/denom, combined["23"]/denom]


def calculate_class_lengths(encoding, target_class):
    """
    Internal function for distribution descriptor.

    Returns the position (1-based) of the first, 25%, 50%, 75%, and 100% occurrence of a class.
    """
    positions = [i for i, char in enumerate(encoding) if char == target_class]
    if not positions:
        return [0] * 5

    n = len(positions)
    return [
        positions[0]+1,
        positions[min(math.ceil(n*0.25)-1, n-1)]+1,
        positions[min(math.ceil(n*0.50)-1, n-1)]+1,
        positions[min(math.ceil(n*0.75)-1, n-1)]+1,
        positions[-1]+1
    ]


def distribution_descriptor(encoding):
    """
    Distribution descriptor: relative positions of key percentiles for each class.

    Args:
        encoding (str): Encoded string.

    Returns:
        list: 15 distribution features.
    """
    features = []
    for c in ["1", "2", "3"]:
        raw = calculate_class_lengths(encoding, c)
        features.extend([r / len(encoding) for r in raw])
    return features


# =========================
# Feature 4: Sequence Order
# =========================

def QSOD(input_folder, fasta_file_name):
    """
    Quasi Sequence Order Descriptor (QSOD), 100 dimensions.

    Args:
        input_folder (str): Path to directory containing FASTA.
        fasta_file_name (str): File name of the FASTA.

    Returns:
        list: QSOD features.
    """
    return extract_protein_feature("QSOrder", input_folder, fasta_file_name, place_protein_id=0)


def SOCD(input_folder, fasta_file_name):
    """
    Sequence Order Coupling Descriptor (SOCD), 60 dimensions.

    Args:
        input_folder (str): Path to directory containing FASTA.
        fasta_file_name (str): File name of the FASTA.

    Returns:
        list: SOCD features.
    """
    return extract_protein_feature("SOCNumber", input_folder, fasta_file_name, place_protein_id=0)


# =========================
# Extra: Fetching Sequences from UniProt
# =========================

def fetch_sequence(uniprot_id):
    """
    Fetch a protein sequence from UniProt using its UniProt ID.

    Args:
        uniprot_id (str): e.g., 'P12345'

    Returns:
        str or None: Protein sequence if found, else None.
    """
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        lines = response.text.split('\n')
        return ''.join(lines[1:])
    else:
        print(f"[Error] Failed to fetch sequence for UniProt ID: {uniprot_id}")
        return None


# Example usage (can be moved to a separate script or __main__ guard)
if __name__ == "__main__":
    example_id = "P12345"
    sequence = fetch_sequence(example_id)
    if sequence:
        print(f"Sequence for UniProt ID {example_id}:\n{sequence}")
    else:
        print("No sequence found.")
