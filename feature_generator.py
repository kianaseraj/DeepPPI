# -*- coding: utf-8 -*-
"""
feature_generator.py

Generates protein features from amino acid sequences using multiple descriptors:
- AAC, DPC
- Physicochemical encodings (CTD)
- Quasi sequence order (QSOrder)
- Sequence order coupling (SOCNumber)
- APAAC (normalized)

Outputs features as tensors using PyTorch.

Author: Kiana Seraj
"""

import os
import numpy as np
import torch
import protpy
from constants import *
from utils import AAC, DPC, composition_descriptor, transition_descriptor, distribution_descriptor

# =========================
# Utility: Min-Max Normalization
# =========================

def norm(x):
    """
    Normalize to range [0, 1].
    """
    min_val, max_val = min(x), max(x)
    return [(i - min_val) / (max_val - min_val) for i in x]

def norm_(x):
    """
    Normalize to range [-1, 1].
    """
    min_val, max_val = min(x), max(x)
    return [2 * ((i - min_val) / (max_val - min_val)) - 1 for i in x]


# =========================
# CTD Descriptor Generator
# =========================

def CTD(sequence):
    """
    Computes CTD (Composition, Transition, Distribution) features
    for all predefined physicochemical properties.

    Returns:
        list: 504-length CTD vector.
    """
    ctd = []
    for descriptor_func in physiochemical_properties:
        encoded = descriptor_func(sequence)
        ctd.extend(composition_descriptor(encoded))
        ctd.extend(distribution_descriptor(encoded))
        ctd.extend(transition_descriptor(encoded))
    return ctd


# =========================
# Load QSOrder and SOCNumber (from precomputed txt files)
# =========================

def QSOrder(protname, feature_dir):
    """
    Load QSOrder descriptor from file.

    Args:
        protname (str): Protein ID.
        feature_dir (str): Directory containing QSOrder output txt.

    Returns:
        list: QSOrder descriptor.
    """
    path = os.path.join(feature_dir, f"{protname}_QSOrder.txt")
    with open(path) as f:
        values = f.read().strip().split("\t")[1:]
    return [float(v) for v in values]


def SOCNumber(protname, feature_dir):
    """
    Load Sequence Order Coupling Number (SOCNumber) from file and normalize Grantham part.

    Args:
        protname (str): Protein ID.
        feature_dir (str): Directory containing SOCNumber output txt.

    Returns:
        tuple: (Schneider SOC [30], normalized Grantham SOC [30])
    """
    path = os.path.join(feature_dir, f"{protname}_SOCNumber.txt")
    with open(path) as f:
        values = f.read().strip().split("\t")[1:]
    soc = [float(v) for v in values]
    return soc[:30], norm(soc[30:])


# =========================
# APAAC Feature Generator
# =========================

def APAAC(sequence):
    """
    Amphiphilic Pseudo Amino Acid Composition with normalization.

    Args:
        sequence (str): Amino acid sequence.

    Returns:
        list: Normalized APAAC features.
    """
    apaac = list(protpy.amphiphilic_pseudo_amino_acid_composition(sequence, lamda=30, weight=0.5).values[0])
    return norm_(apaac)


# =========================
# Main Feature Extraction Pipeline
# =========================

def extract_features(sequence_dir, feature_dir, output_dir):
    """
    Extract features for all sequences in a directory.

    Args:
        sequence_dir (str): Path to directory containing .npy files with protein sequences.
        feature_dir (str): Path to directory containing QSOrder and SOCNumber files.
        output_dir (str): Path to directory to save torch tensor features.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    for filename in os.listdir(sequence_dir):
        if not filename.endswith(".npy"):
            continue
        
        filepath = os.path.join(sequence_dir, filename)
        sequence = np.load(filepath, allow_pickle=True).item()
        protname = filename.split(".")[0]

        aac = AAC(sequence)
        dpc = DPC(sequence)
        ctd = CTD(sequence)
        qsorder = QSOrder(protname, feature_dir)
        soc_schneider, soc_grantham = SOCNumber(protname, feature_dir)
        apaac = APAAC(sequence)

        features = aac + dpc + ctd + qsorder + soc_schneider + soc_grantham + apaac
        torch.save(features, os.path.join(output_dir, f"{protname}.pt"))


# =========================
# Optional Run Block
# =========================

if __name__ == "__main__":
    # Example usage (edit these paths as needed)
    extract_features(
        sequence_dir="PATH/to/sequence_directory",
        feature_dir="PATH/to/feature_extraction_output",
        output_dir="PATH/to/save_tensor_features"
    )
