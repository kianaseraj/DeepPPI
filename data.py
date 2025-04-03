# -*- coding: utf-8 -*-
"""
data.py

Creates a PyTorch Dataset and DataLoader for protein-protein interaction (PPI) modeling.
Each protein is represented by a precomputed feature vector stored as a `.pt` tensor.

Author: Kiana Seraj
"""

import os
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader

# =========================
# Dataset Class
# =========================

class ProtDataset(Dataset):
    """
    Custom Dataset for protein-protein interaction pairs.

    Each row in the input data file should contain:
        [protein_id_1, protein_id_2, label]
    The corresponding features must exist in `feature_dir` as:
        protein_id.pt

    Args:
        feature_dir (str): Directory with protein feature .pt files.
        data_file (np.array): Array with protein pairs and labels.
    """

    def __init__(self, feature_dir: str, data_file: np.ndarray):
        self.feature_dir = feature_dir
        self.prot_1 = data_file[:, 0]
        self.prot_2 = data_file[:, 1]
        self.labels = data_file[:, 2].astype(float).astype(int)
        self.n_samples = data_file.shape[0]

    def __getitem__(self, index):
        """
        Fetch a single sample: (feature1, feature2, label).
        Returns torch tensors.
        """
        prot1_path = os.path.join(self.feature_dir, f"{self.prot_1[index]}.pt")
        prot2_path = os.path.join(self.feature_dir, f"{self.prot_2[index]}.pt")

        prot1 = torch.load(prot1_path)
        prot2 = torch.load(prot2_path)
        label = torch.tensor(self.labels[index])

        return prot1, prot2, label

    def __len__(self):
        return self.n_samples


# =========================
# Main (Example Usage)
# =========================

if __name__ == "__main__":
    # Example placeholder paths (replace with real ones or use argparse)
    train_data_path = "PATH/to/training_dataset.npy"
    val_data_path = "PATH/to/validation_dataset.npy"
    feature_dir = "PATH/to/feature_vectors"

    # Load training/validation pair-label arrays
    train_data = np.load(train_data_path)
    val_data = np.load(val_data_path)

    # Create Dataset objects
    train_dataset = ProtDataset(feature_dir, train_data)
    val_dataset = ProtDataset(feature_dir, val_data)

    # Create DataLoaders
    train_loader = DataLoader(train_dataset, batch_size=256, shuffle=True, num_workers=0)
    val_loader = DataLoader(val_dataset, batch_size=256, shuffle=False, num_workers=0)

    print("Train samples:", len(train_dataset))
    print("Val samples:", len(val_dataset))
