# -*- coding: utf-8 -*-
"""
model.py

Fully-connected neural network for binary classification of protein-protein interactions (PPI).
Processes two protein embeddings in parallel, then merges them for final prediction.

Author: Kiana Seraj
"""

import torch
import torch.nn as nn
import torch.nn.functional as F


class FC(nn.Module):
    """
    Fully connected neural network for PPI classification with two protein inputs.

    Architecture:
        - Each protein: 3 FC + BatchNorm + LeakyReLU + Dropout layers
        - Combined: FC + LeakyReLU + Dropout + Output

    Output:
        - 2D tensor (binary classification via Sigmoid)
    """

    def __init__(self):
        super(FC, self).__init__()

        # Input size is 1164 for each protein (based on features)

        # Protein 1
        self.pro1_fc1 = nn.Linear(1164, 512)
        self.pro1_bn1 = nn.BatchNorm1d(512)
        self.pro1_fc2 = nn.Linear(512, 256)
        self.pro1_bn2 = nn.BatchNorm1d(256)
        self.pro1_fc3 = nn.Linear(256, 128)
        self.pro1_bn3 = nn.BatchNorm1d(128)

        # Protein 2
        self.pro2_fc1 = nn.Linear(1164, 512)
        self.pro2_bn1 = nn.BatchNorm1d(512)
        self.pro2_fc2 = nn.Linear(512, 256)
        self.pro2_bn2 = nn.BatchNorm1d(256)
        self.pro2_fc3 = nn.Linear(256, 128)
        self.pro2_bn3 = nn.BatchNorm1d(128)

        # Combined layer
        self.fc1 = nn.Linear(256, 128)
        self.out = nn.Linear(128, 2)

        # Activation & dropout
        self.relu = nn.LeakyReLU()
        self.sigmoid = nn.Sigmoid()
        self.dropout = nn.Dropout(0.2)

    def forward(self, pro1_data, pro2_data):
        """
        Forward pass for two protein input embeddings.

        Args:
            pro1_data (Tensor): Feature vector for protein 1.
            pro2_data (Tensor): Feature vector for protein 2.

        Returns:
            Tensor: Output prediction of shape (batch_size, 2)
        """
        # Protein 1 path
        x1 = self.pro1_fc1(pro1_data)
        x1 = self.pro1_bn1(x1)
        x1 = self.relu(x1)
        x1 = self.dropout(x1)

        x1 = self.pro1_fc2(x1)
        x1 = self.pro1_bn2(x1)
        x1 = self.relu(x1)
        x1 = self.dropout(x1)

        x1 = self.pro1_fc3(x1)
        x1 = self.pro1_bn3(x1)
        x1 = self.relu(x1)
        x1 = self.dropout(x1)

        # Protein 2 path
        x2 = self.pro2_fc1(pro2_data)
        x2 = self.pro2_bn1(x2)
        x2 = self.relu(x2)
        x2 = self.dropout(x2)

        x2 = self.pro2_fc2(x2)
        x2 = self.pro2_bn2(x2)
        x2 = self.relu(x2)
        x2 = self.dropout(x2)

        x2 = self.pro2_fc3(x2)
        x2 = self.pro2_bn3(x2)
        x2 = self.relu(x2)
        x2 = self.dropout(x2)

        # Merge
        x = torch.cat((x1, x2), dim=1)
        x = self.fc1(x)
        x = self.relu(x)
        x = self.dropout(x)
        x = self.out(x)
        return self.sigmoid(x)
