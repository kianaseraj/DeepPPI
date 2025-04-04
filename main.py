# -*- coding: utf-8 -*-
"""
main.py

Main training script for protein-protein interaction prediction using a fully connected NN.

Author: Kiana Seraj
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.optim.lr_scheduler import MultiStepLR
from tqdm import tqdm
import numpy as np

from src.metrics import get_accuracy, get_mse
from src.model import FC
from src.data import train_loader, val_loader


# Set device
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print("Device:", device)
print("Data size:", len(train_loader), "train batches,", len(val_loader), "val batches")


# ========================
# Training Function
# ========================
def training(model, train_loader, device, optimizer, scheduler):
    model.train()
    loss_func = nn.BCELoss()
    predictions_tr = torch.Tensor()
    labels_tr = torch.Tensor()

    for prot1, prot2, label in train_loader:
        prot1, prot2, label = prot1.to(device), prot2.to(device), label.view(-1, 1).float().to(device)

        optimizer.zero_grad()
        output = model(prot1, prot2)[:, 1].unsqueeze(1)  # Get probability of class 1

        loss = loss_func(output, label)
        loss.backward()
        optimizer.step()

        predictions_tr = torch.cat((predictions_tr, output.cpu()), 0)
        labels_tr = torch.cat((labels_tr, label.cpu()), 0)

    scheduler.step()

    labels_tr = labels_tr.detach().numpy()
    predictions_tr = predictions_tr.detach().numpy()
    acc_tr = get_accuracy(labels_tr, predictions_tr, 0.5)
    print(f"[Train] Loss: {loss.item():.4f} | Accuracy: {acc_tr:.2f}%")


# ========================
# Validation Function
# ========================
def validation(model, val_loader, device):
    model.eval()
    predictions = torch.Tensor()
    labels = torch.Tensor()

    with torch.no_grad():
        for prot1, prot2, label in val_loader:
            prot1, prot2 = prot1.to(device), prot2.to(device)
            output = model(prot1, prot2)[:, 1].unsqueeze(1)  # Take probability for class 1
            predictions = torch.cat((predictions, output.cpu()), 0)
            labels = torch.cat((labels, label.view(-1,1).cpu()), 0)

    return labels.numpy(), predictions.numpy()


# ========================
# Main Training Loop
# ========================
model = FC().to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
scheduler = MultiStepLR(optimizer, milestones=[1, 5], gamma=0.5)

num_epochs = 30

for epoch in range(num_epochs):
    print(f"\n=== Epoch {epoch+1}/{num_epochs} ===")
    training(model, train_loader, device, optimizer, scheduler)

    labels, predictions = validation(model, val_loader, device)
    acc = get_accuracy(labels, predictions, 0.5)
    loss = get_mse(labels, predictions)

    print(f"[Validation] Loss: {loss:.4f} | Accuracy: {acc:.2f}%")
