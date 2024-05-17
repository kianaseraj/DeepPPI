from mtrics import *
from model import FC
from data import train_loader, val_loader
import torch.nn as nn
import torch.nn.functional as F
from torch.optim.lr_scheduler import ReduceLROnPlateau, MultiplicativeLR
import torch_optimizer as optim
import sklearn
from tqdm import tqdm
import mathplotlib.pyplot as plt
import numpy as np
import torch
import math
import pathlib
import os
import glob


device = torch.device("cuda")

print("Data size", len(train_loader), len(val_loader))
     
def training(model, train_loader, val_loader, device, optimizer):
    print(f"training on{len(train_loader)} samples")
    model.train()
    loss_func = nn.MSELoss()
    predictions_tr = torch.Tensor()
    scheduler = MultiStepLR(optimizer, milestones=[1,5], gamma=0.5)
    labels_tr = torch.Tensor()
    for count,(prot1, prot2, label) in enumerate(trainloader):
      prot1 = prot1.to(device)
      prot2 = prot2.to(device)
      optimizer.zero_grad()
      output = model(prot1, prot2)
      predictions_tr = torch.cat((predictions_tr, output.cpu()), 0)
      labels_tr = torch.cat((labels_tr, label.view(-1,1).cpu()), 0)
      loss = loss_func(output, label.view(-1,1).float().to(device))
      loss.backward()
      optimizer.step()
    scheduler.step()
    labels_tr = labels_tr.detach().numpy()
    predictions_tr = predictions_tr.detach().numpy()
    acc_tr = get_accuracy(labels_tr, predictions_tr , 0.5)
    print(f'train_loss : {loss} - train_accuracy : {acc_tr}')
    

def validation(model, val_loader, device):
    model.eval()
    predictions = torch.Tensor()
    labels = torch.Tensor()

    with torch.no_grad():
        for prot1, prot2, label  in val_loader:
            prot1 = prot1.to(device)
            prot2 = prot2.to(device)
            output = model(prot1, prot2)
            predictions = torch.cat((predictions, output.cpu()), 0)
            labels = torch.cat((labels, label.view(-1,1).cpu()), 0)
    labels = labels.detach().numpy()
    predictions = predictions.detach().numpy()
    return labels, predictions



#training

model = FC()
model.to(device)
optimizer =torch.optim.Adam(model.parameters(), lr=0.001)
num_epochs = 30
loss_func = nn.MSELoss()
min_loss = 100
best_accuracy = 0

for epoch in range(num_epochs):
    training(model, train_loader, val_loader, device, optimizer)
    labels, predictions = validation(model, val_loader, device)
    acc = get_accuracy(labels, predictions, 0.5)
    loss = get_mse(labels, predictions)
    print(f'Epoch {epoch}/ {num_epochs} [==============================] - val_loss : {loss} - val_accuracy : {accuracy}')


