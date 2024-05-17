import torch
from torch.utils.data import Dataset, DataLoader
import numpy as np
import os
import glob
from os.path import isfile, join

feature_dir = "/fs/pool/pool-schwille-user/Frohn_Bela/Seraj-Kiana/ProteinEvolution/ProteinEvolver/ppi/deep_ppi/seq_feature"
train_data =  np.load("/fs/pool/pool-schwille-user/Frohn_Bela/Seraj-Kiana/ProteinEvolution/ProteinEvolver/ppi/deep_ppi/ppi_data/training_data_red_1500_0_duplicate.npy")
val_data = np.load("/fs/pool/pool-schwille-user/Frohn_Bela/Seraj-Kiana/ProteinEvolution/ProteinEvolver/ppi/deep_ppi/ppi_data/validation_data_red_1500_0.npy")


class ProtDataset(Dataset):
    def __init__(self, feature_dir : str, data_file : str):
        self.feature_dir = feature_dir
        self.prot1 = data_file[:,0]
        self.prot2 = data_file[:,1]
        self.label = data_file[:,2]
        self.n_samples = data_file.shape[0]


    def __len__(self):
        return self.n_samples

    def __getitem__(self, key):
        prot1 = os.path.join(self.feature_dir, self.prot1[key]+".pt")
        prot2 = os.path.join(self.feature_dir, self.prot2[key]+".pt")
       
        prot1  = torch.load(glob.glob(prot1)[0])
        prot2  = torch.load(glob.glob(prot2)[0])
       
        label = torch.tensor(self.label[key])
      
        return prot1, prot2, prot1_val, prot2_val, label


train_dataset = ProtDataset(feature_dir, train_data)
val_dataset = ProtDataset(feature_dir, val_data)


train_loader = DataLoader(train_dataset, batch_size=16, num_workers=0)
val_loader = DataLoader(val_dataset, batch_size=16, num_workers=0)

print("Data size", len(train_loader), len(val_loader))
        
        


