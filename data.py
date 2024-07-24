import torch
from torch.utils.data import Dataset, DataLoader
import numpy as np
import os


train_data =  np.load("/fs/pool/pool-schwille-user/Frohn_Bela/Seraj-Kiana/ProteinEvolution/ProteinEvolver/ppi/deep_ppi/ppi_data/training_data_red_1500_0_duplicate.npy")
val_data = np.load("/fs/pool/pool-schwille-user/Frohn_Bela/Seraj-Kiana/ProteinEvolution/ProteinEvolver/ppi/deep_ppi/ppi_data/validation_data_red_1500_0.npy")
feature_dir = "/fs/pool/pool-schwille-user/Frohn_Bela/Seraj-Kiana/ProteinEvolution/ProteinEvolver/ppi/deep_ppi/tensor_feature"

class ProtDataset(Dataset):
    def __init__(self, feature_dir : str, data_file : np.array):
        self.feature_dir = feature_dir
        self.prot_1 = data_file[:,0]
        self.prot_2 = data_file[:,1]
        self.label_ = data_file[:,2].astype(float).astype(int)
        self.n_samples = data_file.shape[0]


    def __getitem__(self, index):
        """
        Gets one sample from the dataset and formats the features to be useable by
        the NN, cuda tensor.
        """

        prot1 = os.path.join(self.feature_dir, self.prot_1[index]+".pt")
        prot2 = os.path.join(self.feature_dir, self.prot_2[index]+".pt")

        prot1  = torch.load(prot1)
        prot2  = torch.load(prot2)

        label = torch.tensor(self.label_[index])

        return prot1, prot2, label

    def __len__(self):
        """
        returns the size of dataset
        """
        return self.n_samples


train_dataset = ProtDataset(feature_dir, train_data)
val_dataset = ProtDataset(feature_dir, val_data)


train_loader = DataLoader(train_dataset, batch_size=256, num_workers=0)
val_loader = DataLoader(val_dataset, batch_size=256, num_workers=0)


print("train_data:", train_data.shape,"\n","val_data:", val_data.shape)



