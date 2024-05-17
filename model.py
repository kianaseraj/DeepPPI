import torch
import torch.nn as nn
import torch.nn.functional as F 


class FC(nn.Module):

    def __init__(self):
        super(FC, self).__init__()

        #protein 1
        self.pro1_fc1 = nn.Linear(1164, 512)
        self.pro1_fcd1 = nn.Dropout(0.2)
        self.pro1_fc2 = nn.Linear(512, 256)
        self.pro1_fc2d = nn.Dropout()
        self.pro1_fc3 = nn.Linear(256, 128)
        self.pro1_fc3d = nn.Dropout()


        #protein 2  
        self.pro2_fc1 = nn.Linear(1164, 512)
        self.pro2_fcd1 = nn.Dropout(0.2)
        self.pro2_fc2 = nn.Linear(512, 256)
        self.pro2_fc2d = nn.Dropout()
        self.pro2_fc3 = nn.Linear(256, 128)
        self.pro2_fc3d = nn.Dropout()


        #combined_layer 
        self.fc1 = nn.Linear(256, 128)
        self.fcd1 = nn.Dropout(0.2)
        self.out = nn.Linear(128, 2)



    def forward(self, pro1_data, pro2_data):

        
