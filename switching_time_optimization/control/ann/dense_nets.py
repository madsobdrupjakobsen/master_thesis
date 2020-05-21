import torch
import torch.nn as nn
import torch.nn.functional as F


class simpleNet(nn.Module):

    def __init__(self,n_s_all=6):
        super(simpleNet, self).__init__()

        # an affine operation: y = Wx + b
        self.n_s_all = n_s_all
        self.fc1 = nn.Linear(49, 120)  # 6*6 from image dimension
        self.fc2 = nn.Linear(120, 84)
        self.fc3 = nn.Linear(84, self.n_s_all)

    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = self.fc3(x)
        return x

    # not using in danse net https://pytorch.org/tutorials/beginner/blitz/neural_networks_tutorial.html
    def num_flat_features(self, x):
        size = x.size()[1:]  # all dimensions except the batch dimension
        num_features = 1
        for s in size:
            num_features *= s
        return num_features

class simpleNet_RK(nn.Module):

    def __init__(self,n_s_all=6):
        super(simpleNet_RK, self).__init__()

        # an affine operation: y = Wx + b
        self.n_s_all = n_s_all
        self.fc1 = nn.Linear(49+48, 120)  # 6*6 from image dimension
        self.fc2 = nn.Linear(120, 84)
        self.fc3 = nn.Linear(84, self.n_s_all)

    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = self.fc3(x)
        return x

    # not using in danse net https://pytorch.org/tutorials/beginner/blitz/neural_networks_tutorial.html
    def num_flat_features(self, x):
        size = x.size()[1:]  # all dimensions except the batch dimension
        num_features = 1
        for s in size:
            num_features *= s
        return num_features




class simpleClassNet(nn.Module):

    def __init__(self):
        super(simpleClassNet, self).__init__()

        # an affine operation: y = Wx + b
        self.fc1 = nn.Linear(2, 120) 
        self.fc2 = nn.Linear(120, 84)
        self.fc3 = nn.Linear(84, 3)

    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = self.fc3(x)
        return x

class convClassNet(nn.Module):

    def __init__(self):
        super(simpleClassNet, self).__init__()

        # an affine operation: y = Wx + b
        self.fc1 = nn.Linear(2, 120) 
        self.fc2 = nn.Linear(120, 84)
        self.fc3 = nn.Linear(84, 3)

    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = self.fc3(x)
        return x
