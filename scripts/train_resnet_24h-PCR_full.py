#!/usr/bin/env python
# -*- coding:utf-8 -*-

import numpy as np
import sys

import torch
import torch.nn as nn
import torch.nn.functional as F

from torch.utils.data import DataLoader

from torch.distributions.beta import Beta
from torch.distributions.normal import Normal
from torch.distributions.constraints import positive

ISZ = 87 # With chromatin features.
OSZ = 1
NN = 100

def InverseLinear(x):
   # Inverse-Linear activation function
   return 1.0 / (1.0-x+torch.abs(x)) + x+torch.abs(x)

class MLP(nn.Module):

   def __init__(self):
      super(MLP, self).__init__()
      self.lyrE1 = nn.Sequential(
         nn.Linear(ISZ, NN),
         nn.Dropout(.3),
         nn.BatchNorm1d(NN),
         nn.ReLU()
      )
      self.lyrE2 = nn.Sequential(
         nn.Linear(NN, NN),
         nn.Dropout(.3),
         nn.BatchNorm1d(NN),
         nn.ReLU(),
         nn.Linear(NN, NN),
         #nn.Dropout(.1),
         nn.BatchNorm1d(NN),
      )
      self.a = nn.Linear(NN,OSZ)
      self.b = nn.Linear(NN,OSZ)

   def forward(self, x): 
      x = self.lyrE1(x)
      x = torch.relu(x + self.lyrE2(x))
      a = torch.clamp(InverseLinear(self.a(x)/2), min=.01, max=100)
      b = torch.clamp(InverseLinear(self.b(x)/2), min=.01, max=100)
      return (a,b)

   def optimize(self, train_data, test_data, epochs=100, bsz=256):
      # The initial learning rates are set to avoid the parameters
      # to blow up. If they are higher no learning takes place.
      optimizer = \
            torch.optim.Adam(self.parameters(), lr=0.001)

      batches = DataLoader(dataset=train_data,
            batch_size=bsz, shuffle=True)
      test_set = DataLoader(dataset=test_data,
            batch_size=bsz, shuffle=True)

      best = float('inf')

      for ep in range(epochs):
         batch_loss = 0.0
         self.train()
         for bno,data in enumerate(batches):
            y = torch.clamp(data[:,:1], min=.001, max=.999)
            feat = data[:,1:]
            (a,b) = self(feat)
            loss = -torch.mean(Beta(a,b).log_prob(y))
            batch_loss += float(loss)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

         # Test data.
         self.eval()
         with torch.no_grad():
            test_rcst = 0.0
            for sno,data in enumerate(test_set):
               y = torch.clamp(data[:,:1], min=.001, max=.999)
               feat = data[:,1:]
               (a,b) = self(feat)
               test_rcst -= float(torch.mean(Beta(a,b).log_prob(y)))

         # Print logs on stderr.
         sys.stdout.write('%d\t%f\t%f\n' % \
               (ep, batch_loss / bno, test_rcst / sno))


def make_data(f):
   data = list()
   next(f) # Skip header.
   for line in f:
      items = line.split()
      data.append(torch.tensor([float(x) for x in items[:-4]]))
   data = torch.stack(data)
   # Select time point 24 hours.
   data = data[data[:,5] == 0 ,:]
   # Select method 6xPCR.
   data = data[data[:,6] == 1 ,:]
   # Use conflict mode (True only if 0 or 1).
   data[:,0] = data[:,0] == data[:,0]*data[:,0]
   print(data.shape)
   # Keep 10% for testing.
   rng = torch.rand(data.shape[0])
   return data[rng < .1,:], data[rng >= .1,:]


if __name__ == '__main__':
   with open('barcode_view_all_features_with_mapping.txt') as f:
      test_data, train_data = make_data(f)
   E = MLP().cuda()
   E.optimize(train_data.to('cuda'), test_data.to('cuda'))
