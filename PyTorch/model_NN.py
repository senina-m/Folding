import torch
import torch.nn as nn

batch_size = 4
input_d = 5
output_d = 5
hidden_d = 10

model = torch.nn.Sequential(
    nn.ReLU(input_d),
    nn.ReLU(),
    nn.Linear(hidden_d, output_d),
)
