import argparse
import math
from typing import List

import torch
import torch.nn as nn
import torch.nn.functional as F

class BFieldAutodiffModel(nn.Module):
    """Periodic SiLU B-field model with derivative/divergence training outputs.

    Input is [normalized R, normalized phi, normalized Z]. Output contains
    3 normalized B values, 9 normalized spatial gradients, and 1 weighted
    normalized cylindrical-divergence residual.
    """
    def __init__(self):
        super().__init__()
        self.fc0 = nn.Linear(4, 256)
        self.fc1 = nn.Linear(256, 256)
        self.fc2 = nn.Linear(256, 128)
        self.fc3 = nn.Linear(128, 3)
        self.b0 = nn.Parameter(torch.zeros(256))
        self.b1 = nn.Parameter(torch.zeros(256))
        self.b2 = nn.Parameter(torch.zeros(128))

    def field(self, x):
        phase = 18.0 * math.pi * (x[:, 1] + 1.0)
        features = torch.stack(
            (x[:, 0], torch.sin(phase), torch.cos(phase), x[:, 2]), dim=1
        )
        h = F.silu(self.fc0(features) + self.b0)
        h = F.silu(self.fc1(h) + self.b1)
        h = F.silu(self.fc2(h) + self.b2)
        return self.fc3(h)

    def forward(self, x):
        x = x.requires_grad_(True)
        b = self.field(x)
        gradients = torch.jit.annotate(List[torch.Tensor], [])
        for k in range(3):
            g = torch.autograd.grad(
                [b[:, k].sum()], [x], create_graph=True, retain_graph=True
            )[0]
            if g is None:
                g = torch.zeros_like(x)
            gradients.append(g)

        # Scales are literal because TorchScript cannot capture module globals.
        grad_scale = torch.tensor(
            [[3.55015, 0.52316, 4.50280],
             [0.53320, 0.06598, 0.00752],
             [2.94440, 0.12928, 2.03438]],
            dtype=x.dtype, device=x.device)
        out = torch.jit.annotate(List[torch.Tensor], [])
        for k in range(3):
            out.append(b[:, k:k+1])
            for d in range(3):
                out.append(gradients[k][:, d:d+1] / grad_scale[k, d])

        b_scale = torch.tensor(
            [1.953079, 17.53601, 3.749450], dtype=x.dtype, device=x.device)
        dbr_dR = gradients[0][:, 0] * (b_scale[0] / 2.0)
        bR_over_R = b[:, 0] * b_scale[0] / (x[:, 0] * 2.0 + 6.0)
        dbphi_dphi_over_R = (
            gradients[1][:, 1] * (b_scale[1] / math.pi)
            / (x[:, 0] * 2.0 + 6.0)
        )
        dbz_dZ = gradients[2][:, 2] * (b_scale[2] / 2.0)
        divergence = dbr_dR + bR_over_R + dbphi_dphi_over_R + dbz_dZ
        # sqrt(3 * lambda_div) / div_scale for lambda_div=0.1.
        div_output = divergence * (0.5477225575 / 3.03743887)
        out.append(div_output.unsqueeze(1))
        return torch.cat(out, dim=1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--initial-weights", default=None,
                        help="optional compatible SiLU-periodic state_dict")
    parser.add_argument("--output", default="bfield_B_periodic_dB_init.pt")
    args = parser.parse_args()
    model = BFieldAutodiffModel()
    if args.initial_weights:
        state = torch.load(args.initial_weights, map_location="cpu", weights_only=True)
        model.load_state_dict(state)
    scripted = torch.jit.script(model.eval())
    scripted.save(args.output)
    print("Saved trainable TorchScript model:", args.output)


if __name__ == "__main__":
    main()
