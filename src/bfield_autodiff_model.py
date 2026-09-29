import torch
import torch.nn as nn
from typing import List

class BFieldAutodiff(nn.Module):
    def __init__(self):
        super().__init__()

        self.net = nn.Sequential(
            nn.Linear(3, 256),
            nn.SiLU(),
            nn.Linear(256, 256),
            nn.SiLU(),
            nn.Linear(256, 128),
            nn.SiLU(),
            nn.Linear(128, 3),
        )

    def forward(self, x):
        x = x.requires_grad_(True)

        B = self.net(x)

        grads = []

        grads = torch.jit.annotate(List[torch.Tensor], [])

        for k in range(3):
            g_list = torch.autograd.grad(
                [B[:, k].sum()],
                [x],
                create_graph=True,
                retain_graph=True
            )
            
            g = g_list[0]
            
            if g is None:
                g = torch.zeros_like(x)

            grads.append(g)
            
        dB = torch.cat(grads, dim=1)

        out = torch.cat([
            B[:, 0:1], dB[:, 0:1], dB[:, 1:2], dB[:, 2:3],
            B[:, 1:2], dB[:, 3:4], dB[:, 4:5], dB[:, 5:6],
            B[:, 2:3], dB[:, 6:7], dB[:, 7:8], dB[:, 8:9],
        ], dim=1)

        return out

if __name__ == "__main__":
    model = BFieldAutodiff().cuda()
    model.eval()

    x = torch.randn(1024, 3, device="cuda")

    y = model(x)
    print(y.shape)

    scripted = torch.jit.script(model)
    scripted.save("bfield_autodiff.pt")

    y2 = scripted(x)
    print(y2.shape)    
