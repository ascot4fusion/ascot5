import torch
import torch.nn as nn
from typing import List


class BFieldAutodiffWrapper(nn.Module):
    def __init__(self, base_model):
        super().__init__()
        self.base = base_model

    def forward(self, x):
        x = x.requires_grad_(True)

        B = self.base(x)

        grads = torch.jit.annotate(List[torch.Tensor], [])

        for k in range(3):
            g_list = torch.autograd.grad(
                [B[:, k].sum()],
                [x],
                create_graph=False,
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
    base = torch.jit.load("bfield_B_trained.pt")
    base.eval()

    wrapped = BFieldAutodiffWrapper(base).cuda()
    wrapped.eval()

    x = torch.randn(1024, 3, device="cuda")

    y = wrapped(x)
    print(y.shape)

    scripted = torch.jit.script(wrapped)
    scripted.save("bfield_B_autodiff.pt")

    y2 = scripted(x)
    print(y2.shape)
