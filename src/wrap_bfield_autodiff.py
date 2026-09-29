import torch
import torch.nn as nn
from typing import List


class BFieldMLP(nn.Module):
    """Recreate the TorchFort MLP whose state is saved in the checkpoint."""
    def __init__(self):
        super().__init__()
        self.fc0 = nn.Linear(3, 256)
        self.fc1 = nn.Linear(256, 256)
        self.fc2 = nn.Linear(256, 128)
        self.fc3 = nn.Linear(128, 3)
        self.b0 = nn.Parameter(torch.zeros(256))
        self.b1 = nn.Parameter(torch.zeros(256))
        self.b2 = nn.Parameter(torch.zeros(128))

    def forward(self, x):
        x = torch.relu(self.fc0(x) + self.b0)
        x = torch.relu(self.fc1(x) + self.b1)
        x = torch.relu(self.fc2(x) + self.b2)
        return self.fc3(x)


def load_bfield_model(checkpoint_path):
    """Restore the TorchFort weight-only checkpoint as a callable MLP."""
    checkpoint = torch.jit.load(checkpoint_path, map_location="cpu")
    params = dict(checkpoint.named_parameters())
    model = BFieldMLP()
    with torch.no_grad():
        for name in ("fc0", "fc1", "fc2", "fc3"):
            getattr(model, name).weight.copy_(params[name + ".weight"])
            getattr(model, name).bias.copy_(params[name + ".bias"])
        for name in ("b0", "b1", "b2"):
            getattr(model, name).copy_(params[name])
    return model.eval()


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
                [B[:, k].sum()], [x], create_graph=False, retain_graph=True
            )
            g = g_list[0]
            if g is None:
                g = torch.zeros_like(x)
            grads.append(g)
        dB = torch.cat(grads, dim=1)
        return torch.cat([
            B[:, 0:1], dB[:, 0:1], dB[:, 1:2], dB[:, 2:3],
            B[:, 1:2], dB[:, 3:4], dB[:, 4:5], dB[:, 5:6],
            B[:, 2:3], dB[:, 6:7], dB[:, 7:8], dB[:, 8:9],
        ], dim=1)


if __name__ == "__main__":
    base = load_bfield_model("bfield_B_trained.pt").cuda()
    wrapped = torch.jit.script(BFieldAutodiffWrapper(base).eval())
    x = torch.randn(1024, 3, device="cuda")
    y = wrapped(x)
    print(y.shape)
    wrapped.save("bfield_B_autodiff.pt")
