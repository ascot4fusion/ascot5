import torch
import torch.nn as nn

from bfield_autodiff_model import BFieldAutodiffModel


class BFieldAutodiffInference(nn.Module):
    """Restore physical B and dB from the trained TorchFort model outputs."""
    def __init__(self, trained_model):
        super().__init__()
        self.base = trained_model

    def forward(self, x):
        y = self.base(x)
        b_scale = torch.tensor([1.953079, 17.53601, 3.749450], dtype=x.dtype, device=x.device)
        grad_scale = torch.tensor(
            [[3.55015, 0.52316, 4.50280],
             [0.53320, 0.06598, 0.00752],
             [2.94440, 0.12928, 2.03438]],
            dtype=x.dtype, device=x.device)
        deriv_coord_scale = torch.tensor(
            (0.5, 1.0 / 3.141592653589793, 0.5),
            dtype=x.dtype, device=x.device)
        out = torch.jit.annotate(list[torch.Tensor], [])
        output_base = (0, 4, 8)
        for k in range(3):
            out.append(y[:, output_base[k]:output_base[k] + 1] * b_scale[k])
            grad_base = output_base[k] + 1
            for d in range(3):
                normalized_grad = y[:, grad_base + d:grad_base + d + 1] * grad_scale[k, d]
                physical_grad = normalized_grad * b_scale[k] * deriv_coord_scale[d]
                out.append(physical_grad)
        return torch.cat(out, dim=1)


def load_bfield_model(checkpoint_path):
    """Reconstruct the periodic model from TorchFort's saved parameters."""
    checkpoint = torch.jit.load(checkpoint_path, map_location="cpu")
    params = dict(checkpoint.named_parameters())
    model = BFieldAutodiffModel()
    with torch.no_grad():
        for name, parameter in model.named_parameters():
            parameter.copy_(params[name])
    return model.eval()


if __name__ == "__main__":
    base = load_bfield_model("bfield_B_periodic_dB_trained.pt").cuda()
    wrapped = torch.jit.script(BFieldAutodiffInference(base).eval())
    x = torch.randn(1024, 3, device="cuda")
    y = wrapped(x)
    assert y.shape == (1024, 12)
    print("Inference output shape:", tuple(y.shape))
    wrapped.save("bfield_B_autodiff.pt")
