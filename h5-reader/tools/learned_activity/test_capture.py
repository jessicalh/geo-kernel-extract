import numpy as np
import pytest
import torch
from e3nn import o3
from e3nn.nn import Gate
from torch import nn

from capture import capture_frame


@pytest.fixture(autouse=True)
def double_precision_constants():
    previous = torch.get_default_dtype()
    torch.set_default_dtype(torch.float64)
    yield
    torch.set_default_dtype(previous)


class TinyEncoder(nn.Module):
    def __init__(self):
        super().__init__()
        self.gate = Gate("1x0e", [torch.tanh], "1x0e", [torch.sigmoid], "1x2e")

    def forward(self, value):
        return self.gate(value)


def full_matrix(six):
    xx, xy, xz, yy, yz, zz = six
    return np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])


def test_capture_preserves_prediction_and_rotates_back():
    torch.manual_seed(41)
    model = TinyEncoder().double().eval()
    inputs = torch.randn(3, model.gate.irreps_in.dim, dtype=torch.float64)
    positions = np.zeros((3, 3))
    expected = model(inputs)
    prediction, original = capture_frame(
        model,
        lambda: model(inputs),
        frame_index=0,
        raw_positions=positions,
        atom_indices=[0, 2],
        extraction_to_model_rotation=np.eye(3),
    )
    torch.testing.assert_close(prediction, expected, rtol=0, atol=0)
    rotation = o3.rand_matrix(dtype=torch.float64)
    representation = model.gate.irreps_in.D_from_matrix(rotation)
    _, rotated = capture_frame(
        model,
        lambda: model(inputs @ representation.T),
        frame_index=0,
        raw_positions=positions,
        atom_indices=[0, 2],
        extraction_to_model_rotation=rotation.numpy(),
    )
    np.testing.assert_allclose(
        original["channels"][0]["tensors"], rotated["channels"][0]["tensors"], atol=1e-7
    )
    tensor = full_matrix(original["channels"][0]["tensors"][0])
    assert abs(np.trace(tensor)) < 1e-12
    assert np.linalg.norm(tensor) == pytest.approx(float(expected[0, 1:].norm()))
    assert not model.gate._forward_hooks


def test_failed_forward_removes_hooks():
    model = TinyEncoder().eval()

    def fail():
        raise RuntimeError("deliberate test failure")

    with pytest.raises(RuntimeError, match="deliberate"):
        capture_frame(
            model,
            fail,
            frame_index=0,
            raw_positions=np.zeros((1, 3)),
            atom_indices=[0],
            extraction_to_model_rotation=np.eye(3),
        )
    assert not model.gate._forward_hooks
