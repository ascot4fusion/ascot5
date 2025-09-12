import pytest
import numpy as np
from unyt import m, s, T, kg, J, A

import a5py.physlib.formulas as formulas
from a5py.physlib.formulas import Quantity, resolve_quantities

@pytest.fixture
def isolated_registry():
    old_registry = Quantity.registry.copy()
    Quantity.registry.clear()
    yield
    Quantity.registry.clear()
    Quantity.registry.update(old_registry)


@pytest.fixture
def momentum():
    momentum = Quantity("momentum", "kg*m/s")
    momentum.add_formula(
        lambda m, v: m * v, ["mass", "velocity"]
        )
    momentum.add_formula(
        lambda pr, pphi, pz: pr + pphi + pz, ["pr", "pphi", "pz"]
        )
    return momentum


def test_quantity_variables(isolated_registry, momentum):
    assert momentum.variables == [("mass", "velocity"), ("pr", "pphi", "pz"),]


@pytest.mark.parametrize(
    "variables, expected",
    [
        ({}, None),
        ({"mass": 1*kg}, None),
        ({"velocity": 1*m/s, "pr": 1*kg*m/s}, None),
        ({"mass": 1*kg, "velocity": 2*m/s}, 2*kg*m/s),
        ({"pr": 3*kg*m/s, "pphi": 4*kg*m/s, "pz": 5*kg*m/s}, 12*kg*m/s),
        ({"mass": 2*kg, "pr": 3*kg*m/s, "pphi": 4*kg*m/s, "pz": 5*kg*m/s}, 12*kg*m/s),
    ],
    ids=["no-vars", "one-missing", "mixed-but-missing", "formula1", "formula2",
         "extra-vars",],
)
def test_quantity_compute(isolated_registry, momentum, variables, expected):
    if expected is None:
        with pytest.raises(RuntimeError):
            momentum.compute(**variables)
    else:
        assert momentum.compute(**variables) == expected


@pytest.mark.parametrize(
    "quantity, variables, expected",
    [
        ("pnorm", {"pr": 1*kg*m/s, "pphi": 2*kg*m/s, "pz": 3*kg*m/s}, 3.7416573867739413),
        ("bnorm", {"br": 1*T, "bphi": 2*T, "bz": 3*T}, 3.7416573867739413),
        ("divb", {"br": 1*T, "brdr": 2*T/m, "bphidphi": 3*T, "bzdz": 4*T/m, "r": 5*m}, 6.8),
        ("gradbr", {"br": 1*T, "brdr": 2*T/m, "bphi": 3*T, "bphidr": 4*T/m,
                    "bz": 5*T, "bzdr": 6*T/m,}, 7.437357441610946),
        ("gradbphi", {"br": 1*T, "brdphi": 2*T, "bphi": 3*T, "bphidphi": 4*T,
                      "bz": 5*T, "bzdphi": 6*T, "r": 7*m}, 1.0624796345158494),
        ("gradbz", {"br": 1*T, "brdz": 2*T/m, "bphi": 3*T, "bphidz": 4*T/m,
                    "bz": 5*T, "bzdz": 6*T/m,}, 7.437357441610946),
        ("curlbr", {"bzdphi": 1*T, "bphidz": 2*T/m, "r": 3*m}, -1.66666667),
        ("curlbphi", {"brdz": 1*T/m, "bzdr": 2*T/m}, -1.0),
        ("curlbz", {"bphi": 1*T, "bphidr": 2*T/m, "brdphi": 3*T, "r": 4*m}, 1.5),
        ("jr", {"curlbr": 1*T/m}, 0.79577471545948e6),
        ("jphi", {"curlbphi": 1*T/m}, 0.79577471545948e6),
        ("jz", {"curlbz": 1*T/m}, 0.79577471545948e6),
        ("jnorm", {"jr": 1*A/m**2, "jphi": 2*A/m**2, "jz": 3*A/m**2}, 3.7416573867739413),
    ]
)
def test_values(quantity, variables, expected):
    val = getattr(formulas, quantity).compute(**variables)
    assert np.isclose(val, expected)


@pytest.mark.parametrize(
    "available, requested, needed, computethese",
    [
        (("cake",), ("cake",), {"cake",}, {"cake",}),
        (("cake", "sugar",), ("cake",), {"cake",}, {"cake",}),
        (("cake", "dough", "icing",), ("cake",), {"cake"}, {"cake"}),
        (("batter",), ("cake",), {"batter",}, {"cake", "batter",}),
        (("dough", "icing",), ("cake",), {"dough", "icing"}, {"cake", "dough", "icing"}),
        (("dough", "sugar",), ("cake",), {"dough", "sugar"}, {"cake", "dough", "icing", "sugar"}),
        (("cake", "icing",), ("dough",), {"cake", "icing"}, {"cake", "dough", "icing"}),
        (("cake",), ("dough",), None, None),
    ],
    ids=(
        "cake-can-be-made-from-cake",
        "cake-can-be-made-from-cake-no-sugar-needed",
        "cake-can-be-made-from-cake-even-if-ingredients-available",
        "cake-can-be-made-from-dough-and-icing",
        "cake-can-be-made-from-batter",
        "cake-can-be-made-from-dough-and-icing-made-with-sugar",
        "dough-can-be-made-with-cake-and-icing",
        "cake-cant-be-made-with-just-dough",
        ),
    )
def test_resolve_quantities(
    isolated_registry, available, requested, needed, computethese,
    ):
    """Test resolve_quantities.

    - Cake is made from dough and icing, or from batter.
    - Dough is made from cake and icing.
    - Icing is made from sugar.
    """
    q1 = Quantity("cake", "1")
    q1.add_formula(
        lambda b, c: b + c, ["dough", "icing"]
        )
    q1.add_formula(
        lambda e: e, ["batter"]
        )
    q2 = Quantity("dough", "1")
    q2.add_formula(
        lambda a, c: a + c, ["cake", "icing"]
        )
    q3 = Quantity("icing", "1")
    q3.add_formula(
        lambda d: d, ["sugar"]
        )

    if computethese is not None:
        needed_, computethese_ = resolve_quantities(available, requested)
        assert computethese_ == computethese
        assert needed_ == needed
    else:
        with pytest.raises(RuntimeError):
            resolve_quantities(requested, available)
