"""Test options."""
import os

import pytest
import unyt

from a5py.data import AscotData
from a5py.data.options import (
    constraints, OptionsStruct, SimulationOptions, HistParams, Dimension,
    )
from a5py.data.options.parameters import Simulation, ParameterClass
from a5py.data.access.hdf5io import TreeFileManager

@pytest.fixture()
def file():
    datamanager = TreeFileManager(
        root="root",
        filename="testascot.h5",
        file_exists=False,
        input_categories=[],
        )
    datamanager.set_note(name="options", node="root", note="") # Creates node
    yield datamanager.access_data(name="options", node="root")

    os.unlink("testascot.h5")

@pytest.fixture
def path():
    yield "test.toml"
    try:
        os.unlink("test.toml")
    except FileNotFoundError:
        pass


def test_write_hdf5(file):
    opt = SimulationOptions.from_dict(
        simulation={"mode": "field-line"},
        histograms=[
            {
            "dimensions": [("r", 0, 1, 10)],
            "charge_interval": (1,1),
            }
        ],
        )
    opt._write_hdf5(file)

def test_read_hdf5(file):
    opt = SimulationOptions.from_dict(
        simulation={"mode": "field-line"},
        histograms=[
            {
            "dimensions": [("r", 0, 1, 10)],
            "charge_interval": (1,1),
            }
        ],
        )
    opt._write_hdf5(file)

    opt2 = SimulationOptions.from_hdf5(file)
    assert opt2.simulation.mode == opt.simulation.mode
    assert opt2.histograms[0].dimensions[0] == opt.histograms[0].dimensions[0]

def test_options_attribute_doc():
    doc = Simulation()
    doc = doc._attribute_doc("mode")
    print(doc)


@pytest.mark.parametrize(
    "attr",
    ["single_line", "two_line", "multi_line"],
    ids=["single_line", "two_line", "multi_line"],
)
def test_options_constraints_parse(attr):

    class TestParameters(ParameterClass):
        """For testing parameters.
        
        Attributes
        ----------
        single_line : float
            Single line: (> 0), default=1.
        two_line : float
            Two-line:
            (> 0), default=1.
        multi_line : float
            Multi-line split at
            constraint:
            (> 0),
            default=1.
        """

    obj = TestParameters()
    doc = obj._attribute_doc(attr)
    info = constraints.parse(doc)
    assert info["constraint"] == "(> 0)"
    assert info["default"] == "1"


@pytest.mark.parametrize(
    "constraint, valid_values, invalid_values",
    [
        ("{1, 2, 3, 4}", (1, 2, 3, 4), (0, 5, 1.5)),
        ("(x)", (1,), ()),
        ("(> 0)", (1,), (0,)),
        ("[a, b]", ([1, 2], (1, 2)), ([], [1, 2, 3])),
        ("[a > 0, b > 0]", ([1, 1],), ([0, 1], [1, 0], [0, 0])),
        ("[0 <= a < 360, 0 < b <= 360]", ([1, 2],), ([-1, 2], [1, 0], [360, 1], [1, 361])),
        ("[0 <= a < 360, ...]", (0, [0], [1, 359],), ([0, 361],)),
    ],
    ids=[
        "sequence", "no-constraint", "positive", "2-array-no-limits",
        "2-array-lower-limit", "2-array-both-limits", "unlimited-array-with-limits",
    ]
    )
def test_options_constraints_enforce(constraint, valid_values, invalid_values):
    for value in valid_values:
        assert constraints.enforce(value, constraint)

    for value in invalid_values:
        assert not constraints.enforce(value, constraint)



@pytest.mark.parametrize(
        "params",
        [
            {"mode": "field-line"},
        ],
    ids=["valid-parameters"],
)
def test_options_from_dict(params):
    opt = SimulationOptions.from_dict(
        simulation=params,
    )
    assert opt.simulation.mode == "field-line"

@pytest.mark.parametrize(
        "params",
        [
            {"unknown": {"mode": "field-line"}},
            {"simulation": {"unknown": "field-line"}},
            {"simulation": {"mode": "incorrect"}},
        ],
    ids=["unknown-group", "unknown-parameter", "invalid-value"],
)
def test_options_from_dict_raises(params):
    with pytest.raises(ValueError):
        SimulationOptions.from_dict(**params)


@pytest.mark.parametrize(
        "params, msg",
        [
            ({"simulation": {"unknown": "field-line"}}, ""),
            ({"simulation": {"mode": 1}}, ""),
        ],
    ids=["unknown-group", "invalid-value"],
)
def test_options_set_attr_raises(params, msg):
    opt = SimulationOptions()
    with pytest.raises(ValueError):
        for group, settings in params.items():
            for setting, value in settings.items():
                setattr(getattr(opt, group), setting, value)

def test_options_set_group_raises():
    opt = SimulationOptions()
    with pytest.raises(TypeError):
        opt.simulation = 1

    with pytest.raises(ValueError):
        opt.unknown = 1

def test_options_histograms():
    opt = SimulationOptions()
    assert opt.histograms == ()

    opt = SimulationOptions.from_dict(
        histograms=[
            {
            "dimensions": [("r", 0, 1, 10)],
            "charge_interval": (1,1),
            }
        ],
    )
    assert opt.histograms == (
        HistParams(
            dimensions=(Dimension("r", 0, 1, 10),),
            charge_interval=(1,1),
        ),
    )
    opt.histograms[0].charge_interval = (2,2)
    assert opt.histograms[0].charge_interval == (2,2)


@pytest.mark.parametrize(
    "params, err",
    [
        ({"dimensions": [("unknown", 0., 1., 10)]}, ValueError),
        ({"dimensions": [("charge", 0., 1., 10)]}, ValueError),
        ({"dimensions": [("r", 0., 1., -10)]}, ValueError),
        ({"dimensions": [("r", 2., 1., 10)]}, ValueError),
        ({"dimensions": [("r", 0., 1., 1.1)]}, ValueError),
        ({"dimensions": [("r", 0., 1.*unyt.s, 10)]}, ValueError),
        ({"dimensions": [("r", 0., 1., 10), ("r", 0., 1., 10),]}, ValueError),
        ({"dimensions": []}, ValueError),
        ({"dimensions": [(),]}, TypeError),
        ({"dimensions": [("r", 0., 1., 10),], "charge_interval": (1,0)}, ValueError),
    ]
)
def test_options_histograms_raises(params, err):
    with pytest.raises(err):
        SimulationOptions.from_dict(histograms=[params])


@pytest.mark.parametrize(
    ("params", "field", "expected"),
    [
        ({"simulation": {"mode": "field-line"}}, "simulation_mode", 4),
    ]
)
def test_options_options_struct_from_params(params, field, expected):
    opt = SimulationOptions.from_dict(**params)
    struct = OptionsStruct.from_params(opt)
    assert getattr(struct, field) == expected


def test_options_write_toml(path):
    opt = SimulationOptions.from_dict(
        simulation={"mode": "field-line"},
        histograms=[{
            "dimensions": [("r", 0, 1, 10)],
        },
        {
            "dimensions": [("r", 0, 2, 10)],
        },
        ]
    )
    opt.write_toml(path, include_comments="detailed")
    #assert os.path.exists(path)


def test_options_read_toml(path):
    opt = SimulationOptions.from_dict(
        simulation={"mode": "field-line"},
        histograms=[{
            "dimensions": [("r", 0, 1, 10)],
        }]
    )
    opt.write_toml(path)
    opt = SimulationOptions.from_toml(path)
    assert opt.simulation.mode == "field-line"


def test_options_freeze():
    opt = SimulationOptions.from_dict(
        simulation={"mode": "field-line"},
    )
    assert opt.simulation.mode == "field-line"
    opt.simulation.mode = "hybrid"
    assert opt.simulation.mode == "hybrid"
    opt._freeze()
    with pytest.raises(AttributeError):
        opt.simulation.mode = "field-line"
