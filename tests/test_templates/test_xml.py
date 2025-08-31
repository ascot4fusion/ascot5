import pytest

from a5py.templates import make_simple_type, make_element_block, doc, make_schema


def test_test():
    #print(make_simple_type("IntegerPositive", "xs:integer", min_val=1))
    #print(make_simple_type("IntegerBinary", "xs:integer", min_val=0, max_val=1))
    #print(make_simple_type("FloatPositive", "xs:float", min_val=0.0, min_inclusive=False))

    params = {
    "SIM_MODE": "Integer1234",
    "DT": "FloatPositive",
    }

    #print(make_element_block("Simulation",
    #                        "Simulation mode and time-step",
    #                        params))
    #print(doc("simulation_mode", "FloatPositive"))
    make_schema()
    assert 1 == 2
