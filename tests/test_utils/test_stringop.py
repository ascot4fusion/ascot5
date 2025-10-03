"""Tests for a5py.utils.stringop."""
from datetime import datetime

from a5py import utils


def test_format2universaldate_correct_format():
    """The function should format datetime into ASCOT5's universal format.

    This is "YYYY-MM-DD HH:MM:SS" with zero padding where appropriate.
    """
    dt = datetime(2023, 5, 17, 8, 3, 45)
    result = utils.format2universaldate(dt)
    assert result == "2023-05-17 08:03:45"


def test_decorate():
    """Test decoration."""
    s = "hello"
    result = utils.decorate(s)
    assert utils.undecorate(result) == s
    assert result.endswith("\033[0m")


def test_undecorate():
    """Test undecoration."""
    decorated = "\033[01m\033[32mHello\033[0mWorld\033[04m!\033[0m"
    assert utils.undecorate(decorated) == "HelloWorld!"

def test_combination_of_effects():
    """Test that all decorations are in place."""
    s = "combo"
    result = utils.decorate(s, color="green", bold=True, underline=True)
    assert "\033[32m" in result
    assert "\033[01m" in result
    assert "\033[04m" in result
    assert result.endswith("\033[0m")
    assert utils.undecorate(result) == s
