"""Smoke test: the control-plane dependencies import in the project venv."""


def test_imports():
    import numpy
    import pybind11
    import scipy
    import sgp4
    import skyfield

    for mod in (numpy, scipy, skyfield, sgp4, pybind11):
        assert mod.__name__
