"""
Methods to evaluate quantities from electric field data.

File: libefield.py
"""

class LibEfield():

    quantities = ["er", "ephi", "ez"]

    def evaluate(self, R, phi, z, t, quantity):

        out = None
        if quantity in ["er", "ephi", "ez"]:
            out = self.eval_efield(R, phi, z, t)[quantity]

        assert out is not None, "Unknown quantity"

        return out
