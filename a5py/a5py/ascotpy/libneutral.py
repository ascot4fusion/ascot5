"""
Methods to evaluate quantities from neutral data.

File: libneutral.py
"""

class LibNeutral():

    quantities = ["n0"]

    def evaluate(self, R, phi, z, t, quantity):

        out = None
        if quantity in ["no"]:
            out = self.eval_neutral(R, phi, z, t)[quantity]

        assert out is not None, "Unknown quantity"

        return out
