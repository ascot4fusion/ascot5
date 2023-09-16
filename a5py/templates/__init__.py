"""Package for generating inputs from templates and imported data.
"""
from .analyticalinputs import AnalyticalInputs
from .optiontemplates  import OptionTemplates
from .poincare  import PoincareTemplates
from .importdata import ImportAdas
from. convertascot4 import Ascot4Templates

class InputFactory(AnalyticalInputs, OptionTemplates, PoincareTemplates,
                   ImportAdas, Ascot4Templates):
    """Class for creating input data from templates or imported data.

    The templates are constructed by calling :meth:`construct` and specifying
    the template. The templates are constructed in other (public) methods
    of this class. The name of the template is the name of the method, and
    details of the template can be found from the method description.
    """

    def __init__(self, ascot):
        """Store class:`Ascot` object needed for input creation.
        """
        self._ascot = ascot

    def construct(self, template, **kwargs):
        """Create input from a template or import data.

        This method is just a wrapper for calling other methods that generate
        the inputs. The name of the template should be the name of
        the corresponding method (whitespace instead of underscores and capital
        letters are allowed).

        Parameters
        ----------
        template : str
            Name of the template.
        **kwargs
            Parameters passed to the template.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.

        Raises
        ------
        ValueError
            When given template was not found.
        """
        template = template.lower().replace(" ", "_")
        try:
            template = getattr(self, template)
        except AttributeError:
            raise ValueError("Unknown template: " + template)
        return template(**kwargs)
