"""Package for generating inputs from templates and imported data.
"""
from .analyticalinputs import AnalyticalInputs
from .optiontemplates  import OptionTemplates
from .poincare  import PoincareTemplates
from .importdata import ImportData
from .convertascot4 import Ascot4Templates

class Template(AnalyticalInputs, OptionTemplates, PoincareTemplates,
                ImportData, Ascot4Templates):
    """Class for creating input data from templates or imported data.

    The templates are constructed by calling :meth:`construct` and specifying
    the template. The templates are constructed in other (public) methods
    of this class. The name of the template is the name of the method, and
    details of the template can be found from the method description.
    """

    def __init__(self, ascot) -> None:
        self._ascot = ascot
        def getmethods(cls):
            """Get private methods of a class."""
            return [attr for attr in dir(cls) if not attr.startswith('_')]

        self._templates = getmethods(AnalyticalInputs)
        self._templates += getmethods(OptionTemplates)
        self._templates += getmethods(PoincareTemplates)
        self._templates += getmethods(ImportData)
        self._templates += getmethods(Ascot4Templates)

    def usetemplate(self, template, **kwargs):
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

    def showtemplate(self, template=None):
        """Show information about a template or a list of all templates.

        Parameters
        ----------
        template : str, optional
            Name of the template.

            If None, all templates are shown.
        """
        if template is None:
            print(self._templates)
        elif template in self._templates:
            print(getattr(self, template).__doc__)
        else:
            raise ValueError("Unknown template: " + template)
