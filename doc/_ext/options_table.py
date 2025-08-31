"""Tools to turn options parameters into RST table.

Example
-------
Following creates a table on the RST file with ``simulation`` parameters:

.. code-block::

   .. options-table:: simulation
"""
import textwrap
from a5py.data import options
from sphinx.util.docutils import SphinxDirective

def generate_options_table(optionsdataclass):
    """Generate RST table from given the given options dataclass."""
    cls = getattr(options, optionsdataclass)
    attr = getattr(options.Options, optionsdataclass)

    lines = []
    lines.append(".. list-table::")
    lines.append("   :header-rows: 1")
    lines.append("")
    lines.append(f"   * - ``Options.{cls.__name__}``")
    lines.append(f"     - {attr.__doc__}")

    for name, member in cls.__dict__.items():
        if not isinstance(member, property):
            continue
        indented_doc = textwrap.indent(member.__doc__, "       ").split("\n")
        lines.append(f"   * - ``{cls.__name__}.{name}``")
        # With literal block we can ensure white space in the generated table.
        lines.append(f"     - | " + indented_doc[0][7:])
        for line in indented_doc[2:]:
            lines.append(f"{line}")
            if line == "":
                lines.append(f"       |")

    return "\n".join(lines)


class OptionsTableDirective(SphinxDirective):
    """Takes one of the options dataclasses as an argument and turns its
    properties to table directive.
    """
    required_arguments = 1

    def run(self):
        class_path = self.arguments[0]
        rst = generate_options_table(class_path)
        self.state_machine.insert_input(rst.splitlines(), self.state.document.current_source)
        return []

def setup(app):
    app.add_directive("options-table", OptionsTableDirective)
    return {"version": "0.1", "parallel_read_safe": True}