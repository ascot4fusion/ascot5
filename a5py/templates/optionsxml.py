import textwrap
import json
import warnings
from xml.etree.ElementTree import ElementTree
import xml.etree.ElementTree as ET
from xml.dom import minidom

import xmlschema


from a5py.data import options
from .template import InputTemplate

from dataclasses import dataclass

XS_NS = "http://www.w3.org/2001/XMLSchema"

@dataclass
class SimpleType():
    """SimpleType used in the XML schema."""
    name: str
    """Name of the simpleType e.g. "IntegerPositive"."""
    base: str
    """Base type of the simpleType e.g. "xs:integer"."""
    min_val: int | float = None
    """Minimum value of the simpleType if applicable."""
    max_val: int | float = None
    """Maximum value of the simpleType if applicable."""
    min_inclusive: bool = True
    """Whether the minimum value is inclusive or exclusive."""
    max_inclusive: bool = True
    """Whether the maximum value is inclusive or exclusive."""
    array: str = None
    """If not None, the simpleType is an array with elements of this type."""
    union: str = None
    """If not None, the simpleType is a union of these types."""

custom_types = [
    SimpleType("IntegerPositive", "integer", min_val=1),
    SimpleType("IntegerBinary", "integer", min_val=1, max_val=2),
    SimpleType("Integer012", "integer", min_val=0, max_val=2),
    SimpleType("Integer1234", "integer", min_val=1, max_val=4),
    SimpleType("FloatPositive", "float", min_val=0., min_inclusive=False),
    SimpleType("FloatNonNegative", "float", min_val=0.),
    SimpleType("Float0to360", "float", min_val=0., max_val=360.),
    SimpleType("FloatNNList", "FloatNonNegative", array="FloatNonNegative"),
    SimpleType("FloatOrFloatNNList", "FloatOrFloatNNList", array="FloatOrFloatNNList"),
]
"""List of custom simpleTypes that are used in the XML schema."""

parameter_types = {
    "simulation_mode": "Integer1234",
    "record_mode": "IntegerBinary",
    "use_explicit_fixedstep": "IntegerBinary",
    "explicit_fixedstep": "FloatPositive",
    "gyrodefined_fixedstep": "IntegerPositive",
    "enable_adaptive": "IntegerBinary",
    "adaptive_tolerance_orbit": "FloatPositive",
    "adaptive_tolerance_collisions": "FloatPositive",
    "adaptive_max_drho": "FloatPositive",
    "adaptive_max_dphi": "FloatPositive",
    "enable_orbit_following": "IntegerBinary",
    "enable_coulomb_collisions": "IntegerBinary",
    "enable_mhd": "IntegerBinary",
    "enable_atomic": "Integer012",
    "enable_icrh": "IntegerBinary",
    "enable_aldforce": "IntegerBinary",
    "disable_first_order_gctransformation": "IntegerBinary",
    "disable_energy_ccollision": "IntegerBinary",
    "disable_pitch_ccollision": "IntegerBinary",
    "disable_gcdiff_ccollision": "IntegerBinary",
    "reverse_time": "IntegerBinary",
    "activate_simulation_time_limits": "IntegerBinary",
    "activate_real_time_limits": "IntegerBinary",
    "activate_rho_limit": "IntegerBinary",
    "activate_energy_limits": "IntegerBinary",
    "activate_orbit_limit": "Integer012",
    "activate_neutralization": "IntegerBinary",
    "activate_ionization": "IntegerBinary",
    "lab_time_limit": "xs:float",
    "max_mileage": "FloatPositive",
    "max_real_time": "FloatPositive",
    "rho_coordinate_limits": "2xFloatPositive",
    "min_energy": "FloatPositive",
    "min_local_thermal_energy": "FloatPositive",
    "max_number_of_toroidal_orbits": "FloatPositive",
    "max_number_of_poloidal_orbits": "FloatPositive",
    "collect_dist5d": "IntegerBinary",
    "collect_dist6d": "IntegerBinary",
    "collect_dist5drho": "IntegerBinary",
    "collect_dist6drho": "IntegerBinary",
    "collect_distcom": "IntegerBinary",
    "r_bins": "IntegerPositive",
    "phi_bins": "IntegerPositive",
    "z_bins": "IntegerPositive",
    "rho_bins": "IntegerPositive",
    "theta_bins": "IntegerPositive",
    "ppara_bins": "IntegerPositive",
    "pperp_bins": "IntegerPositive",
    "pr_bins": "IntegerPositive",
    "pphi_bins": "IntegerPositive",
    "pz_bins": "IntegerPositive",
    "time_bins": "IntegerPositive",
    "charge_bins": "IntegerPositive",
    "r_interval": "2xFloatPositive",
    "phi_interval": "2xFloat0to360",
    "z_interval": "2xFloat",
    "rho_interval": "2xFloatPositive",
    "theta_interval": "2xFloat0to360",
    "ppara_interval": "2xFloat",
    "pperp_interval": "2xFloat",
    "pr_interval": "2xFloat",
    "pphi_interval": "2xFloat",
    "pz_interval": "2xFloat",
    "time_interval": "2xFloat",
    "charge_interval": "2xIntegerPositive",
    "collect_orbit": "IntegerBinary",
    "collect_transport_coefficient": "IntegerBinary",
}
"""Mapping from options paramateres to their types in XML.

Unfortunately this cannot be generated automatically.
"""


class ReadOptionsXml(InputTemplate):

    def __init__(self, ascot, file):
        """Convert XML file to :class:`Opt` input.

        Parameters
        ----------
        file : str
            Either path to the XML file or it's contents as a string.
        """
        def load_xml():
            try:
                return ElementTree().parse(file)
            except:
                pass
            try:
                return ElementTree.fromstring(file)
            except:
                pass
            raise ValueError("Could not parse XML file or string.")

        tree = load_xml()
        opt = Opt.get_default()
        for tag in opt.keys():
            val = tree.findall("*"+tag)
            if len(val) == 0:
                warnings.warn(
                    f"Using default value for {tag} as it was not found in xml"
                    )
                continue
            opt[tag] = val[0].text
        return opt


class WriteOptionsXml():
    pass

def schema(opt=None, fnxsd=None, fnxml=None):
    """Produce a schema and an XML file containing ASCOT5 option parameters.

    Parameters
    ----------
    opt : dict, optional
        Options used to fill the XML with parameter values.

        If None, the default options are used.
    fnxsd : str, optional
        None or filename to write XSD data.
    fnxml : str, optional
        None or filename to write XML data.

    Returns
    -------
    schema :
        Schema that contains parameter descriptions and which can validate
        XML files containing ASCOT5 options.
    xml :
        XML file containing the options parameters and their values.
    """
    if xmlschema is None:
        raise ImportError(
            "This feature requires package xmlschema:\n"
            "pip install xmlschema")

    schema = xmlschema.XMLSchema(xsd)
    if opt is None: opt = Opt.get_default()
    opt_hierarchy = {
        "SIMULATION_MODE_AND_TIMESTEP":{}, "END_CONDITIONS":{},
        "PHYSICS":{}, "DISTRIBUTIONS":{}, "ORBIT_WRITE":{},
        "TRANSPORT_COEFFICIENT":{}}
    for k, v in opt.items():
        if k == "SIM_MODE":
            grp = "SIMULATION_MODE_AND_TIMESTEP"
        elif k == "ENDCOND_SIMTIMELIM":
            grp = "END_CONDITIONS"
        elif k == "ENABLE_ORBIT_FOLLOWING":
            grp = "PHYSICS"
        elif k == "ENABLE_DIST_5D":
            grp = "DISTRIBUTIONS"
        elif k == "ENABLE_ORBITWRITE":
            grp = "ORBIT_WRITE"
        elif k == "ENABLE_TRANSCOEF":
            grp = "TRANSPORT_COEFFICIENT"
        opt_hierarchy[grp][k] = v

    data = json.dumps({"parameters":opt_hierarchy})
    xml = xmlschema.from_json(data, schema=schema, preserve_root=True)

    if fnxsd is not None:
        with open(fnxsd, 'w') as f:
            f.writelines(xsd)
    if fnxml is not None:
        ElementTree(xml).write(fnxml)

    return schema, xml


def constraint_to_xsd(constraint: str, xtype: str = "xs:float") -> str:
    kind, spec = parse_constraint(constraint)

    if kind == "set":
        enums = "\n".join(f'    <xs:enumeration value="{v}"/>' for v in spec)
        return f"<xs:restriction base=\"{xtype}\">\n{enums}\n</xs:restriction>"

    if kind == "range":
        m = re.match(r"([<>]=?|==)\s*(-?\d+\.?\d*)", spec)
        if not m:
            raise ValueError("Unsupported range")
        op, val = m.groups()
        tag = {"<": "maxExclusive", "<=": "maxInclusive", ">": "minExclusive", ">=": "minInclusive"}[op]
        return f"<xs:restriction base=\"{xtype}\">\n  <xs:{tag} value=\"{val}\"/>\n</xs:restriction>"

    if kind == "tuple":
        elems = []
        for i, s in enumerate(spec):
            if "..." in s:  # unbounded repetition
                s = spec[0]
            fragment = constraint_to_xsd(f"({s})", xtype)
            elems.append(f"<xs:element name=\"item{i}\" type=\"xs:{xtype[3:]}\"/> <!-- {s} -->")
        return "<xs:sequence>\n" + "\n".join(elems) + "\n</xs:sequence>"

    raise ValueError(f"Unknown constraint: {constraint}")


def doc(name, stype):
    """Helper function to add parameters with docstrings, types, and
    proper indentation.

    Paramaters
    ----------
    name : str
        Name of the parameter as it appears in option attributes
        (without preceding underscore).
    stype : str
        This parameters simple type.

    Returns
    -------
    fstr
        Formatted string that specifies XML element.
    """
    parameter = options.simulation.__dict__[name]
    docstring = parameter.__doc__.replace('<', '&lt;').replace('>', '&gt;')
    text = f"""
        <xs:element name="{name}" type="{stype}">
            <xs:annotation>
            <xs:documentation>
            {textwrap.indent(docstring, "            ")}
            </xs:documentation>
            </xs:annotation>
        </xs:element>"""
    return text


def prettify(elem):
    """Return a pretty-printed XML string for the Element."""
    rough_string = ET.tostring(elem, 'utf-8')
    return minidom.parseString(rough_string).toprettyxml(indent=4*" ")

def make_simple_type(schema: ET.Element, simple_type: SimpleType):
    """Add simpleType element to schema.

    Parameters
    ----------
    schema : Element
        The XML schema to be appended.
    simple_type : SimpleType
        SimpleType to add to schema
    """
    element = ET.SubElement(schema, f"xs:simpleType", name=simple_type.name)
    restriction = ET.SubElement(element, f"xs:restriction", base=f"{simple_type.base}")

    if simple_type.min_val is not None:
        tag = "minInclusive" if simple_type.min_inclusive else "minExclusive"
        ET.SubElement(restriction, f"xs:{tag}", value=str(simple_type.min_val))
    if simple_type.max_val is not None:
        tag = "maxInclusive" if simple_type.max_inclusive else "maxExclusive"
        ET.SubElement(restriction, f"xs:{tag}", value=str(simple_type.max_val))
    if simple_type.array is not None:
        ET.SubElement(element, f"xs:list", itemType=simple_type.array)
    if simple_type.union is not None:
        ET.SubElement(element, f"xs:union", memberTypes=simple_type.union)


def make_element_block(schema, container_name, container_doc=None, params=None):
    """
    Create an XML schema element block with parameters.

    Returns an ElementTree.Element
    """
    element = ET.SubElement(schema, f"xs:element", name=container_name)

    if container_doc is  not None:
        annotation = ET.SubElement(element, f"xs:annotation")
        doc = ET.SubElement(annotation, f"xs:documentation")
        doc.text = container_doc

    complex_type = ET.SubElement(element, f"xs:complexType")
    all_elem = ET.SubElement(complex_type, f"xs:all")

    if params is not None:
        for param in params:
            ET.SubElement(all_elem, f"xs:element", ref=param)

    #return element


def make_schema():
    """Create XML schema dynamically from  the options module."""
    schema = ET.Element(f"xs:schema", attrib={"xmlns:xs":XS_NS})

    for simple_type in custom_types:
        make_simple_type(schema, simple_type)

    sections = [
        name for name, attr in
        vars(options.Options).items() if isinstance(attr, property)
        ]
    make_element_block(schema, "parameters", params=sections)
    for section in sections:
        cls = getattr(options, section)
        params = [
            name for name, attr in
            vars(cls).items() if isinstance(attr, property)
            ]
        make_element_block(schema, section, cls.__doc__, params=params)
        for param in params:
            sub = ET.SubElement(schema, parameter_types[param], name=param)
            annotation = ET.SubElement(sub, f"xs:annotation")
            doc = ET.SubElement(annotation, f"xs:documentation")
            doc.text = getattr(cls, param).__doc__
    xml_string = prettify(schema)
    print(xml_string)
    return xml_string
