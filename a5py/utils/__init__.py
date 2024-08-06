"""Utility tools that are i) not related to physics, ii) not specific to ASCOT5,
and iii) which don't fit anywhere else.
"""
import re

def format2universaldate(date):
    """Convert a datetime object to a string in the format ASCOT5 uses.

    Parameters
    ----------
    date : datetime.datetime
        The datetime object to be converted.

    Returns
    -------
    formatted_date : str
        The datetime object converted to a string in the format
        "YYYY-MM-DD HH:MM:SS".
    """
    return date.strftime("%Y-%m-%d %H:%M:%S")

def decorate(string, color=None, bold=False, underline=False):
    """Make text underlined, bold, and/or colourful.

    Parameters
    ----------
    string : str
        String to be decorated.

    Returns
    -------
    decorated : str
        The decorated string.
    """
    colors = {
        "green":"\033[32m",
        "purple":"\033[35m",
    }
    effects = {
        "reset":"\033[0m",
        "bold":"\033[01m",
        "underline":"\033[04m",
    }

    if color:
        if color not in colors:
            raise ValueError(
                f"Unknown color: {color}. "
                f"Available colors are: {colors.keys()}.")
        color = colors[color]
    else:
        color = ""

    bold = effects["bold"] if bold else ""
    underline = effects["underline"] if underline else ""
    return f"{bold}{underline}{color}{string}{effects["reset"]}"

def undecorate(string):
    """Remove decorations (ANSI escape sequences) from a string.

    Parameters
    ----------
    string : str
        A string to undecorate.

    Returns
    -------
    undecorated : str
        The undecorated string.
    """
    return re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])").sub("", string)
