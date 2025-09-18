import re
from datetime import datetime


def format2universaldate(date: datetime, withms: bool=False) -> str:
    """Convert a datetime object to a string.

    Parameters
    ----------
    date : datetime.datetime
        The datetime object to be converted.
    withms : bool, *optional*
        If True, include milliseconds.

    Returns
    -------
    formatted_date : str
        The datetime object converted to a string in the format
        "YYYY-MM-DD HH:MM:SS" or "YYYY-MM-DD HH:MM:SS.mmm" if milliseconds are
        included.
    """
    if withms:
        return date.strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]
    return date.strftime("%Y-%m-%d %H:%M:%S")


def decorate(
        string: str,
        color: str | None = None,
        bold: bool = False,
        underline: bool = False,
        ) -> str:
    """Make text underlined, bold, and/or colourful.

    Parameters
    ----------
    string : str
        String to be decorated.
    color : {"green", "purple"}, optional
        Color the text.
    bold : bool, optional
        Bold the text.
    underline : bool, optional
        Underline the text.

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
                f"Available colors are: {colors.keys()}."
                )
        color = colors[color]
    else:
        color = ""

    bolded = effects["bold"] if bold else ""
    underlined = effects["underline"] if underline else ""
    return f"{bolded}{underlined}{color}{string}{effects['reset']}"


def undecorate(string: str) -> str:
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
