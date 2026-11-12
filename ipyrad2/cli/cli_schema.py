"""Machine-readable documentation for CLI arguments."""

import argparse
import json
from pathlib import Path

from ipyrad2 import __version__

from .cli_demux import DEMUX_DEPENDENCIES, DEMUX_VALUE_CONSTRAINTS
from .common import intlike


def print_demux_schema(parser: argparse.ArgumentParser) -> None:
    """Print a JSON description of the demux argument parser."""
    types = {int: "integer", intlike: "integer", str: "string", Path: "path"}
    formatter = parser._get_formatter()
    options = {}
    for group in parser._action_groups:
        for action in group._group_actions:
            if action.help == argparse.SUPPRESS:
                continue
            option = {
                "flags": action.option_strings,
                "group": group.title,
                "description": formatter._expand_help(action),
                "required": action.required,
                "default": action.default,
                "nargs": action.nargs if action.nargs is not None else 1,
            }
            if isinstance(action, argparse._HelpAction):
                option.update(type="flag", action="help")
            elif isinstance(action, argparse._StoreTrueAction):
                option.update(type="boolean", action="store_true")
            else:
                option["type"] = types[action.type]
                if action.type is intlike:
                    option["input_format"] = "Numeric text rounded to an integer, e.g., 2e7."
            if action.choices is not None:
                option["choices"] = list(action.choices)
            if action.dest in DEMUX_VALUE_CONSTRAINTS:
                option.update(DEMUX_VALUE_CONSTRAINTS[action.dest][0])
            options[action.dest] = option

    print(json.dumps({
        "schema_version": 1,
        "ipyrad2_version": __version__,
        "command": ["ipyrad2", "demux"],
        "options": options,
        "dependencies": [
            {"option": name, relation: other, "description": message}
            for name, relation, other, message in DEMUX_DEPENDENCIES
        ],
    }, indent=2))
