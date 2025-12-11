#!/usr/bin/env python
import argparse
import ipyrad2 as ip
import os
import tomlkit

from ipyrad2.utils.exceptions import IPyradError
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, List


def read_params(paramsfile: str):
    """

    """
    with open(paramsfile) as infile:
        params = tomlkit.load(infile)

    params = _replace_values(params.unwrap(), -1, None)
    params = RecursiveNamespace(**params)
    return params


def new_params(name: str = None, force: bool = False):
    """ 
    Write out the parameters of this assembly to a file properly
    formatted as input for `ipyrad -p <params.txt>`. A good and
    simple way to share/archive parameter settings for assemblies.
    This is also the function that's used by __main__ to
    generate default params.txt files for `ipyrad -n`
    """
    if name is None:
        outfile = "params-default.txt"
    else:
        outfile = f"params-{name}.txt"

    ## Test if params file already exists?
    ## If not forcing, test for file and bail out if it exists
    if not force:
        if os.path.isfile(outfile):
            raise IPyradError(PARAMS_EXISTS.format(outfile))

    with open(outfile, 'w') as paramsfile:
        ## Write the header. Format to 80 columns
        header = "# ------- ipyrad params file (v.{})".format(ip.__version__)
        header += ("-" * (80 - len(header))) + "\n"
        paramsfile.write(header + "\n")

        doc = tomlkit.document()
        main = {"main":{"name":name,
                        "project_dir":"./",
                        "raw_fastq_path":"/path/to/fastqs/*.gz",
                        "barcodes_path":"/path/to/bcodes.txt",
                        "sorted_fastq_path":"/path/to/sorted_fastqs/*.gz",
                        "reference_sequence":"/path/to/ref.fa",
                        "pop_assign_file":"/path/to/pops.txt",
                       }
               }

        doc.update(main)
        paramsfile.write(tomlkit.dumps(doc) + "\n")

        parser = ip.cli.cli_main.setup_parsers()
        arg_defaults = _get_arg_defaults(parser)

        subcommands = ["demux", "trim", "denovo", "map", "assemble"]
        for command in subcommands :
            paramsfile.write(f"[{command}]\n")
            args = arg_defaults[command]
            # demux - Remove barcodes which is a classic mode main param
            _ = args.pop("barcodes", None)
            # demux/trim/denovo/map - Remove fastqs which either comes in from main
            #   params or is set as default in classic mode
            _ = args.pop("fastqs", None)
            # map/assemble - Remove reference as it is either passed in or 
            #   constructed in classic mode
            _ = args.pop("reference", None)
            # assemble - Remove rad_bams/wgs_bams which are constructed in classic mode
            _ = args.pop("rad_bams", None)
            _ = args.pop("wgs_bams", None)
            # assemble - Remove name which is passed in to classic mode
            _ = args.pop("name", None)
            # Remove all output directories. Classic mode will use force to use
            # the defaults, so all classic mode directories are determined solely
            # by the project_dir and name parameters
            _ = args.pop("out", None)
            # Remove other args that don't make sense in a params file
            _ = args.pop("cores", None)
            _ = args.pop("threads", None)
            _ = args.pop("force", None)
            _ = args.pop("log_level", None)
            _ = args.pop("log_file", None)
            doc = tomlkit.document()
            doc.update(args)

            # Dump the TOMLDocument to a string
            toml_string = tomlkit.dumps(doc)
            paramsfile.write(toml_string + "\n")


class RecursiveNamespace(SimpleNamespace):
    """
    A helper class for reformatting a toml document into a Namespace
    object so that dictionary keys can be accessed as properties, so
    the behavior aligns better with the output of argparse.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        for key, val in kwargs.items():
            if isinstance(val, dict):
                setattr(self, key, RecursiveNamespace(**val))
            elif isinstance(val, list):
                setattr(self, key, [self._map_entry(item) for item in val])

    def _map_entry(self, entry):
        if isinstance(entry, dict):
            return RecursiveNamespace(**entry)
        elif isinstance(entry, list):
            return [self._map_entry(item) for item in entry]
        return entry


def _replace_values(d, old, new):
    """
    Helper function for reading in params and replacing `-1` values
    with None. toml doesn't have a native None type.
    """
    for k, v in d.items():
        if isinstance(v, dict):
            _replace_values(v, old, new)
        elif v == old:
            d[k] = new
    return d


def _get_arg_defaults(parser):
    """
    Recursively extract all argument names and their default values
    from an argparse.ArgumentParser with subparsers.
    """

    arg_defaults = {}
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            for subcmd, subparser in action.choices.items():
                args = {}
                for action in subparser._actions:
                    if not action.option_strings and action.nargs == 0:
                        continue  # skip positional without defaults
                    dest = action.dest
                    if dest not in ('help', '==SUPPRESS=='):
                        args[dest] = action.default
                arg_defaults[subcmd] = args
    return arg_defaults


def _path_encoder(obj):
    """
    Encodes a pathlib.Path object into a TOML string. Also handles
    arguments with default None type by setting to -1.
    """
    if isinstance(obj, Path):
        # Convert the Path object to its string representation
        # and create a tomlkit.String item.
        return tomlkit.items.String.from_raw(str(obj))
    elif obj is None:
        return tomlkit.integer(-1)
    else:
        # If the object is not a Path, raise ConvertError so
        # tomlkit can continue with the next encoder.
        raise tomlkit.exceptions.ConvertError

_ = tomlkit.register_encoder(_path_encoder)


PARAMS_EXISTS = """
    Error: Params file already exists: {}
    Use force argument to overwrite.
    """

if __name__ == "__main__":

        new_params(force=True)

# A different way to parse out params, but i think it's worse
#def new_params(paramsfile: str):
#
#    with open(paramsfile, 'w') as outfile:
#        subcommands = ["demux", "trim", "denovo", "map", "assemble"]
#        parser = setup_parsers()
#        required_args = {"demux":["-d", "./wat", "-b", "bcodes.txt"],
#                         "trim":["-d", "./wat"],
#                         "denovo":["-d", "./wat"],
#                         "map":["-d", "./wat", "-r", "ref/"],
#                         "assemble":["-d", "./wat", "-r", "ref/"]}
#        for command in subcommands :
#            args_list = [command] + required_args[command]
#            print(args_list)
#            args = parser.parse_args(args_list)
#            print(args)
#            outfile.write(f"[{command}]\n")
#            doc = tomlkit.document()
#            doc.update(args.__dict__)
#
#            # Dump the TOMLDocument to a string
#            toml_string = tomlkit.dumps(doc)
#            outfile.write(toml_string + "\n")
