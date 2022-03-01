#!/usr/bin/env python

import re
import sys
from collections import defaultdict

import click
import ruamel.yaml


yaml = ruamel.yaml.YAML()
yaml.default_flow_style = False
yaml.indent(mapping=2, sequence=4, offset=2)


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split("([0-9]+)", key)]
    return sorted(l, key=alphanum_key)


@click.command()
@click.argument("environment_files", nargs=-1, required=True, type=click.Path(exists=True))
def main(environment_files):
    """Merge multiple ENVIRONMENT_FILES into a single one.

    ENVIRONMENT_FILES is the paths of the environment files to be merged.
    """
    # Merge data
    data = defaultdict(set)
    for file in environment_files:
        with open(file, "r") as fd:
            text = fd.read()
        input_data = yaml.load(text)
        for key, values in input_data.items():
            data[key].update(values)

    #  Convert set to list
    data = dict(data)
    for key in data.keys():
        data[key] = natural_sort(data[key])

    yaml.dump(data, sys.stdout)


if __name__ == "__main__":
    main()
