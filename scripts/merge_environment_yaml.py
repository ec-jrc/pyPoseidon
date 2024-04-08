#!/usr/bin/env python

import re
from collections import defaultdict

import click
import srsly


def natural_sort(lst):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split("([0-9]+)", key)]
    return sorted(lst, key=alphanum_key)


@click.command()
@click.argument("environment_files", nargs=-1, required=True, type=click.Path(exists=True))
def main(environment_files):
    """Merge multiple ENVIRONMENT_FILES into a single one.

    ENVIRONMENT_FILES is the paths of the environment files to be merged.
    """
    # Merge data
    data = defaultdict(set)
    for file in environment_files:
        input_data = srsly.read_yaml(file)
        for key, values in input_data.items():
            data[key].update(values)

    #  Convert set to list
    data = dict(data)
    for key in data.keys():
        data[key] = natural_sort(data[key])
    srsly.write_yaml("-", data)


if __name__ == "__main__":
    main()
