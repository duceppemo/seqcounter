#!/usr/bin/env python3

"""
Kept so that `python seqcounter.py` keeps working from a clone, which is how the first version of seqcounter was
run and how the README documented it. The code lives in the seqcounter package, and `seqcounter` is the installed
entry point; both take the same arguments.

The import is inside the guard on purpose: this file has the same name as the package, so a tool that imports it
by path (setuptools reading the version, for one) must not run anything.
"""

__author__ = 'duceppemo'


if __name__ == '__main__':
    from seqcounter.cli import main

    main()
