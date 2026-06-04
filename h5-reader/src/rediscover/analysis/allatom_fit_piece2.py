#!/usr/bin/env python3
"""Compatibility wrapper for the decomposed all-atom fit pipeline."""

from allatom_fit_common import *  # noqa: F401,F403
from allatom_fit_build3 import *  # noqa: F401,F403
from allatom_fit_legacy import *  # noqa: F401,F403


if __name__ == "__main__":
    build2_main()
