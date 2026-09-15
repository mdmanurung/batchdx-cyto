"""batchdx-cyto: batch-effect diagnosis, correction and evaluation for cytometry.

The package deliberately keeps *analysis* in importable, testable Python modules and
leaves *scheduling* to Snakemake.  Nothing in this package implements a second
workflow engine.
"""

from batchdx.cyto._version import __version__

__all__ = ["__version__"]
