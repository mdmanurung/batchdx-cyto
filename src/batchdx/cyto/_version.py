"""Single source of the package version, mirrored from ``pyproject.toml``."""

__version__ = "0.1.0.dev0"

#: Version of the on-disk contract used by every emitted record.  Bump this when a
#: record gains or loses a field so that old result directories are rejected loudly
#: rather than silently misread.
CONTRACT_SCHEMA_VERSION = 1

#: Version of the user-facing configuration file format.
CONFIG_SCHEMA_VERSION = 2
