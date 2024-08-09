try:
    from ._version import version as __version__ # Version file?
except ImportError:
    __version__ = "unknown"
