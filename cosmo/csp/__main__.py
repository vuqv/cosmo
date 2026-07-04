"""Entry point for ``python -m cosmo.csp``.

Forwards to the O'Brien continuous-synthesis (3-stage) runner.
"""
from .protocol import csp

if __name__ == "__main__":
    csp()
