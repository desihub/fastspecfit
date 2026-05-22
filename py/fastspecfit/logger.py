"""
fastspecfit.logger
==================

Unique logger object, distinct from the one used by the DESI libraries,
allocated at startup and used throughout code base.  We need to create
our own logger to avoid the DESI libraries accidentally changing the
log level out from under us.

Having just one logger allows us to initialize its level on startup and have
those changes propagate everywhere.

This needs to be in its own file to prevent circular imports with other code.

"""
import sys
import logging
from logging import DEBUG


class _DynamicStdoutHandler(logging.StreamHandler):
    """StreamHandler that always writes to the current sys.stdout.

    Rebinds self.stream on every emit() so that OS-level stdout redirections
    (e.g. desispec.io.util.stdouterr_redirected) are followed correctly.
    """
    def emit(self, record):
        self.stream = sys.stdout
        super().emit(record)


def getFastspecLogger():
    """Create a fastspecfit-specific logger that writes to stdout.

    Returns the same logger object on every call without resetting its level,
    so external code (e.g. sc_data.initialize) can safely call setLevel() and
    have that setting persist.

    Returns
    -------
    log : :class:`logging.Logger`
        Logger named ``fastspec``.

    """
    log = logging.getLogger('fastspec')

    if not log.handlers:
        ch = _DynamicStdoutHandler()
        fmtfields = ['%(levelname)s', '%(filename)s', '%(lineno)s', '%(funcName)s']
        fmtfields.append(' %(message)s')
        formatter = logging.Formatter(':'.join(fmtfields),
                                      datefmt='%Y-%m-%dT%H:%M:%S')
        ch.setFormatter(formatter)
        log.addHandler(ch)
        log.setLevel(logging.INFO)

    return log


log = getFastspecLogger()
