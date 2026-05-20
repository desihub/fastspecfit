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
from logging import DEBUG
from desiutil.log import get_logger

def getFastspecLogger():
    """Create a fastspecfit-specific logger that writes to stdout.

    Formats messages identically to DESIUtil defaults. Every call returns the
    *same* log object and resets its level to INFO, so this should only be
    called once at program startup.

    Returns
    -------
    log : :class:`logging.Logger`
        Logger named ``fastspec`` with level INFO.

    """
    import sys
    import logging

    root_name = 'fastspec'
    log = logging.getLogger(root_name)

    if not log.handlers:
        ch = logging.StreamHandler(sys.stdout)
        fmtfields = ['%(levelname)s', '%(filename)s', '%(lineno)s', '%(funcName)s']
        fmtfields.append(' %(message)s')
        formatter = logging.Formatter(':'.join(fmtfields),
                                      datefmt='%Y-%m-%dT%H:%M:%S')
        ch.setFormatter(formatter)
        log.addHandler(ch)

    log.setLevel(logging.INFO)

    return log


#log = getFastspecLogger()
log = get_logger()

