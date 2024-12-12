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
    """
    Create a logging object unique to the fastspec tools.
    Configure it to send its log messages to stdout and to
    format them identically to the DESIUtil defaults.

    Note that every call to this function returns the *same*
    log object, so will reset its properties including log level.
    Hence, it should really only be called once at the start
    of the program, as we do here, to create a singleton
    log object.

    Returns a log object with default level INFO.

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

