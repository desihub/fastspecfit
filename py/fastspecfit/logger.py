"""
fastspecfit.logger
==================

Unique logger object, allocated at startup and used throughout code base.

Having just one logger allows us to initialize its level on startup and have
those changes propagate everywhere.

This needs to be in its own file to prevent circular imports with other code.

"""
from desiutil.log import get_logger, DEBUG

log = get_logger()
