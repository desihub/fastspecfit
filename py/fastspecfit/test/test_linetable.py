"""
fastspecfit.test.test_linetable
===============================

"""
import os
import pytest


def test_linetable():
    import astropy
    from fastspecfit.linetable import LineTable
    emline_table = LineTable(emlines_file=None)
    assert(os.path.isfile(emline_table.file))
    assert(type(emline_table.table) == astropy.table.Table)

    with pytest.raises(FileNotFoundError):
        LineTable(emlines_file='doesnotexist.ecsv')
