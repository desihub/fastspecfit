.. _tutorials:

Tutorials
=========

.. contents:: Contents
    :depth: 2

The following resources provide worked examples for common science tasks using
``FastSpecFit`` and the broader DESI data ecosystem.

Current Notebooks
-----------------

These Jupyter notebooks are included in the ``FastSpecFit`` repository under
``doc/nb/`` and are compatible with the latest ``3.4.0`` release:

- `tutorial-fastspecfit.ipynb`_ — running ``fastspec`` and ``fastphot``
  interactively and reading the results.
- `tutorial-kcorrections.ipynb`_ — computing custom K-corrections and
  rest-frame photometry from the best-fitting SED model.
- `tutorial-vmax.ipynb`_ — estimating :math:`V_{\rm max}` for luminosity
  function calculations.

Legacy Notebooks
----------------

These notebooks target specific public data releases and may not be compatible
with the current version of ``FastSpecFit``.

- `tutorial-fastspec-and-observed-spectra.ipynb`_ — comparing ``fastspec``
  model spectra against observed DESI spectra. Written for the
  :ref:`Fuji VAC (EDR)<fuji vac>` (``FastSpecFit`` v2.x).

DESI Tutorials
--------------

For a broader introduction to working with DESI data — including how to access
data, navigate the file system, and work with DESI spectra — see the `DESI
tutorials`_ repository. The `Getting Started`_ and `Digging Deeper`_
collections are particularly relevant for new users.

.. _`tutorial-fastspecfit.ipynb`: https://github.com/desihub/fastspecfit/blob/main/doc/nb/tutorial-fastspecfit.ipynb
.. _`tutorial-kcorrections.ipynb`: https://github.com/desihub/fastspecfit/blob/main/doc/nb/tutorial-kcorrections.ipynb
.. _`tutorial-vmax.ipynb`: https://github.com/desihub/fastspecfit/blob/main/doc/nb/tutorial-vmax.ipynb
.. _`tutorial-fastspec-and-observed-spectra.ipynb`: https://github.com/desihub/fastspecfit/blob/main/doc/nb/tutorial-fastspec-and-observed-spectra.ipynb
.. _`DESI tutorials`: https://github.com/desihub/tutorials
.. _`Getting Started`: https://github.com/desihub/tutorials/tree/main/01_getting_started
.. _`Digging Deeper`: https://github.com/desihub/tutorials/tree/main/02_digging_deeper
