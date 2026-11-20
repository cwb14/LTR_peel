"""LTR_peel - iterative detection of nested LTR retrotransposons.

Each round masks the elements it found and re-runs, so an element sitting inside
another element surfaces one layer at a time. The rounds are then reconciled
into depth buckets: ``depth0`` holds elements with nothing nested inside them,
``depth1`` holds elements with one layer inside, and so on.

The pipeline is driven by the ``ltr_peel`` console script (see
:mod:`ltr_peel.cli`). Every stage is also a module with its own ``--help`` and
can be run standalone:

============================  ==============================================
``python -m ltr_peel.detect``     one detection round on one genome
``python -m ltr_peel.mask``       mask detected elements into the next round's genome
``python -m ltr_peel.reconcile``  pool rounds into depth-bucketed tables
``python -m ltr_peel.annotate``   add ``strand`` and ``family`` columns
``python -m ltr_peel.gff3``       write the pooled GFF3
============================  ==============================================
"""

__version__ = "1.0.1"

__all__ = ["__version__"]
