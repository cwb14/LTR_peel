"""summarize.py helpers: row roles and the copies-for-80%-correct rule."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import summarize as sm  # noqa: E402


def test_roles():
    assert sm.role("ltrquest") == "primary" and sm.role("hite.intact") == "primary"
    assert sm.role("hite.native") == "secondary" and sm.role("ltrquest_native") == "secondary"
    assert sm.role("ltrquest_hc") == "post_hoc" and sm.role("ltrquest_hc.intact") == "post_hoc"


def test_copies_for_80():
    assert sm.copies_for_80({500: 1.0, 1000: 1.0}, True, 2900) == "<=500"
    assert sm.copies_for_80({500: 0.0, 1000: 0.8, 2000: 1.0}, True, 4800) == "1000"
    assert sm.copies_for_80({500: 0.6, 1000: 1.0, 2000: 0.6, 5000: 1.0}, True, 7000) == "5000"  # must hold above k
    assert sm.copies_for_80({500: 0.0}, True, 914) == "all (914)"
    assert sm.copies_for_80({500: 1.0}, False, 914) == "never"
