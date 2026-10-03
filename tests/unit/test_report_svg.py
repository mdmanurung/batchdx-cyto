"""Report figures must keep their bytes across a resumed workflow run."""

import re

import numpy as np

from batchdx.cyto.reporting.plots import density_by_batch


def test_inline_svg_is_repeatable_and_has_figure_specific_ids() -> None:
    values = {"batch_a": np.array([1.0, 2.0, 3.0])}
    first = density_by_batch(values, "marker_a", "expression")
    repeated = density_by_batch(values, "marker_a", "expression")
    other = density_by_batch(values, "marker_b", "expression")

    assert first == repeated
    assert "<dc:date>" not in first
    first_ids = set(re.findall(r'id="(m[0-9a-f]+)"', first))
    other_ids = set(re.findall(r'id="(m[0-9a-f]+)"', other))
    assert first_ids and other_ids
    assert first_ids.isdisjoint(other_ids)
