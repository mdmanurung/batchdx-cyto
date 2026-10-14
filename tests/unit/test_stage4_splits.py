"""Stage 4: donor-level train/tune/validation splits.

The splits exist so that threshold tuning never sees the biospecimens the
validation numbers are computed on. The invariants that matter: donor-level
exclusivity, determinism under the seed, group stratification, and honest
behaviour when there are too few donors to split.
"""

from __future__ import annotations

from batchdx.cyto.sampling.partitions import donor_splits
from batchdx.cyto.schemas.contracts import SampleRecord
from batchdx.cyto.schemas.enums import ReferenceRole


class _FakeSamples:
    def __init__(self, records: dict[str, SampleRecord]):
        self._records = records

    def items(self):
        return self._records.items()


class _FakeDataset:
    def __init__(self, records: list[SampleRecord]):
        self.samples = {r.sample_id: r for r in records}


def _rec(sample_id, batch, group="case", donor=None, reference=None):
    return SampleRecord(
        sample_id=sample_id, batch=batch, group=group, donor_id=donor,
        reference_id=reference,
        reference_role=ReferenceRole("technical_replicate") if reference else ReferenceRole("none"),
    )


def test_every_sample_of_one_donor_lands_in_the_same_split() -> None:
    records = [
        _rec(f"D1_{b}", b, donor="D1") for b in ("B1", "B2", "B3")
    ] + [
        _rec(f"D2_{b}", b, group="control", donor="D2") for b in ("B1", "B2", "B3")
    ] + [
        _rec(f"D3_{b}", b, donor="D3") for b in ("B1", "B2", "B3")
    ]
    splits = donor_splits(_FakeDataset(records), seed=473)
    for donor in ("D1", "D2", "D3"):
        assert sum(donor in splits[s] for s in ("train", "tune", "validation")) == 1


def test_splits_are_deterministic_under_the_seed() -> None:
    records = [
        _rec(f"D{i}_{b}", b, donor=f"D{i}", group="case" if i % 2 else "control")
        for i in range(10) for b in ("B1", "B2")
    ]
    a = donor_splits(_FakeDataset(records), seed=473)
    b = donor_splits(_FakeDataset(records), seed=473)
    assert a == b
    c = donor_splits(_FakeDataset(records), seed=474)
    assert a != c


def test_all_three_splits_are_populated_with_enough_donors() -> None:
    records = [
        _rec(f"D{i}", "B1", donor=f"D{i}", group="case" if i % 2 else "control")
        for i in range(10)
    ]
    splits = donor_splits(_FakeDataset(records), seed=473)
    assert all(len(splits[name]) >= 1 for name in ("train", "tune", "validation"))
    assert len(splits["train"]) >= len(splits["tune"])
    assert len(splits["train"]) >= len(splits["validation"])


def test_reference_replicates_are_excluded_from_splits() -> None:
    records = [
        _rec(f"D{i}", "B1", donor=f"D{i}") for i in range(5)
    ] + [
        _rec(f"ref_{b}", b, group="reference", donor="REFDONOR", reference="REF1")
        for b in ("B1", "B2", "B3")
    ]
    splits = donor_splits(_FakeDataset(records), seed=473)
    everywhere = [d for donors in splits.values() for d in donors]
    assert "REFDONOR" not in everywhere


def test_fewer_than_three_donors_is_reported_unsplit() -> None:
    records = [_rec("D1_B1", "B1", donor="D1"), _rec("D2_B1", "B1", group="control", donor="D2")]
    splits = donor_splits(_FakeDataset(records), seed=473)
    assert splits["unsplit"] == ["D1", "D2"]
    assert splits["train"] == []


def test_without_a_donor_factor_the_sample_is_the_unit() -> None:
    records = [_rec(f"S{i}", "B1") for i in range(6)]
    splits = donor_splits(_FakeDataset(records), seed=473)
    assert sorted(splits["train"] + splits["tune"] + splits["validation"]) == sorted(
        r.sample_id for r in records
    )
    assert splits["unsplit"] == []
