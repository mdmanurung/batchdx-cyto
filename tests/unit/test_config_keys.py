"""Every declared configuration key must do something, or say on its face that it does not.

This is the drift gate for :mod:`batchdx.cyto.config`.  The schema once accepted twelve
keys that no other module ever read, so a user could set one, get a clean validation
message, and receive a run that was not configured the way they asked.  The gate walks
every pydantic model in ``config.py`` and requires each leaf to be either

* read by a non-schema module, or
* listed in :data:`DECLARATIVE_ONLY` below, with a justification and with the schema
  accepting exactly one value for it, so an inert key is inert *by constraint* and
  cannot quietly grow a second, ignored value.

The walk is over the models, not over a hard-coded list of names, so adding a key to
``config.py`` fails this test until it is wired up or justified.  Renaming a key fails
it too.

Two modules are excluded from the consumption scan, and both exclusions are load
bearing:

* ``config.py`` itself, because a key cannot consume itself;
* ``batchdx/cyto/schemas/``, because those models are the emitted *contracts*, and a
  field there is a record of a setting rather than a use of one;
* ``batchdx/cyto/examples.py``, which *writes* a configuration.  It names every key
  because it materialises one, which proves nothing about whether anything reads them.
"""

from __future__ import annotations

import ast
import copy
from pathlib import Path

import pydantic
import pytest
import yaml

from batchdx.cyto import config as config_module
from batchdx.cyto._version import CONFIG_SCHEMA_VERSION
from batchdx.cyto.config import (
    SUPPORTED_VALUES,
    UNSUPPORTED_KEYS,
    Config,
    ConfigError,
    parse_config,
)

PACKAGE_ROOT = Path(config_module.__file__).resolve().parent
#: Modules that write or define a configuration rather than act on one.
NOT_CONSUMERS = {
    PACKAGE_ROOT / "config.py",
    PACKAGE_ROOT / "examples.py",
    PACKAGE_ROOT / "schemas",
}

#: Declared keys that no non-schema module reads, each with the reason it is still
#: declared, and with the schema pinning it to the single value the package delivers.
DECLARATIVE_ONLY: dict[str, str] = {
    "preprocessing.cofactor_by_panel": (
        "This package reads an already-preprocessed baseline and never transforms events. "
        "The per-panel cofactors are required to be declared (there is no universal "
        "cofactor) and are read by the schema to refuse an undeclared one; the cofactor "
        "that actually inverts an FCS export is the one in the processing manifest, "
        "which batchdx.cyto.io.export.plan_inversion reads."
    ),
    "design.reference_identity": (
        "Names the sample-manifest column carrying reference identity. The reader "
        "fixes the column at 'reference_id', so the key is pinned to that name (or "
        "unset); it records the declaration that the design audit and the adapters act "
        "on through SampleRecord.reference_id, without a second way to point elsewhere."
    ),
    "design.reference_role": (
        "Names the sample-manifest column carrying reference role. The reader fixes the "
        "column at 'reference_role', so the key is pinned to that name (or unset). The "
        "role itself is read from SampleRecord.reference_role, never from this key: "
        "'biological_pseudo_control' is not technical evidence whatever it says."
    ),
    "correction.control_policy": (
        "Declares which control samples may act as technical evidence. Pinned to "
        "'measured_only': the alternative would have to promote a "
        "biological_pseudo_control, which docs/usage.md and the reference-role contract "
        "forbid whatever the configuration says, so accepting it would be a lie the "
        "design audit could not honour."
    ),
    "report.offline": (
        "The only renderer in this build is one self-contained offline HTML report, so "
        "'false' would name a mode that does not exist. Pinned to true; the key is kept "
        "because a future served renderer has to be told to honour it."
    ),
    "exports.h5ad": (
        "The corrected H5AD is not optional: the metrics rule reads it and would have "
        "nothing to score. Pinned to true."
    ),
}

#: Keys removed outright.  None of them may reappear in DECLARATIVE_ONLY, which claims
#: the key is accepted.
REFUSED_KEYS = frozenset(UNSUPPORTED_KEYS)

#: Declarative keys whose annotation is a free-form mapping, so there is no set of values
#: to pin.  The gate for these is not a value set but that the schema reads them: each
#: one must be refused when the mapping is empty and the thing it declares is in force.
DECLARATIVE_MAPPINGS = {"preprocessing.cofactor_by_panel"}


# -- the model walk ------------------------------------------------------------------------


def config_leaf_paths(model: type[pydantic.BaseModel] = Config, prefix: str = "") -> list[str]:
    """Every leaf field of ``model``, dotted, recursing into nested models.

    Written against ``model_fields`` rather than a fixed list so that a new key appears
    here the moment it is declared.
    """

    out: list[str] = []
    for name, field in model.model_fields.items():
        path = f"{prefix}{name}"
        annotation = field.annotation
        if isinstance(annotation, type) and issubclass(annotation, pydantic.BaseModel):
            out.extend(config_leaf_paths(annotation, f"{path}."))
        else:
            out.append(path)
    return out


def models_declared_in_config_module() -> list[type[pydantic.BaseModel]]:
    """Every ``_Strict`` subclass ``config.py`` declares, plus the models it imports."""

    found = [
        obj
        for obj in vars(config_module).values()
        if isinstance(obj, type) and issubclass(obj, pydantic.BaseModel)
    ]
    return sorted(found, key=lambda m: m.__name__)


def _field_index(model: type[pydantic.BaseModel] = Config, prefix: str = "") -> dict:
    """Dotted leaf path -> the ``FieldInfo`` that declares it."""

    out: dict = {}
    for name, field in model.model_fields.items():
        path = f"{prefix}{name}"
        annotation = field.annotation
        if isinstance(annotation, type) and issubclass(annotation, pydantic.BaseModel):
            out.update(_field_index(annotation, f"{path}."))
        else:
            out[path] = field
    return out


# -- the consumption scan ------------------------------------------------------------------


def _dotted(node: ast.AST) -> tuple[str, ...] | None:
    parts: list[str] = []
    current = node
    while isinstance(current, ast.Attribute):
        parts.append(current.attr)
        current = current.value
    if not isinstance(current, ast.Name):
        return None
    parts.append(current.id)
    parts.reverse()
    return tuple(parts)


def _suffixes(node: ast.AST) -> frozenset[tuple[str, ...]]:
    full = _dotted(node)
    return frozenset(full[i:] for i in range(len(full))) if full else frozenset()


def _ends_with(suffixes: frozenset[tuple[str, ...]], tail: tuple[str, ...]) -> bool:
    return any(s[-len(tail):] == tail for s in suffixes if len(s) >= len(tail))


def _read_module(path: Path) -> tuple[
    frozenset[tuple[str, ...]], frozenset[str], list[tuple[str, tuple[str, ...]]],
    dict[str, frozenset[tuple[str, ...]]],
]:
    """Attribute chains, string literals, ``(attr, base)`` pairs, and local bindings.

    The bindings map is what keeps the scan honest: ``rep = config.evaluation.
    representation`` followed by ``rep.space`` is a real read of
    ``evaluation.representation.space``, whereas an unrelated ``context.design`` and
    ``s.reference_role`` in the same function are not.
    """

    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    chains: set[tuple[str, ...]] = set()
    strings: set[str] = set()
    accesses: list[tuple[str, tuple[str, ...]]] = []
    bindings: dict[str, set[tuple[str, ...]]] = {}

    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute):
            chains |= _suffixes(node)
            base = _dotted(node.value)
            if base is not None:
                accesses.append((node.attr, base))
        elif isinstance(node, ast.Constant) and isinstance(node.value, str):
            strings.add(node.value)
        elif isinstance(node, ast.Assign):
            bound = _suffixes(node.value) if node.value is not None else frozenset()
            for target in node.targets:
                if isinstance(target, ast.Name) and bound:
                    bindings.setdefault(target.id, set()).update(bound)
        elif isinstance(node, ast.AnnAssign) and node.value is not None:
            bound = _suffixes(node.value)
            if isinstance(node.target, ast.Name) and bound:
                bindings.setdefault(node.target.id, set()).update(bound)

    return (
        frozenset(chains),
        frozenset(strings),
        accesses,
        {name: frozenset(value) for name, value in bindings.items()},
    )


def consuming_modules() -> dict[str, tuple]:
    """Parsed non-schema package modules, keyed by path relative to the package root."""

    out: dict[str, tuple] = {}
    for path in sorted(PACKAGE_ROOT.rglob("*.py")):
        if "__pycache__" in path.parts:
            continue
        if any(path == excluded or excluded in path.parents for excluded in NOT_CONSUMERS):
            continue
        out[str(path.relative_to(PACKAGE_ROOT))] = _read_module(path)
    return out


def consumption_evidence(dotted: str, modules: dict[str, tuple] | None = None) -> str | None:
    """Where, and how, ``dotted`` is read.  ``None`` if no non-schema module reads it."""

    modules = consuming_modules() if modules is None else modules
    parts = tuple(dotted.split("."))
    parent, leaf = parts[:-1], parts[-1]
    for name, (chains, strings, accesses, bindings) in modules.items():
        # 1. read straight off the configuration, e.g. ``config.evaluation.labels``
        if parts in chains:
            return f"{name}: attribute {'.'.join(parts)}"
        if not parent:
            continue
        # 2. read off a local alias of the parent, e.g. ``rep.space``
        for attr, base in accesses:
            if attr != leaf:
                continue
            if _ends_with(frozenset({base}), parent):
                return f"{name}: attribute {leaf!r} on {'.'.join(parent)}"
            if len(base) == 1 and parent in bindings.get(base[0], frozenset()):
                return f"{name}: attribute {leaf!r} on the local alias of {'.'.join(parent)}"
        # 3. reached by name through getattr, e.g. ``getattr(config.detectors, flag)``
        if leaf in strings and _ends_with(chains, parent):
            return f"{name}: name {leaf!r} with {'.'.join(parent)} in scope"
    return None


# -- a valid configuration to mutate ------------------------------------------------------


def _base_config(repo_root: Path) -> dict:
    return yaml.safe_load((repo_root / "config" / "example_paired_reference.yaml").read_text())


def _put(raw: dict, dotted: str, value) -> dict:
    raw = copy.deepcopy(raw)
    parts = dotted.split(".")
    block = raw
    for part in parts[:-1]:
        block = block.setdefault(part, {})
    block[parts[-1]] = value
    return raw


def _impossible_values(annotation, allowed: tuple) -> list:
    """Values of the right *type* that the build must refuse.

    Only type-compatible values are probed: a value of the wrong type is rejected by
    pydantic before any of this package's rules run, and asserting on pydantic's wording
    would test pydantic.
    """

    if annotation is bool:
        return [v for v in (True, False) if v not in allowed]
    if annotation is int or annotation is float:
        return [987654321] if 987654321 not in allowed else []
    # Every remaining pinned key in this schema is string-valued (optionally nullable),
    # and a ``mode="before"`` validator sees the raw value, so a nonsense string is a fair
    # probe for all of them.
    return ["a-value-this-build-cannot-honour"]


def _shipped_configs(repo_root: Path) -> list[Path]:
    return sorted((repo_root / "config").glob("example_*.yaml"))


# -- the gate ------------------------------------------------------------------------------


def test_config_leaf_walk_covers_every_model() -> None:
    """Guard the walk itself: a walk that finds nothing proves nothing."""

    leaves = config_leaf_paths()
    assert len(leaves) >= 40, f"the walk found only {len(leaves)} leaf keys; is it still walking?"
    assert sum("." in leaf for leaf in leaves) >= 8, "nesting is no longer being recursed into"
    assert len(leaves) == len(set(leaves)), "a leaf key is reachable by two paths"
    declared = {m.__name__ for m in models_declared_in_config_module()}
    assert {"Config", "EvaluationConfig", "SamplingConfig"} <= declared, (
        f"the models config.py declares are not importable: {sorted(declared)}"
    )


def test_every_config_leaf_is_read_by_a_non_schema_module_or_declared_declarative() -> None:
    """The drift gate: no key may be accepted and then ignored."""

    modules = consuming_modules()
    assert len(modules) > 20, "the consumption scan found suspiciously few modules"
    unaccounted = {}
    for dotted in config_leaf_paths():
        if consumption_evidence(dotted, modules) is not None:
            continue
        if dotted in DECLARATIVE_ONLY:
            continue
        unaccounted[dotted] = DECLARATIVE_ONLY.get(dotted, "absent from DECLARATIVE_ONLY")
    assert not unaccounted, (
        "these configuration keys are accepted by the schema and read by nothing. Add the "
        "read, refuse the key (config.UNSUPPORTED_KEYS), pin the value "
        "(config.SUPPORTED_VALUES), or add it to DECLARATIVE_ONLY in this file with a "
        f"justification: {sorted(unaccounted)}"
    )


def test_declarative_keys_are_not_also_read_by_a_non_schema_module() -> None:
    """A key is either live or declarative.  Both means the justification is stale."""

    modules = consuming_modules()
    stale = {
        dotted: consumption_evidence(dotted, modules)
        for dotted in DECLARATIVE_ONLY
        if consumption_evidence(dotted, modules) is not None
    }
    assert not stale, (
        "these keys are listed as declarative but are read after all; move them out of "
        f"DECLARATIVE_ONLY: {stale}"
    )


def test_declarative_allow_list_covers_no_other_key() -> None:
    """The allow-list is a list, not a net: nothing may be added to it by accident."""

    leaves = set(config_leaf_paths())
    unknown = sorted(set(DECLARATIVE_ONLY) - leaves)
    assert not unknown, f"DECLARATIVE_ONLY names keys that no model declares: {unknown}"
    for dotted, justification in DECLARATIVE_ONLY.items():
        assert justification.strip(), f"{dotted} needs a justification, not a blank string"
        assert len(justification) > 80, f"{dotted} has a justification too short to be one"
    assert not (set(DECLARATIVE_ONLY) & REFUSED_KEYS), (
        "a key cannot be both removed and listed as accepted-and-declarative: "
        f"{sorted(set(DECLARATIVE_ONLY) & REFUSED_KEYS)}"
    )
    assert not (DECLARATIVE_MAPPINGS - set(DECLARATIVE_ONLY)), (
        "DECLARATIVE_MAPPINGS must be a subset of DECLARATIVE_ONLY: "
        f"{sorted(DECLARATIVE_MAPPINGS - set(DECLARATIVE_ONLY))}"
    )


def test_declarative_keys_accept_only_the_value_the_package_delivers(
    repo_root: Path,
) -> None:
    """An inert key must be inert by constraint, not by accident.

    Every declarative key with a value set is pinned in ``config.SUPPORTED_VALUES``, so a
    future edit cannot widen one of them into a second value that is accepted and never
    acted on.  The free-form mappings have no value set to pin; ``cofactor_by_panel`` is
    gated by :func:`test_declared_cofactor_mapping_is_read_not_ignored` instead.
    """

    for dotted in set(DECLARATIVE_ONLY) - DECLARATIVE_MAPPINGS:
        assert dotted in SUPPORTED_VALUES, (
            f"{dotted} is declarative but not pinned in config.SUPPORTED_VALUES, so any "
            "value would be accepted and ignored"
        )
        allowed, why = SUPPORTED_VALUES[dotted]
        assert allowed, f"{dotted} pins no value"
        assert why.strip(), f"{dotted} has no reason recorded for its pin"


def test_declared_cofactor_mapping_is_read_not_ignored(repo_root: Path) -> None:
    """The one declarative mapping is read: an empty one is refused, not tolerated.

    There is no universal cofactor, so ``preprocessing.cofactor_by_panel`` exists to force
    a per-panel declaration.  If nothing read it, an asinh configuration with no cofactor
    would validate and then be scored on an undeclared coordinate system.
    """

    base = _base_config(repo_root)
    raw = _put(base, "preprocessing.enabled", True)
    raw = _put(raw, "preprocessing.transformation", "asinh")
    raw = _put(raw, "preprocessing.cofactor_by_panel", {})
    with pytest.raises(ConfigError) as excinfo:
        parse_config(raw, source="test.yaml")
    assert "preprocessing.cofactor_by_panel" in str(excinfo.value)
    assert "no universal cofactor" in str(excinfo.value)

    accepted = _put(raw, "preprocessing.cofactor_by_panel", {"panel_A": 5.0})
    config = parse_config(accepted, source="test.yaml")
    assert config.preprocessing.cofactor_by_panel == {"panel_A": 5.0}


def test_removed_keys_are_refused_by_name(repo_root: Path) -> None:
    """A removed key is a named error, not 'unknown key'."""

    assert UNSUPPORTED_KEYS, "the removed-key table is empty; is anything still inert?"
    base = _base_config(repo_root)
    for dotted, why in UNSUPPORTED_KEYS.items():
        leaf = dotted.split(".")[-1]
        raw = copy.deepcopy(base)
        raw.setdefault(dotted.split(".")[0], {})[leaf] = "any-value"
        with pytest.raises(ConfigError) as excinfo:
            parse_config(raw, source="test.yaml")
        message = str(excinfo.value)
        assert dotted in message, f"{dotted} is refused without naming the key: {message}"
        assert why[:40] in message, f"{dotted} is refused without its reason: {message}"


def test_removed_keys_are_absent_from_the_models(repo_root: Path) -> None:
    """If a removed key reappears as a field, it must also be consumed or declarative."""

    leaves = set(config_leaf_paths())
    resurrected = sorted(set(UNSUPPORTED_KEYS) & leaves)
    assert not resurrected, (
        f"{resurrected} are in UNSUPPORTED_KEYS but are declared as model fields again; "
        "delete the field or move the entry out of UNSUPPORTED_KEYS"
    )
    del repo_root


def test_unsupported_values_are_refused_by_name(repo_root: Path) -> None:
    """Every key that can only take some values refuses the rest, naming the key."""

    base = _base_config(repo_root)
    field = _field_index()
    for dotted, (allowed, _why) in SUPPORTED_VALUES.items():
        for bad in _impossible_values(field[dotted].annotation, allowed):
            with pytest.raises(ConfigError) as excinfo:
                parse_config(_put(base, dotted, bad), source="test.yaml")
            message = str(excinfo.value)
            assert dotted in message, (
                f"{dotted}={bad!r} was not accepted but the error does not name the key: "
                f"{message}"
            )
            assert "not supported by this build" in message, message


def test_pinned_value_tables_resolve_to_real_fields() -> None:
    """A typo in either table would otherwise silently disable its own check."""

    leaves = set(config_leaf_paths())
    for dotted in SUPPORTED_VALUES:
        assert dotted in leaves, f"SUPPORTED_VALUES names a field that does not exist: {dotted}"
    for dotted in UNSUPPORTED_KEYS:
        section, _, leaf = dotted.partition(".")
        model = Config.model_fields.get(section)
        assert model is not None, f"UNSUPPORTED_KEYS names an unknown section: {dotted}"
        assert "." in dotted, f"UNSUPPORTED_KEYS must name a leaf, not a section: {dotted}"
        annotation = model.annotation
        assert leaf not in annotation.model_fields, (
            f"{dotted} is in UNSUPPORTED_KEYS but is still a declared field"
        )


def test_honoured_values_are_still_accepted(repo_root: Path) -> None:
    """The pins must not refuse the value the shipped examples declare."""

    base = _base_config(repo_root)
    config = parse_config(copy.deepcopy(base), source="test.yaml")
    for dotted, (allowed, _why) in SUPPORTED_VALUES.items():
        value = config
        for part in dotted.split("."):
            value = getattr(value, part)
        assert value in allowed, f"{dotted}={value!r} is not in its own accepted set {allowed}"


def test_config_schema_version_matches_everything_we_ship(repo_root: Path) -> None:
    """The pinned version must be the one the shipped examples declare and write.

    ``CONFIG_SCHEMA_VERSION`` lives in ``batchdx/cyto/_version.py``, which this package
    does not edit; the value it must hold is reported by
    ``test_config_schema_version_is_reported_for_the_record``.
    """

    for path in _shipped_configs(repo_root):
        raw = yaml.safe_load(path.read_text(encoding="utf-8"))
        assert raw["schema_version"] == CONFIG_SCHEMA_VERSION, (
            f"{path.name} declares schema_version {raw['schema_version']} but "
            f"_version.CONFIG_SCHEMA_VERSION is {CONFIG_SCHEMA_VERSION}"
        )
        # And the shipped file must still validate, pins included.
        parse_config(raw, source=str(path))

    from batchdx.cyto.examples import EXAMPLES, _config_dict

    for name, example in EXAMPLES.items():
        assert _config_dict(name, example)["schema_version"] == CONFIG_SCHEMA_VERSION, (
            f"examples.{name} writes schema_version "
            f"{_config_dict(name, example)['schema_version']}, not {CONFIG_SCHEMA_VERSION}"
        )


def test_config_schema_version_is_reported_for_the_record() -> None:
    """Pin the value in a test so that changing it is a deliberate, visible act.

    Raising this requires a decision that also moves ``config/example_*.yaml`` and
    ``batchdx.cyto.examples._config_dict``, which writes the version into every
    materialised project.  See the work-package report.
    """

    assert CONFIG_SCHEMA_VERSION == 2, (
        f"CONFIG_SCHEMA_VERSION is {CONFIG_SCHEMA_VERSION}. This build refuses "
        "preprocessing.gating/peacoqc/time_channel and pins ten values, which narrows the "
        "accepted language of schema_version 2 without changing any shipped example. A "
        "bump to 3 is the semantically honest choice, but batchdx.cyto.examples._config_dict "
        "hard-codes 2 and writes it into every project that 'batchdx-cyto init' creates, so "
        "the bump is a central change that has to move those three files together."
    )


# -- the two bugs this file exists to keep fixed -------------------------------------------


def test_sampling_unit_cannot_be_declared_as_a_unit_that_is_not_sampled(repo_root: Path) -> None:
    """``evaluation.sampling.unit`` used to be recorded and rendered while ignored."""

    base = _base_config(repo_root)
    for unit in ("donor_id", "biospecimen_id", "aliquot_id", "donor"):
        with pytest.raises(ConfigError) as excinfo:
            parse_config(_put(base, "evaluation.sampling.unit", unit), source="test.yaml")
        assert "evaluation.sampling.unit" in str(excinfo.value)
        assert "per sample" in str(excinfo.value)


def test_sampling_manifest_records_the_unit_that_was_actually_used(frozen) -> None:
    """The manifest must not claim a unit the draw did not use."""

    from batchdx.cyto.io.dataset import SAMPLE_ID
    from batchdx.cyto.sampling.partitions import SAMPLE_SAMPLING_UNIT

    _config, dataset, _design, sampling = frozen
    assert sampling.unit == SAMPLE_SAMPLING_UNIT
    assert set(sampling.per_sample_counts) == set(map(str, dataset.obs[SAMPLE_ID].unique()))
    assert set(sampling.aggregation_weights) == set(sampling.per_sample_counts)


def test_declared_uncertainty_unit_is_surfaced_rather_than_silently_dropped(frozen) -> None:
    """A declared uncertainty unit that is never computed has to say so in the report.

    The shipped examples declare ``evaluation.uncertainty_unit: donor_id``.  Stage 8 is
    not implemented, so the manifest has to state that the declared unit bought nothing,
    in a channel the report renders.
    """

    _config, _dataset, _design, sampling = frozen
    assert any(
        "evaluation.uncertainty_unit" in w and "null" in w
        for w in sampling.population_coverage_warnings
    ), sampling.population_coverage_warnings
