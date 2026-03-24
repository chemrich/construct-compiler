"""
Parametric variant generator for design space exploration.

Given a base YAML spec and a set of axes to vary, generates all
combinations as spec dicts ready for evaluate_batch().

This is the "design agent" half of the harness — it produces the
candidate specs, while harness.evaluate_batch() scores them.

Usage:

    from construct_compiler.validation.variants import vary_spec, DesignAxis

    axes = [
        DesignAxis("rbs", field_path="cassette.0.cistron.expression",
                   values=["high", "medium", "low"]),
        DesignAxis("spacer", field_path="cassette.2.spacer",
                   values=[20, 30, 50, 100]),
        DesignAxis("terminator", field_path="cassette.4.terminator",
                   values=["rrnB_T1", "T7term", "BBa_B0015"]),
    ]

    variants = vary_spec("examples/his_tev_mbp_egfp.yaml", axes)
    # → 3 × 4 × 3 = 36 variant spec dicts

    from construct_compiler.validation.harness import evaluate_batch
    results = evaluate_batch(variants, skip_constraints=True)
"""

from __future__ import annotations

import copy
import itertools
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Union

import yaml


@dataclass
class DesignAxis:
    """
    One axis of variation in the design space.

    Attributes:
        name: human-readable label for this axis (e.g., "RBS strength")
        field_path: dot-separated path into the cassette list, using
                    integer indices for list positions. See _set_nested().
        values: list of values to substitute at that path
    """
    name: str
    field_path: str
    values: list[Any]


@dataclass
class VariantSpec:
    """A generated variant with its parameter choices annotated."""
    spec: dict
    variant_name: str
    parameters: dict[str, Any]  # axis_name → chosen value


def vary_spec(
    base_spec: Union[str, Path, dict],
    axes: list[DesignAxis],
    name_template: str = "{base}_{params}",
) -> list[VariantSpec]:
    """
    Generate all combinations of design axes applied to a base spec.

    Args:
        base_spec: path to YAML file, YAML string, or parsed dict
        axes: list of DesignAxis defining the search space
        name_template: how to name variants. {base} is the original name,
                       {params} is a short string of parameter values.

    Returns:
        List of VariantSpec objects, one per combination
    """
    # Load base spec
    if isinstance(base_spec, (str, Path)):
        path = Path(base_spec)
        if path.exists():
            base = yaml.safe_load(path.read_text())
        else:
            base = yaml.safe_load(base_spec)
    else:
        base = base_spec

    root = base.get("construct", base)
    base_name = root.get("name", "variant")

    # Generate all combinations
    axis_values = [axis.values for axis in axes]
    axis_names = [axis.name for axis in axes]
    field_paths = [axis.field_path for axis in axes]

    variants = []
    for combo in itertools.product(*axis_values):
        # Deep copy and apply substitutions
        variant = copy.deepcopy(base)
        vroot = variant.get("construct", variant)
        params = {}

        for axis, path, value in zip(axes, field_paths, combo):
            _set_nested(vroot, path, value)
            params[axis.name] = value

        # Generate name
        param_str = "_".join(f"{k}={v}" for k, v in params.items())
        vname = name_template.format(base=base_name, params=param_str)
        vroot["name"] = vname

        variants.append(VariantSpec(
            spec=variant,
            variant_name=vname,
            parameters=params,
        ))

    return variants


def vary_spec_dicts(
    base_spec: Union[str, Path, dict],
    axes: list[DesignAxis],
) -> list[dict]:
    """
    Convenience wrapper: returns just the spec dicts (no metadata).
    Suitable for passing directly to evaluate_batch().
    """
    return [v.spec for v in vary_spec(base_spec, axes)]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _set_nested(obj: Any, path: str, value: Any) -> None:
    """
    Set a value in a nested dict/list structure using a dot-separated path.

    Path components that are integers index into lists.
    The special syntax "key=val" matches a dict in a list by key.

    Examples:
        "cassette.0.cistron.expression" → obj["cassette"][0]["cistron"]["expression"]
        "cassette.2.spacer" → obj["cassette"][2]["spacer"]
    """
    parts = path.split(".")
    current = obj

    for i, part in enumerate(parts[:-1]):
        if isinstance(current, list):
            idx = int(part)
            current = current[idx]
        elif isinstance(current, dict):
            # Try integer index first (for dict-wrapped list access)
            try:
                idx = int(part)
                # This is a dict containing lists — shouldn't happen normally
                current = current[part]
            except (ValueError, KeyError):
                current = current[part]
        else:
            raise KeyError(f"Cannot traverse into {type(current)} at path component '{part}'")

    # Set the final value
    final_key = parts[-1]
    if isinstance(current, list):
        current[int(final_key)] = value
    elif isinstance(current, dict):
        current[final_key] = value
    else:
        raise KeyError(f"Cannot set '{final_key}' on {type(current)}")
