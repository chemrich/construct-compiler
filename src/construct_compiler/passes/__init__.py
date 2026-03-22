"""
Compiler passes — each transforms the IR graph, moving parts from
ABSTRACT → RESOLVED → CONCRETE.

Passes are composable functions: pass(graph) -> graph
"""
from .part_resolution import resolve_parts
from .reverse_translation import reverse_translate
from .constraint_resolution import resolve_constraints
from .pipeline import compile_construct
