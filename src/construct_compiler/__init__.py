"""
construct_compiler — A genetic construct design compiler.

Compiles declarative YAML specs into annotated DNA sequences,
assembly plans, and cost estimates.

Usage:
    from construct_compiler import compile_construct
    graph, plan = compile_construct("my_construct.yaml")
"""

from .passes.pipeline import compile_construct
from .frontend.parser import parse_spec
from .core.graph import ConstructGraph
from .backends.genbank import export_genbank

__version__ = "0.1.0"
__all__ = ["compile_construct", "parse_spec", "ConstructGraph", "export_genbank"]
