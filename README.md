# construct-compiler

A genetic construct design compiler for *in silico* DNA design. Write a declarative YAML spec describing what you want to express, and the compiler produces annotated DNA sequences, assembly plans, vendor cost estimates, and GenBank files.

Built for researchers who work with recombinant protein expression in *E. coli* and want to automate the repetitive parts of construct design — choosing codons, picking RBSs, avoiding restriction sites, and figuring out whether to synthesize or clone.

## What it does

```
YAML spec  →  parse  →  IR graph  →  optimization passes  →  assembly plan + GenBank output
```

The compiler takes a high-level description of your construct (which protein, which tags, which host, which assembly method) and lowers it through a series of passes into a concrete, annotated DNA sequence with a costed build plan.

**Passes:**

1. **Part resolution** — fetches protein sequences from UniProt and FPbase, looks up regulatory elements (promoters, RBSs, terminators) and fusion tags from a curated local database.
2. **Reverse translation** — converts amino acid sequences to DNA using host-preferred codons.
3. **Constraint resolution** — uses [DNA Chisel](https://github.com/Edinburgh-Genome-Foundry/DnaChisel) to codon-optimize while simultaneously enforcing constraints: no internal restriction sites (for Golden Gate compatibility), GC content windows, no long homopolymer runs, synthesis vendor compatibility.
4. **Assembly planning** — evaluates strategies (Twist clonal, 2-part Golden Gate, 3-part Golden Gate, IDT gBlock) and recommends the cheapest option, accounting for researcher time at a configurable hourly rate.

## Quick start

```bash
pip install -e ".[dev,web]"
```

There are three ways to use construct-compiler: the **CLI**, the **Python API**, and the **web UI**.

### CLI

```bash
# Compile a spec and get cost comparison + GenBank file
construct-compiler compile my_construct.yaml -o output/

# Just show cost comparison
construct-compiler compile my_construct.yaml --cost-only -q

# Validate a spec without compiling
construct-compiler validate my_construct.yaml

# List available parts
construct-compiler parts --list tags
construct-compiler parts --list promoters
```

### Web UI

Start the server and open `http://localhost:8421`:

```bash
python -m construct_compiler.server
```

The web UI provides a visual construct builder with a circular plasmid map, interactive cost comparison charts, parts list, and GenBank export. No YAML writing needed — configure your construct through the form and click Compile.

### Python API

Write a construct spec:

```yaml
# my_construct.yaml
construct:
  name: "His-TEV-mEGFP"
  host: e_coli_bl21

  backbone:
    resistance: kanamycin
    ori: pBR322

  cassette:
    - promoter: T7lac
    - cistron:
        label: target
        expression: high
        chain:
          - tag: 6xHis
          - cleavage_site: TEV
          - gene:
              id: mEGFP
              source: fpbase
    - terminator: rrnB_T1

  constraints:
    assembly: golden_gate
    enzyme: BsaI
```

Compile it:

```python
from construct_compiler import compile_construct, export_genbank

graph, plan = compile_construct("my_construct.yaml")
export_genbank(graph, "my_construct.gb")
print(plan.summary())
```

## Example output

For a His-TEV-MBP-mEGFP fusion with a polycistronic mScarlet reporter (see `examples/`):

```
Assembly Plan: HisTEV-MBP-mEGFP_with_mScarlet_reporter
Insert length: 3045 bp

Strategy: Twist Clonal Gene *** RECOMMENDED ***
  Synthesis cost: $274.05
  Reagents (with overhead): $48.00
  Researcher time: 1.5 hrs ($225.00)
  TOTAL COST: $547.05
  Wall clock: 13-20 business days

Strategy: Synthesis + 2-Part Golden Gate
  Synthesis cost: $213.15
  Researcher time: 3.5 hrs ($525.00)
  TOTAL COST: $849.56

Strategy: IDT gBlock + 2-Part Golden Gate
  Synthesis cost: $106.64
  TOTAL COST: $743.05
  Wall clock: 6-9 business days (fastest)
```

## Construct spec reference

### Backbone

```yaml
backbone:
  source: addgene          # or "local"
  addgene_id: 26094        # optional
  resistance: kanamycin    # amp, kan, chlor, spec, etc.
  ori: pBR322              # pBR322, p15A, ColA, CloDF13, SC101
```

### Promoters

Built-in: `T7`, `T7lac`, `tac`, `araBAD`, `lacUV5`, `J23100` (constitutive), `J23106` (constitutive, medium).

### RBS / expression levels

Use `expression: high|medium|low|very_low` in a cistron block. The compiler maps these to context-insensitive bicistronic design (BCD) elements from Mutalik et al. 2013, or you can specify a part name directly:

```yaml
cistron:
  rbs: BBa_B0034       # explicit part name
  # or
  expression: medium    # auto-selects BCD12
```

| Level | Default part | Relative strength |
|-------|-------------|-------------------|
| high | BCD2 | 1.0 |
| medium | BCD12 | 0.2 |
| low | BCD22 | 0.05 |
| very_low | BBa_B0033 | 0.01 |

### Fusion tags and cleavage sites

Tags and cleavage sites can be specified in a `chain` (ordered) or as `n_tag`/`c_tag` shorthand:

```yaml
# Chain syntax (ordered, explicit)
chain:
  - tag: 6xHis
  - cleavage_site: TEV
  - solubility_tag: MBP
  - linker: {type: GS_flexible, repeats: 3}
  - gene: {id: mEGFP, source: fpbase}

# Shorthand syntax
n_tag: [6xHis, TEV]
gene: {id: mEGFP, source: fpbase}
c_tag: Strep-II
```

**Purification tags:** `6xHis`, `8xHis`, `Strep-II`, `Twin-Strep`, `FLAG`, `HA`

**Solubility tags:** `MBP`, `GST`, `SUMO`, `Trx`

**Cleavage sites:** `TEV`, `3C`, `Factor_Xa`, `Thrombin`, `Enterokinase`

**Linkers:** `GS_flexible` (GGGGS)n, `rigid_EAAAK` (EAAAK)n, `short_GS`

### Polycistronic designs

Add multiple `cistron` blocks under a single promoter. Use `fused: false` for independently translated reporters:

```yaml
cassette:
  - promoter: T7lac
  - cistron:
      label: target
      expression: high
      chain:
        - tag: 6xHis
        - gene: {id: P12345, source: uniprot}
  - spacer: 30
  - cistron:
      label: reporter
      expression: low
      fused: false
      gene: {id: mScarlet-I, source: fpbase}
  - terminator: rrnB_T1
```

### Constraints

```yaml
constraints:
  assembly: golden_gate       # only golden_gate currently
  enzyme: BsaI                # BsaI, BpiI, BbsI
  codon_optimization: local   # local (DNA Chisel) or vendor (Twist/IDT)
  gc_window: [0.35, 0.65]    # min/max GC content
  max_homopolymer: 6          # max consecutive identical bases
```

### Synthesis preferences

```yaml
synthesis:
  vendor: twist
  budget_per_gene: 400
  max_turnaround_days: 21
```

## Cost model

The assembly planner compares strategies using a fully-loaded cost model:

| Parameter | Default |
|-----------|---------|
| Researcher time | $150/hr ($300k/yr) |
| Overhead multiplier | 1.5x on consumables |
| Twist gene synthesis | $0.07/bp |
| Twist clonal gene | $0.09/bp |
| IDT gBlock | $0.08/bp |
| Golden Gate reagents | ~$40/reaction |
| Assembly success rate | 90% (2-part), 80% (3-part) |

All parameters are configurable via `CostParams`:

```python
from construct_compiler.passes.assembly_planning import CostParams

params = CostParams(
    researcher_hourly_rate=100.0,    # adjust for your lab
    twist_gene_per_bp=0.06,          # volume discount
    overhead_multiplier=1.65,        # your institution's rate
)
graph, plan = compile_construct("spec.yaml", cost_params=params)
```

## Vendor integration

Twist and IDT vendor plugins support screening, codon optimization, and ordering via their APIs. Set credentials as environment variables:

```bash
# Twist Bioscience TAPI
export TWIST_API_KEY=your_key
export TWIST_API_SECRET=your_secret

# IDT SciTools Plus
export IDT_CLIENT_ID=your_id
export IDT_CLIENT_SECRET=your_secret
```

Without credentials, the plugins run in mock mode with heuristic feasibility checks and pricing estimates.

## Architecture

The compiler follows a frontend → IR → passes → backend architecture:

```
┌─────────────┐     ┌──────────────────┐     ┌─────────────┐
│  YAML spec  │────▶│   IR: Typed DAG   │────▶│  GenBank    │
│  (frontend) │     │  of genetic parts │     │  (backend)  │
└─────────────┘     └────────┬─────────┘     └─────────────┘
                             │
                    ┌────────┴─────────┐
                    │  Compiler passes  │
                    │  1. Part resolve  │
                    │  2. Rev translate │
                    │  3. Constraints   │
                    │  4. Assembly plan │
                    └──────────────────┘
```

The IR is a directed graph where nodes are genetic parts with typed ports. Port types (TRANSCRIPTION, TRANSLATION_INIT, PEPTIDE_CHAIN, etc.) enforce biological validity at composition time — putting a terminator after a promoter with no coding sequence in between is a type error.

## Project structure

```
construct_compiler/
├── src/construct_compiler/
│   ├── __main__.py     # CLI (Click)
│   ├── server.py       # FastAPI server for web UI
│   ├── core/           # IR: types, parts, graph
│   ├── frontend/       # YAML parser (spec → IR)
│   ├── passes/         # Compiler passes + pipeline
│   ├── backends/       # GenBank export (+ future SBOL3, protocols)
│   ├── vendors/        # Twist, IDT plugin stubs
│   ├── data/           # Curated parts DB, codon tables, overhang sets
│   └── plugins/        # Plugin system (future)
├── frontend/           # React web UI (single-page app)
│   └── index.html
├── examples/           # Example YAML specs and runner scripts
├── pyproject.toml
└── README.md
```

## Roadmap

- [ ] Variant library fan-out (compile N constructs from parameterized specs)
- [ ] Live Twist/IDT API integration (screening + vendor codon optimization)
- [ ] Protocol generation backend (human-readable assembly instructions)
- [ ] Primer design backend (primer3-py for Golden Gate primers with overhangs)
- [ ] SBOL3 export via pySBOL3
- [ ] Addgene backbone fetching (auto-download backbone sequences)
- [ ] Salis RBS Calculator integration for computed translation initiation rates
- [ ] Verification target (expected digest fragments, colony PCR bands, Sanger references)

## License

MIT
