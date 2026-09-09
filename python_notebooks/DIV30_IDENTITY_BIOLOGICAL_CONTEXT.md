# Biological context for the fixed DIV30 identity analysis

These sources supply interpretation, not new target labels or evidence that a
DIV30 cell becomes a particular DIV90 cell. The frozen four-gene PV module is
MEF2C/KCNC1/KCNC2/TAC1. Its published biological motivation is separate from
statistical independence: direct and ranked expression of these genes already
contribute to the frozen DIV30 consensus. The same distinction applies to later
reference labels and shared mapping features.

## GRIA2

GRIA2 encodes the GluA2 AMPA receptor subunit. Expression is compatible with
neuronal glutamatergic input and receptor development, but does not specify
PV-producing fate. In lineage-traced **mouse hippocampal** interneurons, Matta
and colleagues found that MGE-derived synapses predominantly used receptors
lacking GluA2, whereas CGE-derived synapses used GluA2-containing receptors.
This is an important counterexample to equating GRIA2-high with mature PV
identity. It cannot be used to relabel human DIV30 GRIA2-high cells as CGE:
species, region, developmental stage, RNA, and functional synaptic receptor
composition differ. [Matta et al., Nature Neuroscience 2013](https://www.nature.com/articles/nn.3459).

Human GRIA2 variant experiments support a role for GluA2 in neuronal receptor
function and neurodevelopment, rather than a selective PV-lineage marker.
[Salpietro et al., Nature Communications 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6626132/).

## OPCML

OPCML is an IgLON cell-adhesion molecule. Embryonic mouse expression studies
show spatially and temporally regulated IgLON expression, including different
OPCML promoter isoforms. This supports interpretation as a developmental
adhesion/recognition axis, without demonstrating human PV specificity.
RNA abundance does not resolve isoform-specific extracellular abundance.
[Jagomäe et al., Frontiers in Neuroscience 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8268470/),
[Spatiotemporal expression of IgLON family members, Scientific Reports 2021](https://www.nature.com/articles/s41598-021-97768-5).

## NOTCH1

Low NOTCH1 is biologically compatible with departure from a progenitor program.
Perturbation experiments in chick neural tube and mouse neural progenitors
showed that Prox1 suppresses Notch1 and permits cell-cycle exit and neuronal
differentiation. This supports the direction of the fixed exclusion clause,
but low receptor RNA is neither a direct Notch-activity measurement nor a
PV-fate assay. [Kaltezioti et al., PLOS Biology 2010](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1000565).

## SLC6A1

SLC6A1 encodes GAT-1, a GABA transporter. Studies in neurons, astrocytes and
patient-derived cells link transporter localization and function to GABA
uptake. SLC6A1 RNA is therefore relevant to GABA-handling development, but is
not exclusive to inhibitory neurons, does not establish transport activity,
and cannot by itself explain PVALB protein production.
[Mermer et al., Brain 2021](https://academic.oup.com/brain/article/144/8/2237/6324634),
[Common molecular mechanisms in astrocytes and neurons, 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8418336/).

SLC6A1 was not in the explicit four-gene PV module but **was already included**
in the frozen mapping features and sparse model (coefficient +0.010314168).
Thus, measuring it as a separate axis does not make its association with the
old consensus fully independent. The new omit-SLC6A1 mapping is a diagnostic
sensitivity with fixed labels, not a replacement target.

## External-protocol premise

The comparison cannot start from an assumption that other organoid protocols
never generate PV protein. Walsh and colleagues report PVALB-positive,
fast-spiking human interneurons in forebrain assembloids. Their primary study
is also represented by a locally available transcriptomic object. The local
object's age, cell selection and sequencing sensitivity must be respected
when interpreting its resemblance to the DIV30 population.
[Walsh et al., Neuron 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12447774/).

## Prospective sorting remains conditional

The identity report determines whether a PV-directed surface-reagent review
is justified. No antibody purchase is recommended from RNA separation alone.
The exact numerical thresholds are in log-normalized RNA units and are not
fluorescence cutoffs. Protein accessibility, abundance on dissociated viable
human cells, and any perturbation caused by antibody binding would need
direct validation if that next stage is justified by the biological results.

Sources checked 2026-09-09. All mechanistic extrapolations above are explicitly
limited to the experimental species and setting reported by the source.
