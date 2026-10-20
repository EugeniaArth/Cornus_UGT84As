# UGT84 reference sequences and candidate-screening framework for *Cornus*

## Abstract

This project assembles a literature- and accession-supported protein reference set for UGT84 sequence-family analysis and prepares a reproducible workflow for subsequent candidate recovery in *Cornus*. The dataset contains 32 amino-acid sequences: 24 named UGT84 A/B/J references, six UGT74 proteins used as Group-L comparators, and two additional family-1 context sequences. The panel distinguishes broad phylogenetic-group membership, UGT84 family/subfamily assignment, NCBI product labels, and experimental activity evidence. It is designed to support candidate comparison, not to assign enzyme function from sequence similarity alone.

The 32-sequence reference dataset and its sequence-level manifest are present. Candidate-screening, annotation-mapping, alignment, and tree-inference scripts are prepared as separate steps, but no Cornus candidate set, alignment, or phylogeny has yet been generated in this project. *C. suecica* analysis is deferred.

## Background and objectives

Plant family-1 UDP-dependent glycosyltransferases (UGTs) are a large and functionally diverse enzyme family. UGT84 proteins include experimentally characterized gallic-acid glucosyltransferases involved in the formation of β-glucogallin, but UGT84 sequence labels and reported activities do not apply uniformly to all family members. A sequence search, conserved motif, NCBI product description, or phylogenetic placement is therefore candidate evidence rather than proof of biochemical activity.

The project has three aims: (1) establish a compact, documented set of literature-named UGT84 proteins; (2) retain relevant neighboring and family-1 context sequences without conflating them with UGT84 references; and (3) prepare a transparent, stepwise method for recovering and evaluating *Cornus* candidates once the species-specific protein and annotation inputs are selected.

## Dataset description

The protein panel contains 24 UGT84 references: 17 UGT84A, three UGT84B, and four UGT84J sequences. The sequence accessions, source headers, amino-acid lengths, source evidence, and available functional evidence are recorded in `Results/ugt84_reference_context_panel_manifest.tsv`. The combined FASTA is `Results/ugt84_reference_context_panel_32.faa`.

| Accession | Taxon | Gene label | UGT84 subfamily | Length (aa) |
|---|---|---|---|---:|
| AHA54051.1 | *Quercus robur* | UGT84A13 | A | 510 |
| ALO19890.1 | *Camellia sinensis* | CsUGT84A22 | A | 487 |
| ALO19891.1 | *Camellia sinensis* | CsUGT84J2 | J | 477 |
| ANN02875.1 | *Punica granatum* | UGT84A23 | A | 506 |
| ANN02877.1 | *Punica granatum* | UGT84A24 | A | 517 |
| BBB21213.1 | *Eucalyptus camaldulensis* | UGT84A25a | A | 499 |
| BBB21212.1 | *Eucalyptus camaldulensis* | UGT84A25b | A | 499 |
| BBB21215.1 | *Eucalyptus camaldulensis* | UGT84A26a | A | 511 |
| BBB21214.1 | *Eucalyptus camaldulensis* | UGT84A26b | A | 511 |
| BBB21216.1 | *Eucalyptus camaldulensis* | UGT84J3 | J | 497 |
| BBB21217.1 | *Eucalyptus camaldulensis* | UGT84J4 | J | 497 |
| BBB21218.1 | *Eucalyptus camaldulensis* | UGT84J5 | J | 469 |
| AEW31187.1 | *Vitis vinifera* | VvGT1 / UGT84A44 | A | 478 |
| AAU09443.1 | *Fragaria × ananassa* | FaGT2 / UGT84A6 | A | 555 |
| QZM06937.1 | *Canarium album* | UGT84A77 | A | 503 |
| XLV46658.1 | *Litsea cubeba* | LcUGT84B1 | B | 489 |
| BAA93039.1 | *Citrus unshiu* | CitLGT / UGT84A5 | A | 511 |
| AAF98390.1 | *Brassica napus* | UGT84A9 | A | 497 |
| CAJ77650.1 | *Brassica napus* | UGT84A10 | A | 472 |
| AII32448.1 | *Populus fremontii × P. angustifolia* | UGT84A17 | A | 502 |
| NP_193284.1 | *Arabidopsis thaliana* | UGT84A3 | A | 479 |
| AEE83611.1 | *Arabidopsis thaliana* | UGT84A4 | A | 475 |
| AAB87119.1 | *Arabidopsis thaliana* | UGT84B1 | B | 456 |
| AAB87106.1 | *Arabidopsis thaliana* | UGT84B2 | B | 438 |

Six UGT74 proteins provide neighboring Group-L comparison. Their labels and accessions are retained exactly as given in the saved protein records; they are not designated UGT84 proteins.

| Accession | Taxon | Gene label | Family / group | Length (aa) |
|---|---|---|---|---:|
| AEC10333.1 | *Arabidopsis thaliana* | UGT74F1 | UGT74 / L | 449 |
| AEC10331.1 | *Arabidopsis thaliana* | UGT74F2 | UGT74 / L | 449 |
| ALO19895.1 | *Camellia sinensis* | CsUGT74B4 | UGT74 / L | 459 |
| ALO19896.1 | *Camellia sinensis* | CsUGT74B5 | UGT74 / L | 471 |
| ALO19897.1 | *Camellia sinensis* | CsUGT74AF1 | UGT74 / L | 468 |
| ALO19898.1 | *Camellia sinensis* | CsUGT74AG1 | UGT74 / L | 455 |

Two additional family-1 proteins are included as context sequences. Medicago UGT71G1 (`AAW56092.1`; 465 aa) belongs to Group E. Stevia UGT76G1 (`AAR06912.1`; 458 aa) belongs to Group H. The grouping scheme in Hoffmann et al. (2023, Table 2) maps UGT71/72/88 to Group E, UGT74/75/84 to Group L, and UGT76 to Group H. Thus the six UGT74 proteins and the UGT84 references share broad Group-L classification; that classification does not by itself establish a direct sister relationship between every pair. Only one Group-E sequence and one Group-H sequence are included; no UGT72 or UGT88 representatives are currently part of the dataset.

The grape accession AEW31187.1 is named VvGT1 / UGT84A44 in Table 1 of the 2024 hydrolyzable-tannin review. The source FASTA header retained in the manifest contains the earlier spelling `VvgGT1`; the display name here follows the cited table. VvGT2 and VvGT3 are excluded because that table does not give them formal UGT84 names.

## Sequence sources and evidence records

The UGT84 reference sequences were assembled from published studies and cited accession records. The saved evidence manifest records sequence-family membership support separately from activity evidence. For example, some entries have direct β-glucogallin assay evidence, some were assayed without reported activity for a tested substrate panel, and others are included as family references without a functional claim. These distinctions are retained at the individual-sequence level; the 24-reference count must not be read as 24 experimentally demonstrated β-glucogallin enzymes.

The core literature sources cover oak UGT84A13 (PMID 24412325), tea UGT84A22/UGT84J2 and Arabidopsis references (PMID 26941235), pomegranate UGT84A23/A24 (PMID 27227328), *Eucalyptus camaldulensis* UGT84A/J proteins (PMID 29775866), the tea UGT catalogue and group classification (PMID 37346124), tannin-pathway reconstitution (PMID 37526750), and the 2024 hydrolyzable-tannin review (PMID 40115765). Additional accession-linked literature identifiers are recorded in the sequence manifest, including PMIDs 22090445, 26884604, 33719298, 34361803, and 40665580. Full citations, stable identifiers, and the names of locally included source records are given in `References/ugt84_expanded/SOURCE_CATALOG.md`.

The project includes the curated source FASTAs and metadata, accession FASTAs for the 32 selected proteins, saved PubMed records, available PMC full-text XML, compressed supplements, and extracted text selections. It does not include the five complete comparison proteomes or their assembly/GFF bundles because those data are not inputs to construction of this 32-sequence dataset. PF00201 profile files are also not included: the precise profile build and redistribution terms remain unverified.

## Candidate-recovery and phylogenetic methods

The reference dataset was constructed from four local source files: the 24 formal UGT84 protein sequences and their evidence manifest, six UGT74 comparator sequences, and the two family-1 context sequences. `Scripts/00_stage_reference_sources.py` copies these inputs into the project’s `Inputs/reference_panel/` area and writes a JSON manifest containing their sizes and SHA-256 hashes. `Scripts/01_build_reference_panel.py` combines the staged inputs into the 32-sequence FASTA and matching sequence manifest. `Scripts/02_validate_reference_panel.py` checks source hashes, sequence and metadata correspondence, expected group counts, sequence lengths, and exclusion of VvGT2/VvGT3.

The prepared candidate workflow consists of separate, non-overwriting steps. A broad PF00201 HMMER search recovers UDP-glycosyltransferase-family candidates from a supplied annotated protein FASTA. Its raw hit table is converted to a candidate FASTA; BLASTP then compares those sequences against the 24 formal UGT84 references. Candidate selection uses explicit identity, query-coverage, and E-value settings, followed by protein-to-GFF3 mapping, locus-level representative selection, and assembly of a combined reference-plus-Cornus input. PF00201 is a broad UDP-glycosyltransferase-family profile, not a UGT84 classifier. BLAST similarity likewise supplies candidate-generation evidence only.

MAFFT alignment and IQ-TREE inference are prepared as the final two stages and are intentionally not run until the Cornus candidate set is finalized. The intended tree will be an exploratory sequence-family analysis; any rooting choice and support labels must be documented separately. No phylogenetic placement or enzymatic function is assigned in the present dataset.

## Reproducibility and current status

From within the project directory, the reference-source staging, panel construction, and validation commands are:

```bash
python3 Scripts/00_stage_reference_sources.py
python3 Scripts/01_build_reference_panel.py
python3 Scripts/02_validate_reference_panel.py
```

The scripts refuse to overwrite existing outputs. Staged source inputs are copied byte-for-byte from the included curated sequence records and verified by SHA-256. The resulting panel contains 32 unique accessions: 24 UGT84, six UGT74 Group-L comparators, and two family-1 context sequences. Validation expects 30 Group-L records, one Group-E record, and one Group-H record.

The panel counts and per-sequence descriptions can be recalculated from the included TSV without a separate report generator:

```python
import csv
from collections import Counter
from pathlib import Path

rows = list(csv.DictReader(
    Path("Results/ugt84_reference_context_panel_manifest.tsv").open(),
    delimiter="\t"))
print(Counter(row["panel_role"] for row in rows))
print(Counter(row["phylogenetic_group"] for row in rows))
```

The reference inputs and output dataset are present and locally traceable. The candidate-recovery, GFF mapping, alignment, and tree scripts have not been applied to Cornus data here. In particular, no *C. florida* or *C. suecica* candidate is included in the 32-sequence dataset, and no final tree has been produced.

## Limitations and interpretation

UGT family and group labels are classification evidence, not a guarantee of biochemical specificity. NCBI product names, sequence similarity, conserved motifs, and future tree placement cannot establish substrate activity. Functional statements must remain limited to the individual proteins and assay conditions reported in the cited studies.

The Group-E representation is limited to UGT71G1; the panel does not contain UGT72 or UGT88 examples. UGT76G1 is a Group-H context sequence, not a Group-E member. The PF00201 profile build, complete citation metadata for several accession-linked PMIDs, and source-specific redistribution terms require further verification before a publication or external data release. *C. suecica* analysis remains deferred.
