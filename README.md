# UGT84 evolution in dogwoods (*Cornus*)

**Project title:** Evolution of the UGT84 glycosyltransferase subfamily in *Cornus* (Cornaceae): gene family structure and predicted relevance to galloylated metabolite biosynthesis  
**Short title:** UGT84 evolution in dogwoods (*Cornus*)

## Project aim

Characterize the UGT84 glycosyltransferase subfamily across available *Cornus* genomes. UGT84 enzymes include proteins that form glucose esters of phenolic acids; in other plants, β-glucogallin formation is an entry step in hydrolysable tannin biosynthesis. The project will reconstruct UGT84 phylogeny and duplication history, assess selection on copies, and predict whether *Cornus* proteins resemble characterized gallic-acid-active enzymes. Functional conclusions will remain predictions unless supported by enzyme assays.

## Hypotheses

1. *Cornus* species share a conserved set of orthologous UGT84A copies, alongside some lineage-specific duplications.
2. Test whether *Cornus* UGT84A proteins cluster with characterized β-glucogallin-forming UGT84A enzymes, while evaluating UGT84B and UGT84J as comparators. Do not assume UGT84B is uniformly inactive: *Litsea cubeba* UGT84B1 has been reported to form β-glucogallin.
3. Duplicated copies show relaxed or divergent selection relative to conserved core copies.

## Data currently present

The full project workspace contains six NCBI assembly archives and extracted genome files. This lightweight GitHub project copy intentionally excludes those large assemblies and the curated reference FASTA files. The small reference manifest and literature report are included under `references/`.

| Species | Assembly | Genome FASTA in data workspace | Separate annotation files observed |
|---|---|---|---|
| *Cornus controversa* | GCA_026545225.1 | Yes | No separate annotation bundle observed |
| *Cornus florida* | GCF_030987335.1 | Yes | GFF, GTF, CDS, protein, and RNA files |
| *Cornus kousa* | GCA_057556605.1 | Yes | No separate annotation bundle observed |
| *Cornus sanguinea* | GCA_984573785.1 | Yes | No separate annotation bundle observed |
| *Cornus sericea* | GCA_986885575.1 | Yes | No separate annotation bundle observed |
| *Cornus suecica* | GCA_963922005.1 | Yes | GFF3 and GTF added from Ensembl Darwin Tree of Life project |

A preliminary screening of the first GenBank flatfile record found gene/CDS/transcript features for *C. florida* but none for the other five assemblies. *C. suecica* now also has separate GFF3/GTF annotation files. A full annotation audit and validation of those files remain pending; do not infer annotation absence solely from the flatfiles. GenBank flatfiles are present for all six assemblies.

Some assembly reports list 2026 release dates (*C. kousa*, *C. sanguinea*, and *C. sericea*). These dates require verification against the source records before use.

## Workflow and methods

### 1. Data inventory and assembly quality control

- Inventory assemblies, annotation files, archives, and their provenance.
- Calculate sequence counts, total FASTA length, N50, GC percentage over canonical A/C/G/T bases, and ambiguous-base counts from the FASTA files; compare FASTA sequence lengths with the NCBI sequence reports.
- Check annotation status in both separate annotation bundles and GenBank flatfiles.
- Run BUSCO with a consistent plant lineage and keep genome-mode and protein-mode results distinct. BUSCO is installed in the active environment, but no BUSCO lineage download or external data access has been authorized; BUSCO results are therefore pending.

### 2. Trusted reference set

Start from the selected high-confidence UGT84A sequence set in `references/curated_ugt84/manifest.tsv` and the evidence review in `references/ugt84_beta_glucogallin_evidence.md`. Include UGT84B1 as a known β-glucogallin-forming comparator and use verified UGT84J and UGT74/UGT75 representatives as appropriate. Record accessions, sources, reported activities, and evidence quality. Gallic-acid activity alone is not evidence of a tannin-specific role.

### 3. Candidate identification in annotated *C. florida*

Search the annotated proteins with HMMER using PF00201, validate domains with InterProScan, assess PSPG motifs and sequence integrity, and classify candidates using phylogeny and curated references. PF00201 is a broad UDP-glycosyltransferase-domain screen, not a UGT84-specific classifier.

### 4. Candidate recovery in every genome

Use complementary methods:

- Protein-to-genome mapping with miniprot to identify loci and produce splice-aware candidate models.
- Nucleotide homology searches using trusted CDS/transcript sequences, followed by locus expansion, splice-aware exon resolution, and manual review.
- Compare the methods rather than blindly combining their hits. Inspect ORF integrity, splice boundaries, motif preservation, alignment coverage, and genomic context.
- Retain partial models and possible pseudogenes with explicit flags. Report intact-gene counts separately from broader candidate counts.

Do not rely on de novo prediction alone: missed or incorrectly modeled orthologs are a known risk. De novo predictions may be used as supporting evidence where useful.

### 5. Recovery benchmark in *C. florida*

Run both homology-recovery approaches on *C. florida* and compare models with its annotation. Treat annotation as a benchmark reference, not absolute ground truth; annotation errors or omissions may affect apparent precision and recall.

### 6. Alignment, phylogeny, and orthology

Align curated proteins with MAFFT, inspect and trim conservatively, and verify that motif interpretation and candidate placement are not artifacts of trimming. Select models and infer trees with IQ-TREE using UFBoot and SH-aLRT. Use AU tests for relevant topology constraints. Root and interpret the tree with UGT74/UGT75 outgroups, checking sensitivity to outgroup choice and alignment quality. Assign orthology groups only after reviewing support and model quality.

### 7. Duplication and selection

In *C. florida*, classify tandem and dispersed copies using genomic coordinates; use MCScanX/collinearity evidence for segmental duplication where the assembly and annotation support it. Build codon alignments from curated coding sequences before selection analyses. Use RELAX or codeml only on defensible ortholog/paralog sets and report limitations from small sample sizes, uncertain models, or low power.

### 8. Functional predictions and reporting

Compare candidate residues with characterized active and inactive references; present results as predictions, not demonstrated biochemical function. Planned outputs include a candidate copy-number table with quality flags, an annotated phylogeny, the *C. florida* recovery benchmark, duplication and selection results, candidate enzyme predictions, and a limitations section.

## Step 1 files and reproducibility status

- `scripts/step1_inventory_qc.py` — preliminary local inventory/QC script; it is **not yet validated**.
- `references/curated_ugt84/manifest.tsv` — accessions, selection rationale, lengths, and translation checks for the selected reference set. The protein/CDS FASTA files are deliberately not included in this lightweight copy.
- `references/ugt84_beta_glucogallin_evidence.md` and `.html` — literature review and accession cross-checks.

The assembly metadata and QC output are not included in this package. Step 1 validation is pending: the current workspace metadata/QC files need correction before the script can be relied on. The script requires the local genome FASTAs, sequence reports, and a corrected `genome_metadata.csv`; those inputs are not bundled here. Do not treat the earlier `assembly_qc.tsv` as a validated result.

## Limitations and interpretation

- There are only six assemblies, and assembly level and quality differ among them.
- Missing separate annotation bundles do not by themselves prove that GenBank flatfiles lack annotations; full feature inspection is needed. The *C. suecica* GFF3/GTF files also require provenance/version and assembly-coordinate validation before use.
- Homology-derived models may be partial or biased, especially in unannotated assemblies; copy counts should be presented as lower bounds with confidence flags.
- No enzyme assays in *Cornus* are included, so substrate specificity and biological role remain predictions.
- The set of characterized reference enzymes is limited, and gallic-acid activity is not by itself diagnostic of a tannin-specific role.
- Selection tests and duplication classifications depend on accurate models, alignments, and genomic context; weakly supported conclusions should be labeled accordingly.
