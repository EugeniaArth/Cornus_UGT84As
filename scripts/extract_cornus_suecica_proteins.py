#!/usr/bin/env python3
"""Validate local Suecica gene models and translate protein-coding GFF3 CDSs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import unquote

PROJECT = Path(__file__).resolve().parents[1]
INPUTS = PROJECT / "Inputs" / "cornus_suecica"
ASSEMBLY = "GCA_963922005.1"
GENOME = INPUTS / "GCA_963922005.1_drCorSuec2.1_genomic.fna"
GFF3 = INPUTS / "cornus_suecica_gca963922005v1.gff3"
GTF = INPUTS / "cornus_suecica_gca963922005v1.gtf"
SEQ_REPORT = INPUTS / "sequence_report.jsonl"
DEFAULT_OUTDIR = PROJECT / "Results" / "cornus_suecica" / "proteome_validated_v3"

CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}
STOP_CODONS = {"TAA", "TAG", "TGA"}
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_gff_attributes(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in text.strip().strip(";").split(";"):
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[unquote(key)] = unquote(value)
    return attrs


def parse_gtf_attributes(text: str) -> dict[str, str]:
    return {key: value for key, value in re.findall(r'([^;\s]+)\s+"([^"]*)"', text)}


def norm_id(value: str, prefix: str) -> str:
    return value[len(prefix):] if value.startswith(prefix) else value


def fasta_dict(path: Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    current: str | None = None
    with path.open(encoding="ascii", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if not current or current in sequences:
                    raise ValueError(f"Empty or duplicate genome FASTA ID: {current!r}")
                sequences[current] = []
            else:
                if current is None:
                    raise ValueError("Genome sequence occurs before first FASTA header")
                sequences[current].append(line.upper())
    return {key: "".join(parts) for key, parts in sequences.items()}


def merged_intervals(rows: list[dict[str, object]]) -> list[tuple[str, int, int, str]]:
    grouped: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["seqid"]), str(row["strand"]))].append((int(row["start"]), int(row["end"])))
    result: list[tuple[str, int, int, str]] = []
    for (seqid, strand), intervals in sorted(grouped.items()):
        current_start = current_end = None
        for start, end in sorted(intervals):
            if current_start is None:
                current_start, current_end = start, end
            elif start <= int(current_end) + 1:
                current_end = max(int(current_end), end)
            else:
                result.append((seqid, int(current_start), int(current_end), strand))
                current_start, current_end = start, end
        if current_start is not None:
            result.append((seqid, int(current_start), int(current_end), strand))
    return result


def feature_coords(rows: list[dict[str, object]]) -> list[tuple[str, int, int, str]]:
    return sorted((str(r["seqid"]), int(r["start"]), int(r["end"]), str(r["strand"])) for r in rows)


def coding_order(rows: list[dict[str, object]], strand: str) -> list[dict[str, object]]:
    return sorted(rows, key=lambda row: (int(row["start"]), int(row["end"])), reverse=(strand == "-"))


def extract_segments(rows: list[dict[str, object]], strand: str, genome: dict[str, str], aliases: dict[str, str]) -> str:
    chunks: list[str] = []
    for row in coding_order(rows, strand):
        seqid = str(row["seqid"])
        accession = aliases[seqid]
        piece = genome[accession][int(row["start"]) - 1:int(row["end"])]
        chunks.append(piece.translate(COMPLEMENT)[::-1] if strand == "-" else piece)
    return "".join(chunks)


def translate(dna: str) -> tuple[str, int]:
    amino: list[str] = []
    ambiguous = 0
    for start in range(0, len(dna) - 2, 3):
        codon = dna[start:start + 3]
        residue = CODON_TABLE.get(codon, "X")
        if residue == "X":
            ambiguous += 1
        amino.append(residue)
    return "".join(amino), ambiguous


def write_fasta_record(handle, identifier: str, sequence: str) -> None:
    handle.write(f">{identifier}\n")
    for start in range(0, len(sequence), 60):
        handle.write(sequence[start:start + 60] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR,
                        help="New, non-existing output directory")
    args = parser.parse_args()
    for path in (GENOME, GFF3, GTF, SEQ_REPORT):
        if not path.is_file():
            raise FileNotFoundError(path)
    outdir = args.outdir.resolve()
    if outdir.exists():
        raise FileExistsError(f"Refusing to overwrite extraction outputs: {outdir}")

    report_rows = [json.loads(line) for line in SEQ_REPORT.read_text(encoding="utf-8").splitlines() if line.strip()]
    aliases: dict[str, str] = {}
    report_by_accession: dict[str, dict[str, object]] = {}
    for row in report_rows:
        accession = str(row.get("genbankAccession", ""))
        if not accession:
            continue
        report_by_accession[accession] = row
        alias_values = [str(row.get("sequenceName", "")), accession]
        if row.get("role") == "assembled-molecule":
            alias_values.append(str(row.get("chrName", "")))
        for alias in alias_values:
            if not alias:
                continue
            previous = aliases.get(alias)
            if previous and previous != accession:
                raise ValueError(f"Ambiguous sequence-report alias {alias}: {previous} vs {accession}")
            aliases[alias] = accession

    genome = fasta_dict(GENOME)
    if set(genome) != set(report_by_accession):
        raise ValueError(
            "Genome FASTA and sequence report record IDs differ: "
            f"FASTA-only={sorted(set(genome)-set(report_by_accession))[:10]}, "
            f"report-only={sorted(set(report_by_accession)-set(genome))[:10]}"
        )
    length_mismatches = []
    for accession, row in report_by_accession.items():
        expected = int(row["length"])
        if len(genome[accession]) != expected:
            length_mismatches.append((accession, len(genome[accession]), expected))
    if length_mismatches:
        raise ValueError(f"Genome/sequence-report length mismatches: {length_mismatches[:10]}")

    feature_counts_gff: Counter[str] = Counter()
    all_gff_ids: set[str] = set()
    gff_id_kinds: dict[str, str] = {}
    duplicate_cds_id_rows = 0
    all_gff_parents: list[tuple[str, str]] = []
    gene_meta: dict[str, dict[str, object]] = {}
    transcript_meta: dict[str, dict[str, object]] = {}
    gff_by_transcript: dict[str, dict[str, list[dict[str, object]]]] = defaultdict(lambda: defaultdict(list))
    gff_seqids: set[str] = set()
    gff_bound_errors: list[str] = []

    with GFF3.open(encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) != 9:
                raise ValueError(f"Malformed GFF3 line {line_no}: expected 9 columns")
            seqid, source, kind, start_s, end_s, score, strand, phase, attr_text = columns
            try:
                start, end = int(start_s), int(end_s)
            except ValueError as exc:
                raise ValueError(f"Non-integer GFF3 coordinates at line {line_no}") from exc
            feature_counts_gff[kind] += 1
            gff_seqids.add(seqid)
            if start < 1 or end < start:
                gff_bound_errors.append(f"line {line_no}: {seqid}:{start}-{end}")
            if seqid not in aliases:
                raise ValueError(f"GFF3 sequence ID absent from sequence report: {seqid}")
            accession = aliases[seqid]
            if end > len(genome[accession]):
                gff_bound_errors.append(f"line {line_no}: {seqid}:{start}-{end} exceeds {len(genome[accession])}")
            attrs = parse_gff_attributes(attr_text)
            feature_id = attrs.get("ID", "")
            if feature_id:
                if feature_id in all_gff_ids:
                    if kind != "CDS" or gff_id_kinds.get(feature_id) != "CDS":
                        raise ValueError(f"Duplicate non-discontinuous GFF3 feature ID {feature_id} at line {line_no}")
                    duplicate_cds_id_rows += 1
                else:
                    all_gff_ids.add(feature_id)
                    gff_id_kinds[feature_id] = kind
            parents = [value for value in attrs.get("Parent", "").split(",") if value]
            all_gff_parents.extend((feature_id, parent) for parent in parents)
            if kind == "gene":
                gene_id = norm_id(attrs.get("gene_id", feature_id), "gene:")
                if gene_id:
                    gene_meta[gene_id] = {
                        "name": attrs.get("Name", ""),
                        "biotype": attrs.get("biotype", ""),
                        "description": attrs.get("description", attrs.get("product", "")),
                        "seqid": seqid, "start": start, "end": end, "strand": strand,
                    }
            if kind in {"mRNA", "transcript"}:
                transcript_id = norm_id(attrs.get("transcript_id", feature_id), "transcript:")
                if transcript_id:
                    if transcript_id in transcript_meta:
                        raise ValueError(f"Duplicate GFF3 transcript identifier {transcript_id}")
                    parent_genes = [norm_id(p, "gene:") for p in parents]
                    transcript_meta[transcript_id] = {
                        "protein_coding": attrs.get("biotype", "") == "protein_coding",
                        "gene_id": attrs.get("gene_id", parent_genes[0] if len(parent_genes) == 1 else ""),
                        "seqid": seqid, "start": start, "end": end, "strand": strand,
                        "tag": attrs.get("tag", ""), "attrs": attrs,
                    }
            if kind in {"CDS", "exon", "five_prime_UTR", "three_prime_UTR"}:
                row = {"seqid": seqid, "start": start, "end": end, "strand": strand,
                       "phase": phase, "attrs": attrs}
                for parent in parents:
                    transcript_id = norm_id(parent, "transcript:")
                    gff_by_transcript[transcript_id][kind].append(row)

    if gff_bound_errors:
        raise ValueError(f"GFF3 coordinate errors: {gff_bound_errors[:10]}")
    for child, parent in all_gff_parents:
        if parent not in all_gff_ids:
            raise ValueError(f"GFF3 Parent ID absent from feature IDs: {parent} (child {child})")

    gtf_counts: Counter[str] = Counter()
    gtf_meta: dict[str, dict[str, object]] = {}
    gtf_by_transcript: dict[str, dict[str, list[dict[str, object]]]] = defaultdict(lambda: defaultdict(list))
    gtf_seqids: set[str] = set()
    gtf_bound_errors: list[str] = []
    with GTF.open(encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) != 9:
                raise ValueError(f"Malformed GTF line {line_no}: expected 9 columns")
            seqid, source, kind, start_s, end_s, score, strand, phase, attr_text = columns
            try:
                start, end = int(start_s), int(end_s)
            except ValueError as exc:
                raise ValueError(f"Non-integer GTF coordinates at line {line_no}") from exc
            gtf_counts[kind] += 1
            gtf_seqids.add(seqid)
            if start < 1 or end < start:
                gtf_bound_errors.append(f"line {line_no}: {seqid}:{start}-{end}")
            if seqid not in aliases:
                raise ValueError(f"GTF sequence ID absent from sequence report: {seqid}")
            accession = aliases[seqid]
            if end > len(genome[accession]):
                gtf_bound_errors.append(f"line {line_no}: {seqid}:{start}-{end} exceeds {len(genome[accession])}")
            attrs = parse_gtf_attributes(attr_text)
            transcript_id = attrs.get("transcript_id", "")
            if kind == "transcript" and transcript_id:
                if transcript_id in gtf_meta:
                    raise ValueError(f"Duplicate GTF transcript row {transcript_id}")
                gtf_meta[transcript_id] = {
                    "gene_id": attrs.get("gene_id", ""),
                    "biotype": attrs.get("transcript_biotype", ""),
                    "seqid": seqid, "start": start, "end": end, "strand": strand,
                }
            if transcript_id and kind in {"CDS", "exon", "five_prime_utr", "three_prime_utr", "start_codon", "stop_codon"}:
                row = {"seqid": seqid, "start": start, "end": end, "strand": strand,
                       "phase": phase, "attrs": attrs}
                gtf_by_transcript[transcript_id][kind].append(row)
    if gtf_bound_errors:
        raise ValueError(f"GTF coordinate errors: {gtf_bound_errors[:10]}")
    if gff_seqids != gtf_seqids:
        raise ValueError(f"GFF3/GTF sequence ID sets differ: GFF3-only={sorted(gff_seqids-gtf_seqids)}, GTF-only={sorted(gtf_seqids-gff_seqids)}")

    coding_ids = {tid for tid, meta in transcript_meta.items() if meta["protein_coding"]}
    gtf_coding_ids = {tid for tid, meta in gtf_meta.items() if meta["biotype"] == "protein_coding"}
    if coding_ids != gtf_coding_ids:
        raise ValueError(
            f"GFF3/GTF protein-coding transcript sets differ: GFF3={len(coding_ids)}, GTF={len(gtf_coding_ids)}, "
            f"GFF3-only={sorted(coding_ids-gtf_coding_ids)[:10]}, GTF-only={sorted(gtf_coding_ids-coding_ids)[:10]}"
        )

    mismatch_counts: Counter[str] = Counter()
    mismatch_examples: dict[str, list[str]] = defaultdict(list)
    for transcript_id in sorted(coding_ids):
        g3 = gff_by_transcript[transcript_id]
        g4 = gtf_by_transcript[transcript_id]
        gm = transcript_meta[transcript_id]
        tm = gtf_meta[transcript_id]
        checks = {
            "transcript_span_match": (gm["seqid"], gm["start"], gm["end"], gm["strand"]) ==
                                     (tm["seqid"], tm["start"], tm["end"], tm["strand"]),
            "gene_id_match": norm_id(str(gm["gene_id"]), "gene:") == tm["gene_id"],
            "exon_coordinates_match": feature_coords(g3["exon"]) == feature_coords(g4["exon"]),
            "five_prime_utr_match": feature_coords(g3["five_prime_UTR"]) == feature_coords(g4["five_prime_utr"]),
            "three_prime_utr_match": feature_coords(g3["three_prime_UTR"]) == feature_coords(g4["three_prime_utr"]),
            "cds_coverage_match_after_gtf_start_stop": merged_intervals(g3["CDS"]) ==
                merged_intervals(g4["CDS"] + g4["start_codon"] + g4["stop_codon"]),
        }
        for check, passed in checks.items():
            if not passed:
                mismatch_counts[check] += 1
                if len(mismatch_examples[check]) < 10:
                    mismatch_examples[check].append(transcript_id)
    if mismatch_counts:
        raise ValueError(f"GFF3/GTF cross-check failed for coding transcripts: {dict(mismatch_counts)}; examples={dict(mismatch_examples)}")

    phase_mismatches_total = 0
    phase_mismatch_ids: list[str] = []
    parent_gene_issues: list[str] = []
    rows: list[dict[str, str]] = []
    proteins: list[tuple[str, str]] = []
    seen_protein_ids: set[str] = set()
    quality_counts: Counter[str] = Counter()
    first_phase_counts: Counter[str] = Counter()
    gff3_coding_id_set = set(coding_ids)

    for transcript_id in sorted(coding_ids):
        meta = transcript_meta[transcript_id]
        g3 = gff_by_transcript[transcript_id]
        g4 = gtf_by_transcript[transcript_id]
        parts = g3["CDS"]
        if not parts:
            raise ValueError(f"Protein-coding transcript has no GFF3 CDS: {transcript_id}")
        strand = str(meta["strand"])
        seqid = str(meta["seqid"])
        if any(row["seqid"] != seqid or row["strand"] != strand for row in parts):
            raise ValueError(f"CDS segments disagree with transcript sequence/strand: {transcript_id}")
        ordered = coding_order(parts, strand)
        protein_ids = {str(row["attrs"].get("protein_id", "")) for row in ordered if row["attrs"].get("protein_id")}
        if len(protein_ids) > 1:
            raise ValueError(f"Multiple protein IDs on CDS segments for {transcript_id}: {sorted(protein_ids)}")
        protein_id = next(iter(protein_ids), transcript_id)
        if protein_id in seen_protein_ids:
            raise ValueError(f"Duplicate translated protein identifier: {protein_id}")
        seen_protein_ids.add(protein_id)

        first_phase = int(ordered[0]["phase"]) if str(ordered[0]["phase"]) in {"0", "1", "2"} else 0
        first_phase_counts[str(first_phase)] += 1
        cumulative_nt = 0
        phase_consistent = True
        for index, row in enumerate(ordered):
            phase_text = str(row["phase"])
            if phase_text not in {"0", "1", "2"}:
                phase_consistent = False
                continue
            phase = int(phase_text)
            if index > 0:
                expected_phase = (3 - (cumulative_nt % 3)) % 3
                if phase != expected_phase:
                    phase_consistent = False
            cumulative_nt += int(row["end"]) - int(row["start"]) + 1
        if not phase_consistent:
            phase_mismatches_total += 1
            if len(phase_mismatch_ids) < 25:
                phase_mismatch_ids.append(transcript_id)

        exons = g3["exon"]
        cds_exon_containment = all(
            any(exon["seqid"] == part["seqid"] and exon["strand"] == part["strand"] and
                int(exon["start"]) <= int(part["start"]) and int(exon["end"]) >= int(part["end"])
                for exon in exons)
            for part in ordered
        )
        utrs = g3["five_prime_UTR"] + g3["three_prime_UTR"]
        utr_exon_containment = all(
            any(exon["seqid"] == utr["seqid"] and exon["strand"] == utr["strand"] and
                int(exon["start"]) <= int(utr["start"]) and int(exon["end"]) >= int(utr["end"])
                for exon in exons)
            for utr in utrs
        )
        utr_cds_overlap = any(
            utr["seqid"] == cds["seqid"] and max(int(utr["start"]), int(cds["start"])) <= min(int(utr["end"]), int(cds["end"]))
            for utr in utrs for cds in ordered
        )
        exon_spans = sorted((int(row["start"]), int(row["end"])) for row in exons)
        exon_overlap = any(exon_spans[i][0] <= exon_spans[i - 1][1] for i in range(1, len(exon_spans)))
        intron_count = max(0, len(exon_spans) - 1)
        intron_nt = sum(max(0, exon_spans[i][0] - exon_spans[i - 1][1] - 1) for i in range(1, len(exon_spans)))
        transcript_span_ok = all(
            str(row["seqid"]) == seqid and int(meta["start"]) <= int(row["start"]) <= int(row["end"]) <= int(meta["end"])
            for row in exons + ordered + utrs
        )
        if not cds_exon_containment: quality_counts["CDS_not_contained_in_exon"] += 1
        if not utr_exon_containment: quality_counts["UTR_not_contained_in_exon"] += 1
        if utr_cds_overlap: quality_counts["UTR_overlaps_CDS"] += 1
        if exon_overlap: quality_counts["overlapping_exons"] += 1
        if not transcript_span_ok: quality_counts["features_outside_transcript_span"] += 1

        cds_nt = extract_segments(ordered, strand, genome, aliases)
        coding_nt = cds_nt[first_phase:]
        terminal_remainder = len(coding_nt) % 3
        complete_codon_nt = coding_nt[:-terminal_remainder] if terminal_remainder else coding_nt
        terminal_codon = complete_codon_nt[-3:] if len(complete_codon_nt) >= 3 else ""
        terminal_stop_present = terminal_codon in STOP_CODONS
        protein, ambiguous_codons = translate(complete_codon_nt)
        if protein.endswith("*"):
            protein = protein[:-1]
        internal_stops = protein.count("*")
        if internal_stops:
            quality_counts["transcripts_with_internal_stop_codons"] += 1
            protein = protein.replace("*", "X")
        if terminal_remainder:
            quality_counts["terminal_partial_codon"] += 1
        if not terminal_stop_present:
            quality_counts["terminal_stop_absent_from_cds"] += 1
        if first_phase:
            quality_counts["partial_5prime_phase"] += 1
        if ambiguous_codons:
            quality_counts["ambiguous_codons_in_translation"] += 1

        start_feature_seq = extract_segments(g4["start_codon"], strand, genome, aliases) if g4["start_codon"] else ""
        stop_feature_seq = extract_segments(g4["stop_codon"], strand, genome, aliases) if g4["stop_codon"] else ""
        start_codon_match = bool(start_feature_seq) and coding_nt.startswith(start_feature_seq)
        stop_codon_match = bool(stop_feature_seq) and coding_nt.endswith(stop_feature_seq)
        if g4["start_codon"] and not start_codon_match:
            quality_counts["start_codon_sequence_mismatch"] += 1
        if g4["stop_codon"] and not stop_codon_match:
            quality_counts["stop_codon_sequence_mismatch"] += 1
        if not g4["start_codon"]:
            quality_counts["start_codon_feature_missing"] += 1
        if not g4["stop_codon"]:
            quality_counts["stop_codon_feature_missing"] += 1

        gene_id = norm_id(str(meta["gene_id"]), "gene:")
        gene = gene_meta.get(gene_id, {})
        if not gene:
            parent_gene_issues.append(transcript_id)
        if not protein:
            raise ValueError(f"Empty translated protein sequence for {transcript_id}")
        fasta_header = f"{protein_id} transcript_id={transcript_id} gene_id={gene_id}"
        proteins.append((fasta_header, protein))

        flags = []
        if not phase_consistent: flags.append("phase_mismatch")
        if first_phase: flags.append("partial_5prime_phase")
        if terminal_remainder: flags.append("terminal_partial_codon")
        if internal_stops: flags.append("internal_stop")
        if not terminal_stop_present: flags.append("no_terminal_stop_in_CDS")
        if g4["start_codon"] and not start_codon_match: flags.append("start_codon_mismatch")
        if g4["stop_codon"] and not stop_codon_match: flags.append("stop_codon_mismatch")
        if not g4["start_codon"]: flags.append("start_codon_feature_missing")
        if not g4["stop_codon"]: flags.append("stop_codon_feature_missing")
        if ambiguous_codons: flags.append("ambiguous_codons")
        if not cds_exon_containment: flags.append("CDS_outside_exon")
        if not utr_exon_containment: flags.append("UTR_outside_exon")
        if utr_cds_overlap: flags.append("UTR_CDS_overlap")
        if exon_overlap: flags.append("overlapping_exons")
        if not transcript_span_ok: flags.append("outside_transcript_span")
        if not gene: flags.append("gene_parent_missing")
        quality_counts["models_with_any_flag" if flags else "models_without_flags"] += 1

        rows.append({
            "assembly_accession": ASSEMBLY,
            "protein_accession": protein_id,
            "transcript_id": transcript_id,
            "gene_id": gene_id,
            "gene_name": str(gene.get("name", "")),
            "gene_description": str(gene.get("description", "")),
            "gff3_seqid": seqid,
            "genbank_accession": aliases[seqid],
            "strand": strand,
            "gene_start": str(gene.get("start", "")),
            "gene_end": str(gene.get("end", "")),
            "transcript_start": str(meta["start"]),
            "transcript_end": str(meta["end"]),
            "exon_count": str(len(exons)),
            "inferred_intron_count": str(intron_count),
            "inferred_intron_total_nt": str(intron_nt),
            "five_prime_UTR_segments": str(len(g3["five_prime_UTR"])),
            "three_prime_UTR_segments": str(len(g3["three_prime_UTR"])),
            "CDS_segment_count": str(len(ordered)),
            "GTF_exons_match_GFF3": "yes",
            "GTF_UTRs_match_GFF3": "yes",
            "GTF_CDS_plus_start_stop_match_GFF3": "yes",
            "CDS_phase_consistent": "yes" if phase_consistent else "no",
            "first_CDS_phase": str(first_phase),
            "CDS_length_nt_before_phase": str(len(cds_nt)),
            "CDS_length_nt_after_initial_phase": str(len(coding_nt)),
            "terminal_partial_nt_trimmed_for_translation": str(terminal_remainder),
            "GTF_start_codon_feature_present": "yes" if g4["start_codon"] else "no",
            "GTF_start_codon_matches_CDS": "yes" if start_codon_match else "no",
            "GTF_stop_codon_feature_present": "yes" if g4["stop_codon"] else "no",
            "GTF_stop_codon_matches_CDS": "yes" if stop_codon_match else "no",
            "terminal_stop_codon_in_CDS": "yes" if terminal_stop_present else "no",
            "internal_stop_codons_converted_to_X": str(internal_stops),
            "ambiguous_codons_translated_as_X": str(ambiguous_codons),
            "protein_length_aa": str(len(protein)),
            "model_flags": ";".join(flags),
        })

    if phase_mismatches_total:
        raise ValueError(f"CDS phase continuity failed for {phase_mismatches_total} transcripts; examples={phase_mismatch_ids}")
    if parent_gene_issues:
        raise ValueError(f"Protein-coding transcripts with missing gene parents: {len(parent_gene_issues)}; examples={parent_gene_issues[:10]}")
    if len(proteins) != len(coding_ids):
        raise AssertionError(f"Translated {len(proteins)} proteins, expected {len(coding_ids)}")

    gff_unannotated_records = sorted(set(report_by_accession) - {aliases[s] for s in gff_seqids})
    feature_audit = {
        "assembly_accession": ASSEMBLY,
        "genome_fasta_records": len(genome),
        "sequence_report_records": len(report_by_accession),
        "genome_sequence_report_lengths_match": True,
        "gff3_sequence_ids": sorted(gff_seqids),
        "gtf_sequence_ids": sorted(gtf_seqids),
        "gff3_gtf_sequence_id_sets_match": gff_seqids == gtf_seqids,
        "annotated_genome_accessions": sorted({aliases[s] for s in gff_seqids}),
        "sequence_report_records_without_GFF3_annotation": gff_unannotated_records,
        "GFF3_feature_counts": dict(feature_counts_gff),
        "GTF_feature_counts": dict(gtf_counts),
        "repeated_GFF3_CDS_ID_rows_for_discontinuous_features": duplicate_cds_id_rows,
        "GFF3_protein_coding_transcripts": len(coding_ids),
        "GTF_protein_coding_transcripts": len(gtf_coding_ids),
        "GFF3_GTF_coding_transcript_IDs_match": True,
        "GFF3_GTF_exon_and_UTR_coordinates_match": True,
        "GFF3_GTF_CDS_coverage_matches_after_adding_GTF_start_stop_features": True,
        "GFF3_GTF_transcript_and_gene_IDs_match": True,
        "GFF3_GTF_mismatches": dict(mismatch_counts),
        "first_CDS_phase_counts": dict(first_phase_counts),
        "transcripts_with_internal_phase_mismatches": phase_mismatches_total,
        "translation_quality_counts": dict(quality_counts),
        "interpretation": "Proteins represent protein-coding transcript models present in the GFF3. GFF3/GTF cover the 11 chromosome seqids; the sequence report has additional unlocalized/unplaced and organelle records without annotation in these files.",
    }

    outdir.mkdir(parents=True, exist_ok=False)
    protein_path = outdir / "suecica_annotated_proteins.faa"
    manifest_path = outdir / "suecica_protein_transcript_manifest.tsv"
    audit_path = outdir / "annotation_feature_audit.json"
    run_path = outdir / "protein_extraction_manifest.json"
    with protein_path.open("x", encoding="ascii") as handle:
        for header, sequence in proteins:
            write_fasta_record(handle, header, sequence)
    fields = list(rows[0])
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with audit_path.open("x", encoding="utf-8") as handle:
        json.dump(feature_audit, handle, indent=2)
        handle.write("\n")

    run_manifest = {
        "assembly_accession": ASSEMBLY,
        "annotation_primary": str(GFF3.relative_to(PROJECT)),
        "annotation_crosscheck": str(GTF.relative_to(PROJECT)),
        "translation_method": "Concatenate GFF3 CDS segments in transcriptional order; trim only the first CDS phase for a partial 5-prime model and an incomplete terminal codon; translate standard nuclear genetic code; remove one terminal stop; encode internal stop/ambiguous codons as X and flag them in the transcript manifest.",
        "sequence_report_mapping": "GFF3/GTF chromosome labels are mapped through sequence_report.jsonl genbankAccession entries to genome FASTA identifiers.",
        "input_sha256": {
            str(path.relative_to(PROJECT)): sha256(path) for path in (GENOME, GFF3, GTF, SEQ_REPORT)
        },
        "protein_coding_transcript_count": len(coding_ids),
        "translated_protein_count": len(proteins),
        "protein_fasta": str(protein_path.relative_to(PROJECT)),
        "protein_fasta_sha256": sha256(protein_path),
        "transcript_manifest": str(manifest_path.relative_to(PROJECT)),
        "transcript_manifest_sha256": sha256(manifest_path),
        "feature_audit": str(audit_path.relative_to(PROJECT)),
        "feature_audit_sha256": sha256(audit_path),
        "python_version": sys.version.split()[0],
        "external_data_transfer": False,
        "limitations": feature_audit["interpretation"],
    }
    with run_path.open("x", encoding="utf-8") as handle:
        json.dump(run_manifest, handle, indent=2)
        handle.write("\n")

    print(f"Translated {len(proteins)} protein-coding transcript models.")
    print(f"GFF3/GTF cross-checks: exon/UTR coordinates exact; CDS coverage exact after adding GTF start/stop features; phase mismatches={phase_mismatches_total}.")
    print(f"Sequence report records without annotation in these files: {len(gff_unannotated_records)}.")
    print(f"Protein FASTA and QC manifests written under {outdir}")


if __name__ == "__main__":
    main()
