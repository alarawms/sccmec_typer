#!/usr/bin/env python3
import argparse
import sys
import os

# Add lib to path
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from lib.aligner import run_minimap2
from lib.parser import parse_paf
from lib.classifier import classify_sccmec

_DEFAULT_BEST_FIT_THRESHOLD = 0.75


def apply_best_fit(result, threshold=_DEFAULT_BEST_FIT_THRESHOLD):
    """Promote the estimator's top candidate when classification returns Unknown.

    When ``--best-fit`` is set and the standard classifier cannot confirm an
    SCCmec type (sccmec_type == "Unknown"), this function uses the estimator's
    best_guess as the final classification, labelled as "Type X (est.)".

    Samples below *threshold* remain "Unknown" — do not lower the threshold
    without understanding the biological evidence, as short-read assemblies of
    CC97/SCCmec V frequently produce split contigs that fool the orientation
    check even when the type is unambiguous from gene content alone.

    Args:
        result: enriched classifier result dict (must include "estimation" key).
        threshold: minimum best_guess_score required to promote (default 0.75).

    Returns:
        Tuple (result, promoted) where *promoted* is True if a promotion
        occurred.  The result dict is modified in-place.
    """
    if result.get("sccmec_type") != "Unknown":
        return result, False

    estimation = result.get("estimation")
    if not estimation:
        return result, False

    score = estimation.get("best_guess_score", 0.0)
    best_guess = estimation.get("best_guess", "")

    if not best_guess or score < threshold:
        return result, False

    original_type = result["sccmec_type"]
    result["sccmec_type"] = f"{best_guess} (est.)"
    result.setdefault("warnings", []).append(
        f"Best-fit estimation promoted from Unknown: "
        f"{best_guess} (score={score:.2f}, threshold={threshold:.2f})"
    )
    result["best_fit_applied"] = {
        "original": original_type,
        "promoted_to": result["sccmec_type"],
        "score": round(score, 3),
        "threshold": threshold,
        "ruling": estimation.get("best_guess_ruling", ""),
    }

    return result, True


def main():
    parser = argparse.ArgumentParser(description="SCCmec Typer: A minimap2-based tool for typing SCCmec elements.")
    parser.add_argument("--1", dest="input1", required=True, help="Input genome assembly (FASTA) or forward reads (FASTQ)")
    parser.add_argument("--2", dest="input2", help="Input reverse reads (FASTQ) for paired-end data")
    parser.add_argument("-d", "--db", required=True, help="Path to reference database (mfa)")
    parser.add_argument("-o", "--output", default="sccmec_results", help="Output prefix")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("--no-viz", action="store_true", help="Skip SVG/HTML visualization generation")
    parser.add_argument(
        "--best-fit",
        dest="best_fit",
        action="store_true",
        default=False,
        help=(
            "When classification returns Unknown (e.g. due to a split assembly "
            "breaking IS431 orientation detection), promote the estimator's top "
            "candidate as the final SCCmec type, labelled 'Type X (est.)'. "
            "Combine with --min-estimate-score to control the confidence floor."
        ),
    )
    parser.add_argument(
        "--min-estimate-score",
        dest="min_estimate_score",
        type=float,
        default=_DEFAULT_BEST_FIT_THRESHOLD,
        metavar="SCORE",
        help=(
            "Minimum estimator score (0–1) required to promote a best-fit type "
            "when --best-fit is active. Candidates below this threshold remain "
            "Unknown. Default: %(default)s."
        ),
    )

    args = parser.parse_args()
    
    print(f"Starting SCCmec Typer on {args.input1}...")
    
    # Determine input type and preset.
    # asm20 (≤20% divergence) improves sensitivity on fragmented assemblies
    # and for mecC/mecB which are genuinely more divergent than mecA.
    preset = "asm20"
    is_reads = False
    
    if args.input2:
        print("Note: Paired-end inputs detected. Switching to Read Mode (Short Reads).")
        preset = "sr"
        is_reads = True
    elif args.input1.lower().endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
        print("Note: Input appears to be FASTQ. Switching to Read Mode.")
        preset = "map-ont" # Default to ONT, could be sr if user specified but for now auto-detect extension? 
        # Ideally we'd have a flag for short reads vs long reads. 
        # Assuming long reads for now as per "sccmec typing" often uses nanopore.
        is_reads = True

    # 1. Alignment
    paf_file = f"{args.output}.paf"
    
    if is_reads:
        # Read Mode: Map Reads (Query) -> Genes (Target)
        run_minimap2(args.input1, args.db, paf_file, args.threads, preset=preset, target_is_db=True, input_file_2=args.input2)
        
        # 2. Parsing (Coverage Calculation)
        from lib.coverage import calculate_gene_coverage
        hits, soft_hits = calculate_gene_coverage(paf_file)
        
    else:
        # Assembly Mode: Map Genes (Query) -> Genome (Target)
        run_minimap2(args.input1, args.db, paf_file, args.threads, preset=preset, target_is_db=False)
        
        # 2. Parsing
        hits, soft_hits = parse_paf(paf_file)
    
    # 3. Classification
    result = classify_sccmec(hits)

    from lib.confidence import enrich_result_with_confidence
    result = enrich_result_with_confidence(result, soft_hits=soft_hits)

    # 4. Best-fit promotion (optional)
    best_fit_applied = False
    if args.best_fit:
        result, best_fit_applied = apply_best_fit(result, threshold=args.min_estimate_score)
        if best_fit_applied:
            bfi = result["best_fit_applied"]
            print(
                f"Best-fit applied: Unknown → {bfi['promoted_to']} "
                f"(score={bfi['score']:.2f})"
            )

    # Output results
    import json
    import csv

    # JSON Output
    json_output_file = f"{args.output}.json"
    with open(json_output_file, 'w') as f:
        json.dump(result, f, indent=4)
    print(f"JSON results written to {json_output_file}")

    # CSV Output (Detailed Elements)
    csv_output_file = f"{args.output}_elements.csv"
    with open(csv_output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        # Get headers from first hit if available, else default
        headers = ["gene", "contig", "start", "end", "strand", "id_pct", "cov_pct", "len", "aln_len"]
        writer.writerow(headers)
        
        for hit in result.get("hits_summary", []):
            writer.writerow([
                hit.get("gene", ""),
                hit.get("contig", ""),
                hit.get("start", ""),
                hit.get("end", ""),
                hit.get("strand", ""),
                f"{hit.get('id_pct', 0):.2f}",
                f"{hit.get('cov_pct', 0):.2f}",
                hit.get("len", ""),
                hit.get("aln_len", "")
            ])
    print(f"Detailed elements CSV written to {csv_output_file}")

    # TSV Output (Summary)
    tsv_output_file = f"{args.output}.tsv"
    with open(tsv_output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        # Header
        writer.writerow(['Sample', 'Status', 'mecA_Present', 'Mec_Gene', 'SCCmec_Type', 'Mec_Complex', 'Ccr_Complex', 'Genes_Detected', 'Warnings', 'Estimated_Type', 'Estimation_Score', 'Best_Fit_Applied'])
        # Data
        sample_name = os.path.basename(args.input1)
        genes_str = ",".join(result.get('genes_detected', []))
        warnings_str = "; ".join(result.get('warnings', []))
        est_type = ""
        est_score = ""
        if result.get("estimation"):
            est_type = result["estimation"].get("best_guess", "")
            est_score = f"{result['estimation'].get('best_guess_score', 0):.2f}"

        # Derive mec gene from element-level hits so mecA/mecC/mecB are reported
        # even when the full cassette cannot be typed (e.g. fragmented assembly).
        all_hit_genes = {h['gene'] for h in result.get('hits_summary', [])}
        all_hit_genes |= {h['gene'] for h in soft_hits}
        if 'mecA' in all_hit_genes:
            mec_gene_detected = 'mecA'
        elif 'mecC' in all_hit_genes:
            mec_gene_detected = 'mecC'
        elif 'mecB' in all_hit_genes:
            mec_gene_detected = 'mecB'
        else:
            mec_gene_detected = 'None'

        meca_present = result.get('mecA_present', False) or mec_gene_detected == 'mecA'

        writer.writerow([
            sample_name,
            result.get('status', 'Unknown'),
            "Yes" if meca_present else "No",
            mec_gene_detected,
            result.get('sccmec_type', 'Unknown'),
            result.get('mec_complex', 'Unknown'),
            result.get('ccr_complex', 'Unknown'),
            genes_str,
            warnings_str,
            est_type,
            est_score,
            "Yes" if best_fit_applied else "No",
        ])
    print(f"TSV summary written to {tsv_output_file}")

    # Visualization output
    if not args.no_viz:
        from lib.visualizer import write_visualization
        sample_name = os.path.basename(args.input1)
        svg_path, html_path = write_visualization(result, args.output, sample_name=sample_name)
        print(f"SVG cassette diagram written to {svg_path}")
        print(f"HTML report written to {html_path}")

if __name__ == "__main__":
    main()
