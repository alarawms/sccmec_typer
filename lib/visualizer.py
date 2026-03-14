import os
from jinja2 import Environment, FileSystemLoader

# IWG-SCC color scheme
GENE_COLORS = {
    "mecA": "#DC143C", "mecC": "#DC143C", "mecB": "#DC143C",
    "mecR1": "#E8A0BF", "mecI": "#E8A0BF",
    "blaZ": "#8E44AD",
}
CCR_COLOR = "#7B68EE"
IS_COLOR = "#95A5A6"
TIER_COLORS = {"High": "#27AE60", "Medium": "#F39C12", "Low": "#E74C3C", "n/a": "#95A5A6"}

# Canonical gene order in SCCmec cassette
# Left to right: DR - J3 - [ccr genes] - J2 - [mec genes + IS] - J1 - DR
MEC_GENES = {"mecA", "mecC", "mecB", "mecR1", "mecI", "IS431", "IS431_1", "IS431_2", "IS1272", "blaZ"}
CCR_GENES_PREFIX = "ccr"


def _get_gene_color(gene_name):
    if gene_name in GENE_COLORS:
        return GENE_COLORS[gene_name]
    if gene_name.startswith("ccr"):
        return CCR_COLOR
    if gene_name.startswith("IS"):
        return IS_COLOR
    return "#95A5A6"


def _get_template_env():
    templates_dir = os.path.join(os.path.dirname(__file__), "../templates")
    return Environment(loader=FileSystemLoader(templates_dir), autoescape=False)


def _layout_cassette_genes(result):
    """Compute x/y positions for genes in the schematic cassette diagram.

    Places genes in canonical SCCmec order:
    DR(left) - J3 - [ccr genes] - J2 - [mec complex genes] - J1 - DR(right)

    Returns a list of gene dicts with layout fields added:
    {name, detected, color, opacity, x, y, w, h, strand, score, tier, tier_color,
     contig_badge, arrow_points}
    """
    genes_detected = set(result.get("genes_detected", []))
    per_gene = result.get("confidence", {}).get("per_gene", {})
    is_split = result.get("assembly", {}).get("is_split", False)

    # Build location lookup for contig badges
    gene_contigs = {}
    for loc in result.get("assembly", {}).get("gene_locations", []):
        gene_contigs[loc["gene"]] = loc["contig"]

    # Collect unique contigs for badge abbreviation
    unique_contigs = sorted(set(gene_contigs.values()))
    contig_abbrev = {c: f"c{i+1}" for i, c in enumerate(unique_contigs)}

    # Separate ccr and mec genes
    ccr_genes_list = sorted([g for g in genes_detected if g.startswith("ccr")])
    mec_genes_list = sorted([g for g in genes_detected if g in MEC_GENES])

    # Also include expected but missing genes
    mec_comp = result.get("confidence", {}).get("mec_component")
    ccr_comp = result.get("confidence", {}).get("ccr_component")

    if mec_comp and mec_comp.get("genes_missing"):
        for g in mec_comp["genes_missing"]:
            if g not in mec_genes_list:
                mec_genes_list.append(g)

    if ccr_comp and ccr_comp.get("genes_missing"):
        for g in ccr_comp["genes_missing"]:
            if g not in ccr_genes_list:
                ccr_genes_list.append(g)

    # Layout constants
    svg_width = 900
    gene_h = 24
    gene_y = 90
    dr_w = 30
    j_w = 40
    gene_spacing = 8
    margin = 20

    # Calculate available width for gene blocks
    structural_w = 2 * dr_w + 3 * j_w + 6 * gene_spacing
    available_w = svg_width - 2 * margin - structural_w
    total_genes = len(ccr_genes_list) + len(mec_genes_list)
    gene_w = min(80, available_w / max(total_genes, 1))

    # Position elements left to right
    x = margin
    layout = []

    # DR left
    dr_left = {"points": f"{x},{gene_y - 12} {x + dr_w},{gene_y} {x},{gene_y + 12}"}
    x += dr_w + gene_spacing

    # J3 region
    j3 = {"x": x, "y": gene_y - 15, "w": j_w, "h": 30, "label": "J3"}
    x += j_w + gene_spacing

    # ccr genes
    for gene_name in ccr_genes_list:
        detected = gene_name in genes_detected
        gene_info = per_gene.get(gene_name, {})
        score = gene_info.get("score", 0)
        tier = gene_info.get("tier", "Low")
        opacity = max(0.4, score) if detected else 1.0

        contig_badge = None
        if is_split and gene_name in gene_contigs:
            contig_badge = contig_abbrev.get(gene_contigs[gene_name])

        gx = x
        # Arrow points: right-pointing arrow
        arrow_points = f"{gx},{gene_y - 12} {gx + gene_w - 6},{gene_y - 12} {gx + gene_w},{gene_y} {gx + gene_w - 6},{gene_y + 12} {gx},{gene_y + 12}"

        layout.append({
            "name": gene_name, "detected": detected,
            "color": _get_gene_color(gene_name), "opacity": opacity,
            "x": gx, "y": gene_y - 12, "w": gene_w, "h": gene_h,
            "strand": "+", "score": score, "tier": tier,
            "tier_color": TIER_COLORS.get(tier, "#95A5A6"),
            "contig_badge": contig_badge, "arrow_points": arrow_points,
        })
        x += gene_w + gene_spacing

    # J2 region
    j2 = {"x": x, "y": gene_y - 15, "w": j_w, "h": 30, "label": "J2"}
    x += j_w + gene_spacing

    # mec genes
    for gene_name in mec_genes_list:
        detected = gene_name in genes_detected
        gene_info = per_gene.get(gene_name, {})
        score = gene_info.get("score", 0)
        tier = gene_info.get("tier", "Low")
        opacity = max(0.4, score) if detected else 1.0

        contig_badge = None
        if is_split and gene_name in gene_contigs:
            contig_badge = contig_abbrev.get(gene_contigs[gene_name])

        gx = x
        arrow_points = f"{gx},{gene_y - 12} {gx + gene_w - 6},{gene_y - 12} {gx + gene_w},{gene_y} {gx + gene_w - 6},{gene_y + 12} {gx},{gene_y + 12}"

        layout.append({
            "name": gene_name, "detected": detected,
            "color": _get_gene_color(gene_name), "opacity": opacity,
            "x": gx, "y": gene_y - 12, "w": gene_w, "h": gene_h,
            "strand": "+", "score": score, "tier": tier,
            "tier_color": TIER_COLORS.get(tier, "#95A5A6"),
            "contig_badge": contig_badge, "arrow_points": arrow_points,
        })
        x += gene_w + gene_spacing

    # J1 region
    j1 = {"x": x, "y": gene_y - 15, "w": j_w, "h": 30, "label": "J1"}
    x += j_w + gene_spacing

    # DR right
    dr_right = {"points": f"{x},{gene_y - 12} {x + dr_w},{gene_y} {x},{gene_y + 12}"}

    j_regions = [j3, j2, j1]
    dr_markers = [dr_left, dr_right]

    return layout, j_regions, dr_markers, svg_width


def _layout_coordinate_tracks(result):
    """Compute coordinate map tracks from gene_locations.

    Groups genes by contig, sorts by start position, computes
    proportional x-positions within each track.

    Returns (tracks, total_height, is_split).
    """
    locations = result.get("assembly", {}).get("gene_locations", [])
    if not locations:
        return [], 0, False

    from collections import defaultdict
    contig_genes = defaultdict(list)
    for loc in locations:
        contig_genes[loc["contig"]].append(loc)

    svg_width = 900
    margin = 60
    track_height = 60
    track_spacing = 40
    tracks = []
    y = 30

    for contig in sorted(contig_genes.keys()):
        genes = sorted(contig_genes[contig], key=lambda g: g["start"])
        if not genes:
            continue

        # Find the range of positions on this contig
        min_pos = min(g["start"] for g in genes)
        max_pos = max(g["end"] for g in genes)
        contig_span = max_pos - min_pos if max_pos > min_pos else 1

        track_w = svg_width - 2 * margin
        per_gene = result.get("confidence", {}).get("per_gene", {})

        track_genes = []
        for g in genes:
            gene_info = per_gene.get(g["gene"], {})
            score = gene_info.get("score", 0.5)
            rel_start = (g["start"] - min_pos) / contig_span
            rel_end = (g["end"] - min_pos) / contig_span
            gx = margin + rel_start * track_w
            gw = max(10, (rel_end - rel_start) * track_w)

            track_genes.append({
                "name": g["gene"],
                "start": g["start"],
                "end": g["end"],
                "strand": g.get("strand", "+"),
                "color": _get_gene_color(g["gene"]),
                "opacity": max(0.4, score),
                "x": gx,
                "w": gw,
            })

        tracks.append({
            "contig": contig,
            "length": max_pos,
            "x": margin,
            "y": y + 20,
            "w": track_w,
            "genes": track_genes,
        })
        y += track_height + track_spacing

    is_split = len(tracks) > 1
    total_height = y + (30 if is_split else 10)
    return tracks, total_height, is_split


def generate_cassette_svg(result):
    """Generate the schematic cassette diagram SVG string."""
    env = _get_template_env()

    if result.get("status") == "Negative":
        # Minimal negative SVG
        return env.get_template("cassette.svg.j2").render(
            width=900, height=120,
            sccmec_type="Negative", designation=None,
            mec_complex="Negative", ccr_complex="Negative",
            type_confidence={"score": 0, "tier": "n/a", "color": "#95A5A6"},
            genes=[], j_regions=[], dr_markers=[], component_badges=[],
            warnings=["No SCCmec detected"],
            confidence_note=None,
        )

    genes, j_regions, dr_markers, svg_width = _layout_cassette_genes(result)

    type_conf = result.get("confidence", {}).get("type_level", {"score": 0, "tier": "n/a"})
    type_conf["color"] = TIER_COLORS.get(type_conf["tier"], "#95A5A6")

    # Component badges
    badges = []
    mec_comp = result.get("confidence", {}).get("mec_component")
    ccr_comp = result.get("confidence", {}).get("ccr_component")

    if mec_comp:
        badges.append({
            "x": svg_width / 2 + 100, "y": 140,
            "label": f"mec ({mec_comp['class']})",
            "score": mec_comp["score"], "tier": mec_comp["tier"],
            "color": TIER_COLORS.get(mec_comp["tier"], "#95A5A6"),
        })
    if ccr_comp and ccr_comp.get("type") != "n/a":
        badges.append({
            "x": 100, "y": 140,
            "label": f"ccr ({ccr_comp['type']})",
            "score": ccr_comp["score"], "tier": ccr_comp["tier"],
            "color": TIER_COLORS.get(ccr_comp["tier"], "#95A5A6"),
        })

    conf_mode = result.get("confidence", {}).get("mode", "full")
    conf_note = "Confidence based on coverage only (reads mode)" if conf_mode == "coverage_only" else None

    height = 170 + len(result.get("warnings", [])) * 15
    return env.get_template("cassette.svg.j2").render(
        width=svg_width, height=height,
        sccmec_type=result.get("sccmec_type", "Unknown"),
        designation=result.get("designation"),
        mec_complex=result.get("mec_complex", "Unknown"),
        ccr_complex=result.get("ccr_complex", "Unknown"),
        type_confidence=type_conf,
        genes=genes, j_regions=j_regions, dr_markers=dr_markers,
        component_badges=badges,
        warnings=result.get("warnings", []),
        confidence_note=conf_note,
    )


def generate_coordinate_map_svg(result):
    """Generate the genomic coordinate map SVG string. Returns None for reads mode."""
    if result.get("assembly", {}).get("contigs") == ["Reads"]:
        return None
    if not result.get("assembly", {}).get("gene_locations"):
        return None

    env = _get_template_env()
    tracks, total_height, is_split = _layout_coordinate_tracks(result)

    if not tracks:
        return None

    split_line_y = total_height - 40 if is_split else 0

    return env.get_template("coordinate_map.svg.j2").render(
        width=900, height=total_height,
        tracks=tracks, is_split=is_split, split_line_y=split_line_y,
    )


def generate_report_html(result, sample_name="sample"):
    """Generate the full HTML report with embedded SVGs and gene table."""
    env = _get_template_env()

    cassette_svg = generate_cassette_svg(result)
    coord_svg = generate_coordinate_map_svg(result)

    # Build gene table from per_gene confidence data
    gene_table = []
    per_gene = result.get("confidence", {}).get("per_gene", {})
    locations = {loc["gene"]: loc for loc in result.get("assembly", {}).get("gene_locations", [])}

    for gene_name in sorted(per_gene.keys()):
        info = per_gene[gene_name]
        loc = locations.get(gene_name, {})
        gene_table.append({
            "name": gene_name,
            "id_pct": info.get("identity_pct", 0),
            "cov_pct": info.get("coverage_pct", 0),
            "strand": loc.get("strand", "n/a"),
            "contig": loc.get("contig", "n/a"),
            "start": loc.get("start", 0),
            "end": loc.get("end", 0),
            "score": info.get("score", 0),
            "tier": info.get("tier", "n/a"),
        })

    type_conf = result.get("confidence", {}).get("type_level", {"score": 0, "tier": "n/a"})
    tier_class = type_conf["tier"].lower() if type_conf["tier"] != "n/a" else "na"

    return env.get_template("report.html.j2").render(
        sample_name=sample_name,
        status=result.get("status", "Unknown"),
        sccmec_type=result.get("sccmec_type", "Unknown"),
        designation=result.get("designation"),
        mec_complex=result.get("mec_complex", "Unknown"),
        ccr_complex=result.get("ccr_complex", "Unknown"),
        mecA_present=result.get("mecA_present", False),
        type_confidence=type_conf,
        confidence_tier_class=tier_class,
        confidence_mode=result.get("confidence", {}).get("mode", "full"),
        cassette_svg=cassette_svg,
        coordinate_map_svg=coord_svg,
        gene_table=gene_table if gene_table else None,
        warnings=result.get("warnings", []),
    )


def write_visualization(result, output_prefix, sample_name="sample"):
    """Write SVG and HTML report files.

    Produces:
      {output_prefix}_map.svg    -- Schematic cassette diagram
      {output_prefix}_report.html -- Full HTML report
    """
    svg = generate_cassette_svg(result)
    svg_path = f"{output_prefix}_map.svg"
    with open(svg_path, "w") as f:
        f.write(svg)

    html = generate_report_html(result, sample_name=sample_name)
    html_path = f"{output_prefix}_report.html"
    with open(html_path, "w") as f:
        f.write(html)

    return svg_path, html_path
