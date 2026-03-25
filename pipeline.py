# ============================================
# Mini Bioinformatics Pipeline
# Author: Padma Shree
# Phase 1 - Project 5
# Combines: DNA Analysis + Transcription +
#           Mutation Detection + FASTA Parsing
# ============================================

import os
import sys
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from datetime import datetime

# ─────────────────────────────────────────────
# PROGRESS BAR
# ─────────────────────────────────────────────

def progress_bar(step, total, label=""):
    percent  = int((step / total) * 100)
    filled   = int(percent / 5)
    bar      = "█" * filled + "░" * (20 - filled)
    print(f"\r  [{bar}] {percent}% — {label}", end="", flush=True)
    if step == total:
        print()

# ─────────────────────────────────────────────
# MODULE 1 — DNA ANALYSIS
# ─────────────────────────────────────────────

def analyze_dna(sequence):
    sequence = sequence.upper()
    bases    = {
        "A": sequence.count("A"),
        "T": sequence.count("T"),
        "G": sequence.count("G"),
        "C": sequence.count("C")
    }
    gc       = round(((bases["G"] + bases["C"]) / len(sequence)) * 100, 2)
    at_gc    = round((bases["A"] + bases["T"]) / max(bases["G"] + bases["C"], 1), 2)
    comp     = sequence.translate(str.maketrans("ATGC", "TACG"))
    rev_comp = comp[::-1]
    return {
        "sequence"   : sequence,
        "length"     : len(sequence),
        "bases"      : bases,
        "gc"         : gc,
        "at_gc"      : at_gc,
        "complement" : comp,
        "rev_comp"   : rev_comp,
        "interpretation": "GC-rich" if gc > 60 else "AT-rich" if gc < 40 else "Balanced"
    }

# ─────────────────────────────────────────────
# MODULE 2 — TRANSCRIPTION & TRANSLATION
# ─────────────────────────────────────────────

CODON_TABLE = {
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met(START)",
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "STOP", "UAG": "STOP",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "UGU": "Cys", "UGC": "Cys", "UGA": "STOP", "UGG": "Trp",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}

def transcribe_translate(dna):
    rna     = dna.upper().replace("T", "U")
    protein = []
    for i in range(0, len(rna) - 2, 3):
        codon = rna[i:i+3]
        if len(codon) == 3:
            aa = CODON_TABLE.get(codon, "?")
            protein.append((codon, aa))
            if aa == "STOP":
                break
    return {
        "rna"        : rna,
        "protein"    : protein,
        "stop_found" : any(aa == "STOP" for _, aa in protein),
        "codon_count": len(protein)
    }

# ─────────────────────────────────────────────
# MODULE 3 — MUTATION DETECTION
# ─────────────────────────────────────────────

def detect_mutations(seq1, seq2):
    seq1   = seq1.upper()
    seq2   = seq2.upper()
    snps   = []
    length = min(len(seq1), len(seq2))
    purines     = set("AG")
    pyrimidines = set("CT")

    for i in range(length):
        if seq1[i] != seq2[i]:
            o, m = seq1[i], seq2[i]
            if (o in purines and m in purines) or (o in pyrimidines and m in pyrimidines):
                mtype = "Transition"
            else:
                mtype = "Transversion"
            snps.append({"position": i+1, "original": o, "mutated": m, "type": mtype})

    indel = len(seq1) - len(seq2)
    return {
        "snps"          : snps,
        "total_snps"    : len(snps),
        "transitions"   : sum(1 for s in snps if s["type"] == "Transition"),
        "transversions" : sum(1 for s in snps if s["type"] == "Transversion"),
        "indel"         : indel,
        "mutation_rate" : round((len(snps) / length) * 100, 2) if length > 0 else 0,
        "similarity"    : round(((length - len(snps)) / length) * 100, 2) if length > 0 else 0
    }

# ─────────────────────────────────────────────
# MODULE 4 — FASTA PARSER
# ─────────────────────────────────────────────

def parse_fasta(filepath):
    sequences = []
    header    = None
    seq_lines = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    sequences.append((header, "".join(seq_lines)))
                header    = line[1:]
                seq_lines = []
            elif line:
                seq_lines.append(line.upper())
        if header is not None:
            sequences.append((header, "".join(seq_lines)))
    return sequences

# ─────────────────────────────────────────────
# VISUALIZATION
# ─────────────────────────────────────────────

def generate_charts(dna_result, trans_result, mut_result, output_dir):
    fig = plt.figure(figsize=(18, 12))
    fig.suptitle("Mini Bioinformatics Pipeline — Analysis Dashboard",
                fontsize=16, fontweight="bold")
    gs  = gridspec.GridSpec(2, 3, figure=fig)

    # Chart 1 — Base composition
    ax1 = fig.add_subplot(gs[0, 0])
    bases   = dna_result["bases"]
    colors  = ["#4CAF50", "#2196F3", "#FF9800", "#E91E63"]
    bars    = ax1.bar(bases.keys(), bases.values(), color=colors, edgecolor="black")
    for bar, val in zip(bars, bases.values()):
        ax1.text(bar.get_x() + bar.get_width()/2,
                bar.get_height() + 0.2, str(val),
                ha="center", fontweight="bold")
    ax1.set_title("Base Composition", fontweight="bold")
    ax1.set_ylabel("Count")

    # Chart 2 — GC content gauge
    ax2 = fig.add_subplot(gs[0, 1])
    gc  = dna_result["gc"]
    color = "#E91E63" if gc > 60 else "#4CAF50" if gc < 40 else "#FF9800"
    ax2.bar(["GC Content"], [gc], color=color, edgecolor="black", width=0.4)
    ax2.bar(["GC Content"], [100 - gc], bottom=[gc],
           color="#f0f0f0", edgecolor="black", width=0.4)
    ax2.text(0, gc/2, f"{gc}%", ha="center", va="center",
            fontsize=16, fontweight="bold", color="white")
    ax2.set_title("GC Content", fontweight="bold")
    ax2.set_ylim(0, 110)
    ax2.set_ylabel("Percentage (%)")

    # Chart 3 — Amino acid distribution
    ax3 = fig.add_subplot(gs[0, 2])
    if trans_result["protein"]:
        aa_counts = {}
        for _, aa in trans_result["protein"]:
            if aa not in ("STOP", "?"):
                aa_counts[aa] = aa_counts.get(aa, 0) + 1
        if aa_counts:
            ax3.bar(aa_counts.keys(), aa_counts.values(),
                   color=plt.cm.Set3.colors[:len(aa_counts)],
                   edgecolor="black")
            plt.setp(ax3.xaxis.get_majorticklabels(), rotation=45, ha="right")
    ax3.set_title("Amino Acid Frequency", fontweight="bold")
    ax3.set_ylabel("Count")

    # Chart 4 — Mutation positions
    ax4 = fig.add_subplot(gs[1, 0])
    if mut_result and mut_result["snps"]:
        positions = [s["position"] for s in mut_result["snps"]]
        colors_m  = ["#E91E63" if s["type"] == "Transversion"
                    else "#2196F3" for s in mut_result["snps"]]
        ax4.scatter(positions, [1]*len(positions),
                   c=colors_m, s=120, edgecolors="black", zorder=5)
        ax4.set_xlim(0, max(positions) + 2)
    else:
        ax4.text(0.5, 0.5, "No mutations / not compared",
                ha="center", va="center", fontsize=10)
    ax4.set_title("Mutation Positions", fontweight="bold")
    ax4.set_xlabel("Position")
    ax4.set_yticks([])

    # Chart 5 — Similarity vs mutation rate
    ax5 = fig.add_subplot(gs[1, 1])
    if mut_result:
        ax5.bar(["Similarity", "Mutation Rate"],
               [mut_result["similarity"], mut_result["mutation_rate"]],
               color=["#4CAF50", "#F44336"], edgecolor="black")
        for i, val in enumerate([mut_result["similarity"], mut_result["mutation_rate"]]):
            ax5.text(i, val + 0.5, f"{val}%", ha="center", fontweight="bold")
        ax5.set_ylim(0, 110)
    else:
        ax5.text(0.5, 0.5, "No reference sequence provided",
                ha="center", va="center", fontsize=10)
    ax5.set_title("Similarity vs Mutation Rate", fontweight="bold")
    ax5.set_ylabel("Percentage (%)")

    # Chart 6 — Pipeline summary
    ax6 = fig.add_subplot(gs[1, 2])
    modules = ["DNA\nAnalysis", "Transcription", "Translation", "Mutation\nDetection"]
    values  = [1, 1, 1, 1 if mut_result else 0]
    colors6 = ["#4CAF50" if v == 1 else "#E0E0E0" for v in values]
    ax6.bar(modules, values, color=colors6, edgecolor="black")
    ax6.set_title("Pipeline Modules Completed", fontweight="bold")
    ax6.set_ylim(0, 1.5)
    ax6.set_yticks([])
    for i, (mod, val) in enumerate(zip(modules, values)):
        ax6.text(i, 0.5, "Done" if val == 1 else "Skipped",
                ha="center", va="center",
                fontsize=10, fontweight="bold", color="white" if val == 1 else "gray")

    plt.tight_layout()
    chart_file = os.path.join(output_dir, "pipeline_dashboard.png")
    plt.savefig(chart_file, dpi=150, bbox_inches="tight")
    plt.close()
    return chart_file

# ─────────────────────────────────────────────
# HTML REPORT
# ─────────────────────────────────────────────

def generate_html_report(dna_result, trans_result, mut_result, output_dir, sequence):
    timestamp = datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    snp_rows  = ""
    if mut_result:
        for s in mut_result["snps"]:
            snp_rows += f"""
            <tr>
                <td>{s['position']}</td>
                <td>{s['original']}</td>
                <td>{s['mutated']}</td>
                <td><span class="badge {'transition' if s['type']=='Transition' else 'transversion'}">{s['type']}</span></td>
            </tr>"""

    protein_rows = ""
    for codon, aa in trans_result["protein"]:
        cls = "stop" if aa == "STOP" else "start" if "START" in aa else ""
        protein_rows += f"<span class='codon {cls}' title='{aa}'>{codon}</span>"

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Bioinformatics Pipeline Report</title>
<style>
  body {{ font-family: Arial, sans-serif; max-width: 1000px; margin: 0 auto; padding: 20px; background: #f5f5f5; }}
  h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
  h2 {{ color: #2980b9; margin-top: 30px; }}
  .card {{ background: white; border-radius: 8px; padding: 20px; margin: 15px 0; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
  .stats-grid {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 15px; }}
  .stat-box {{ background: #3498db; color: white; padding: 15px; border-radius: 8px; text-align: center; }}
  .stat-box.green {{ background: #27ae60; }}
  .stat-box.orange {{ background: #e67e22; }}
  .stat-box.red {{ background: #e74c3c; }}
  .stat-value {{ font-size: 28px; font-weight: bold; }}
  .stat-label {{ font-size: 12px; opacity: 0.9; margin-top: 5px; }}
  table {{ width: 100%; border-collapse: collapse; margin-top: 10px; }}
  th {{ background: #2c3e50; color: white; padding: 10px; text-align: left; }}
  td {{ padding: 8px 10px; border-bottom: 1px solid #eee; }}
  tr:hover {{ background: #f9f9f9; }}
  .badge {{ padding: 3px 8px; border-radius: 4px; font-size: 12px; color: white; }}
  .transition {{ background: #2196F3; }}
  .transversion {{ background: #E91E63; }}
  .codon {{ display: inline-block; padding: 3px 6px; margin: 2px; border-radius: 4px; background: #ecf0f1; font-family: monospace; font-size: 13px; }}
  .codon.start {{ background: #27ae60; color: white; }}
  .codon.stop {{ background: #e74c3c; color: white; }}
  .sequence-box {{ font-family: monospace; word-break: break-all; background: #2c3e50; color: #2ecc71; padding: 15px; border-radius: 8px; font-size: 13px; line-height: 1.8; }}
  .footer {{ text-align: center; color: #7f8c8d; margin-top: 30px; font-size: 13px; }}
  img {{ max-width: 100%; border-radius: 8px; margin-top: 10px; }}
</style>
</head>
<body>
<h1>Mini Bioinformatics Pipeline Report</h1>
<p><strong>Generated:</strong> {timestamp} | <strong>Author:</strong> Padma Shree</p>

<div class="card">
  <h2>Pipeline Summary</h2>
  <div class="stats-grid">
    <div class="stat-box">
      <div class="stat-value">{dna_result['length']}</div>
      <div class="stat-label">Total Bases</div>
    </div>
    <div class="stat-box green">
      <div class="stat-value">{dna_result['gc']}%</div>
      <div class="stat-label">GC Content</div>
    </div>
    <div class="stat-box orange">
      <div class="stat-value">{trans_result['codon_count']}</div>
      <div class="stat-label">Codons Translated</div>
    </div>
    <div class="stat-box {'red' if mut_result and mut_result['total_snps'] > 0 else 'green'}">
      <div class="stat-value">{mut_result['total_snps'] if mut_result else 'N/A'}</div>
      <div class="stat-label">Mutations Detected</div>
    </div>
    <div class="stat-box">
      <div class="stat-value">{mut_result['similarity'] if mut_result else 'N/A'}{'%' if mut_result else ''}</div>
      <div class="stat-label">Sequence Similarity</div>
    </div>
    <div class="stat-box {'red' if mut_result and mut_result['mutation_rate'] > 10 else 'green'}">
      <div class="stat-value">{mut_result['mutation_rate'] if mut_result else 'N/A'}{'%' if mut_result else ''}</div>
      <div class="stat-label">Mutation Rate</div>
    </div>
  </div>
</div>

<div class="card">
  <h2>Module 1 — DNA Sequence Analysis</h2>
  <div class="sequence-box">{dna_result['sequence']}</div>
  <table>
    <tr><th>Property</th><th>Value</th></tr>
    <tr><td>Length</td><td>{dna_result['length']} bases</td></tr>
    <tr><td>GC Content</td><td>{dna_result['gc']}%</td></tr>
    <tr><td>AT/GC Ratio</td><td>{dna_result['at_gc']}</td></tr>
    <tr><td>A count</td><td>{dna_result['bases']['A']}</td></tr>
    <tr><td>T count</td><td>{dna_result['bases']['T']}</td></tr>
    <tr><td>G count</td><td>{dna_result['bases']['G']}</td></tr>
    <tr><td>C count</td><td>{dna_result['bases']['C']}</td></tr>
    <tr><td>Complement</td><td style="font-family:monospace">{dna_result['complement']}</td></tr>
    <tr><td>Reverse Complement</td><td style="font-family:monospace">{dna_result['rev_comp']}</td></tr>
    <tr><td>Interpretation</td><td>{dna_result['interpretation']}</td></tr>
  </table>
</div>

<div class="card">
  <h2>Module 2 — Transcription and Translation</h2>
  <p><strong>RNA Transcript:</strong></p>
  <div class="sequence-box">{trans_result['rna']}</div>
  <p><strong>Protein Codons:</strong></p>
  <div>{protein_rows}</div>
  <p>STOP codon found: <strong>{'Yes' if trans_result['stop_found'] else 'No'}</strong> | Total codons: <strong>{trans_result['codon_count']}</strong></p>
</div>

<div class="card">
  <h2>Module 3 — Mutation Detection</h2>
  {'<p>No reference sequence provided for mutation analysis.</p>' if not mut_result else f"""
  <div class="stats-grid">
    <div class="stat-box"><div class="stat-value">{mut_result['total_snps']}</div><div class="stat-label">Total SNPs</div></div>
    <div class="stat-box green"><div class="stat-value">{mut_result['transitions']}</div><div class="stat-label">Transitions</div></div>
    <div class="stat-box red"><div class="stat-value">{mut_result['transversions']}</div><div class="stat-label">Transversions</div></div>
  </div>
  <table>
    <tr><th>Position</th><th>Original</th><th>Mutated</th><th>Type</th></tr>
    {snp_rows if snp_rows else '<tr><td colspan="4">No SNPs detected</td></tr>'}
  </table>"""}
</div>

<div class="card">
  <h2>Analysis Dashboard</h2>
  <img src="pipeline_dashboard.png" alt="Pipeline Dashboard">
</div>

<div class="footer">
  <p>Generated by Mini Bioinformatics Pipeline | Author: Padma Shree | Built with guidance from Claude AI (Anthropic)</p>
</div>
</body>
</html>"""

    report_file = os.path.join(output_dir, "pipeline_report.html")
    with open(report_file, "w") as f:
        f.write(html)
    return report_file

# ─────────────────────────────────────────────
# MAIN PIPELINE
# ─────────────────────────────────────────────

def run_pipeline(sequence, reference=None, fasta_file=None):
    timestamp  = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = f"pipeline_output_{timestamp}"
    os.makedirs(output_dir, exist_ok=True)

    total_steps = 6
    print("\n" + "="*60)
    print("  Mini Bioinformatics Pipeline")
    print("  Author: Padma Shree | Phase 1 - Project 5")
    print("="*60)
    print(f"\n  Output folder: {output_dir}")
    print("\n  Starting pipeline...\n")

    # Step 1
    progress_bar(1, total_steps, "Validating sequence...")
    sequence = sequence.upper().strip()

    # Step 2
    progress_bar(2, total_steps, "Running DNA analysis...")
    dna_result = analyze_dna(sequence)

    # Step 3
    progress_bar(3, total_steps, "Transcribing and translating...")
    trans_result = transcribe_translate(sequence)

    # Step 4
    progress_bar(4, total_steps, "Detecting mutations...")
    mut_result = detect_mutations(sequence, reference) if reference else None

    # Step 5
    progress_bar(5, total_steps, "Generating charts...")
    chart_file = generate_charts(dna_result, trans_result, mut_result, output_dir)

    # Step 6
    progress_bar(6, total_steps, "Generating HTML report...")
    report_file = generate_html_report(dna_result, trans_result,
                                       mut_result, output_dir, sequence)

    # Print summary
    print("\n" + "="*60)
    print("  PIPELINE COMPLETE!")
    print("="*60)
    print(f"\n  Sequence length    : {dna_result['length']} bases")
    print(f"  GC content         : {dna_result['gc']}%")
    print(f"  Interpretation     : {dna_result['interpretation']}")
    print(f"  RNA transcript     : {trans_result['rna'][:40]}...")
    print(f"  Codons translated  : {trans_result['codon_count']}")
    print(f"  STOP codon found   : {'Yes' if trans_result['stop_found'] else 'No'}")
    if mut_result:
        print(f"  Mutations detected : {mut_result['total_snps']}")
        print(f"  Mutation rate      : {mut_result['mutation_rate']}%")
        print(f"  Similarity         : {mut_result['similarity']}%")
    print(f"\n  Chart saved        : {chart_file}")
    print(f"  HTML report saved  : {report_file}")
    print("\n" + "="*60)
    print(f"  All files saved in: {output_dir}/")
    print("="*60)

    return report_file

# ─────────────────────────────────────────────
# ARGUMENT PARSER
# ─────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Mini Bioinformatics Pipeline by Padma Shree"
    )
    parser.add_argument("--sequence",  "-s", help="Input DNA sequence")
    parser.add_argument("--reference", "-r", help="Reference sequence for mutation detection")
    parser.add_argument("--fasta",     "-f", help="Input FASTA file path")
    args = parser.parse_args()

    if args.sequence:
        run_pipeline(args.sequence, args.reference, args.fasta)
    else:
        print("\n" + "="*60)
        print("  Mini Bioinformatics Pipeline")
        print("  Author: Padma Shree | Phase 1 - Project 5")
        print("="*60)
        print("\nEnter DNA sequence:")
        sequence  = input("Sequence  >>> ").strip()
        print("Enter reference sequence for mutation detection (or press Enter to skip):")
        reference = input("Reference >>> ").strip()
        run_pipeline(sequence, reference if reference else None)