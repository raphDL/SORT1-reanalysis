"""Clean-room AlphaGenome reconstruction for Figure 3."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .common import Audit, api_key, download
from .figure2 import BASES, CHROM, GENE_TABLE, GENES, LIVER, RS_POS, SEQ_LEN, _ag, _mean_rna
from .figure1 import _save_svg

SCAN_START, SCAN_END = -250, 250
HOTSPOT_LENGTH, HOTSPOT_MARGIN = 12, 8
CEBP_START, CEBP_END = -1, 8


def _native_t(fasta_path: Path, interval: Any) -> str:
    import pysam
    fasta=pysam.FastaFile(str(fasta_path)); sequence=fasta.fetch(CHROM,int(interval.start),int(interval.end)).upper()
    index=RS_POS-1-int(interval.start)
    if sequence[index] != "G": raise ValueError("GRCh38 reference mismatch at rs12740374")
    return sequence[:index]+"T"+sequence[index+1:]


def _summarize(output: Any, interval: Any) -> dict[str,float]:
    table=_mean_rna(output.rna_seq,interval,GENE_TABLE[GENE_TABLE.gene_symbol.isin(GENES)])
    return table.groupby("gene_symbol").rna_mean_tss_pm2kb.mean().to_dict()


def _hotspots(position: pd.DataFrame) -> pd.DataFrame:
    candidates=[]
    for start in range(SCAN_START,SCAN_END-HOTSPOT_LENGTH+2):
        end=start+HOTSPOT_LENGTH-1
        if start<=CEBP_END and end>=CEBP_START: continue
        x=position[position.edit_offset.between(start,end)]
        if len(x)!=HOTSPOT_LENGTH: continue
        candidates.append(dict(start_offset=start,end_offset=end,mean_positive_max_loss=x.positive_max_loss.mean(),sum_positive_max_loss=x.positive_max_loss.sum(),peak_positive_max_loss=x.positive_max_loss.max(),peak_offset=int(x.loc[x.positive_max_loss.idxmax(),"edit_offset"])))
    ranked=pd.DataFrame(candidates).sort_values(["mean_positive_max_loss","peak_positive_max_loss","start_offset"],ascending=[False,False,True])
    selected=[]; occupied=[]
    for _,row in ranked.iterrows():
        if any(int(row.start_offset)<=b+HOTSPOT_MARGIN and int(row.end_offset)>=a-HOTSPOT_MARGIN for a,b in occupied): continue
        selected.append(row); occupied.append((int(row.start_offset),int(row.end_offset)))
        if len(selected)==6: break
    out=pd.DataFrame(selected).sort_values("start_offset").reset_index(drop=True); out.insert(0,"hotspot",[f"H{i}" for i in range(1,7)])
    for c in ("start_offset","end_offset","peak_offset"): out[c]=out[c].astype(int)
    return out


def run_fig3b(run_dir: Path, audit: Audit, fasta_path: Path, *, batch_size: int, max_workers: int) -> None:
    genome,dna_client,dna_model,_=_ag(); interval=genome.Interval(CHROM,RS_POS,RS_POS).resize(SEQ_LEN)
    native=_native_t(fasta_path,interval); client=dna_client.create(api_key(),model_version=dna_model.ModelVersion.ALL_FOLDS,timeout=300)
    root=run_dir/"predictions/Figure3B_native_501bp_ism"; root.mkdir(parents=True,exist_ok=True)
    baseline_path=root/"native_T_baseline.tsv"
    if baseline_path.exists(): baseline=pd.read_csv(baseline_path,sep="\t").set_index("gene").rna_tss.to_dict()
    else:
        output=client.predict_sequence(sequence=native,requested_outputs={dna_client.OutputType.RNA_SEQ},ontology_terms=[LIVER],interval=interval)
        baseline=_summarize(output,interval); pd.DataFrame([{"gene":g,"rna_tss":baseline[g]} for g in GENES]).to_csv(baseline_path,sep="\t",index=False); audit.add_api_calls("3B",1); audit.add_api_requests("3B",1)
    checkpoint=root/"scores.tsv"; done=pd.read_csv(checkpoint,sep="\t") if checkpoint.exists() else pd.DataFrame(); completed=set(done.state_id.astype(str)) if not done.empty else set()
    states=[]
    for offset in range(SCAN_START,SCAN_END+1):
        pos=RS_POS+offset; idx=pos-1-int(interval.start); ref=native[idx]
        for alt in BASES:
            if alt==ref: continue
            state=f"{CHROM}:{pos}:{ref}>{alt}|rs12740374_T_background"
            if state not in completed: states.append((state,offset,pos,ref,alt,native[:idx]+alt+native[idx+1:]))
    parts=[] if done.empty else [done]
    for start in range(0,len(states),batch_size):
        batch=states[start:start+batch_size]
        with audit.step(f"3B: score substitutions {start+1}-{start+len(batch)}"):
            outputs=client.predict_sequences(sequences=[x[5] for x in batch],requested_outputs={dna_client.OutputType.RNA_SEQ},ontology_terms=[LIVER],intervals=[interval]*len(batch),progress_bar=False,max_workers=max_workers)
            rows=[]
            for state,output in zip(batch,outputs,strict=True):
                sid,off,pos,ref,alt,seq=state; values=_summarize(output,interval)
                for gene in GENES: rows.append(dict(state_id=sid,edit_offset=off,edit_pos_hg38=pos,native_base=ref,alt_base=alt,sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),gene=gene,intact_T_rna_tss=baseline[gene],edited_rna_tss=values[gene],delta_vs_intact_T=values[gene]-baseline[gene],loss_vs_intact_T=baseline[gene]-values[gene]))
            parts.append(pd.DataFrame(rows)); all_rows=pd.concat(parts,ignore_index=True).drop_duplicates(["state_id","gene"],keep="last"); all_rows.to_csv(checkpoint,sep="\t",index=False,float_format="%.12g"); parts=[all_rows]
            audit.add_api_calls("3B",len(batch)); audit.add_api_requests("3B",1)
    all_rows=pd.concat(parts,ignore_index=True); mean=all_rows.groupby(["state_id","edit_offset","edit_pos_hg38","native_base","alt_base","sequence_sha256"],as_index=False).loss_vs_intact_T.mean().rename(columns={"loss_vs_intact_T":"loss_vs_intact_T"})
    position=mean.groupby(["edit_offset","edit_pos_hg38","native_base"],as_index=False).loss_vs_intact_T.agg(max_loss="max",min_loss="min",mean_loss="mean"); position["positive_max_loss"]=position.max_loss.clip(lower=0)
    out=run_dir/"derived/Figure3B_native_501bp_ism"; out.mkdir(parents=True,exist_ok=True); mean.to_csv(out/"native_locus_501bp_three_gene_mean_scores.tsv",sep="\t",index=False,float_format="%.12g"); position.to_csv(out/"native_locus_501bp_three_gene_mean_position_summary.tsv",sep="\t",index=False,float_format="%.12g"); _hotspots(position).to_csv(out/"native_locus_501bp_SORT1_hotspots.tsv",sep="\t",index=False,float_format="%.12g")
    # Per-gene (not gene-averaged) scores: the only extra artifact Figure 3C's
    # PWM-compatibility scan needs from 3B (see run_fig3c / REPRODUCIBILITY_NEXT_STEPS.md).
    all_rows.to_csv(out/"native_locus_501bp_all_gene_scores.tsv",sep="\t",index=False,float_format="%.12g")
    _render3b(position,run_dir/"figures/Figure3B.svg")


def _render3b(position: pd.DataFrame, output: Path) -> None:
    fig,ax=plt.subplots(figsize=(6.2,2.2)); ax.plot(position.edit_offset,position.positive_max_loss,color="#2478b5",lw=.8); ax.axvspan(-1,8,color="#d62728",alpha=.15); ax.set(xlabel="Position from rs12740374 (bp)",ylabel="Maximum 3-gene RNA loss"); _save_svg(fig,output)


# --- Figure 3C: native-locus motif-family PWM-compatibility scan ---------
#
# Pure local computation, no AlphaGenome calls: scans JASPAR PWMs across the
# 501 bp native-T locus and combines native-motif compatibility with the
# largest predicted 3-gene RNA loss from Figure 3B (run_fig3b) at each motif
# position, within six prespecified sequence-sensitive windows plus the
# known C/EBP positive control. Ported from the working archive's
# report/figure3_restructured/panel_C_motif_family_disruption/make_panel_C.py
# and sort1_pwm_motif_analysis.py (not part of this repository) -- see
# REPRODUCIBILITY_NEXT_STEPS.md R019.

JASPAR_URL = "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
_FIG3C_START_OFFSET, _FIG3C_END_OFFSET = -250, 250
_FIG3C_MIN_NATIVE_SCORE_FRACTION = 0.70
_FIG3C_N_DISPLAY_FAMILIES = 8
_BASE_TO_IDX = {b: i for i, b in enumerate("ACGT")}
_COMPLEMENT = str.maketrans({"A": "T", "C": "G", "G": "C", "T": "A"})


def _reverse_complement(seq: str) -> str:
    return seq.upper().translate(_COMPLEMENT)[::-1]


class _PWM:
    def __init__(self, name: str, matrix_id: str, log_odds: np.ndarray) -> None:
        self.name, self.matrix_id, self.log_odds = name, matrix_id, log_odds

    @property
    def length(self) -> int:
        return int(self.log_odds.shape[0])

    def score(self, seq: str) -> float:
        if len(seq) != self.length:
            raise ValueError("Sequence length does not match PWM length.")
        total = 0.0
        for i, base in enumerate(seq.upper()):
            idx = _BASE_TO_IDX.get(base)
            if idx is None:
                return float("-inf")
            total += float(self.log_odds[i, idx])
        return total


class _MotifHit:
    def __init__(self, motif: str, matrix_id: str, start: int, end: int, strand: str, score: float, matched_seq: str) -> None:
        self.motif, self.matrix_id, self.start, self.end = motif, matrix_id, start, end
        self.strand, self.score, self.matched_seq = strand, score, matched_seq


def _parse_jaspar_pfm_file(path: Path) -> dict[str, dict[str, object]]:
    import re

    records: dict[str, dict[str, object]] = {}
    current_id: str | None = None
    current_name: str | None = None
    current_pfm: dict[str, list[int]] = {}
    base_pat = re.compile(r"^([ACGT])\s*\[([^\]]+)\]\s*$")
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None and set(current_pfm) == {"A", "C", "G", "T"}:
                    records[current_id] = {"matrix_id": current_id, "name": current_name or current_id, "pfm": current_pfm}
                parts = line[1:].strip().split()
                current_id = parts[0]
                current_name = " ".join(parts[1:]) if len(parts) > 1 else current_id
                current_pfm = {}
                continue
            match = base_pat.match(line)
            if match is None or current_id is None:
                continue
            current_pfm[match.group(1)] = [int(float(x)) for x in match.group(2).split()]
    if current_id is not None and set(current_pfm) == {"A", "C", "G", "T"}:
        records[current_id] = {"matrix_id": current_id, "name": current_name or current_id, "pfm": current_pfm}
    return records


def _pfm_to_pwm(*, motif: str, matrix_id: str, pfm: dict[str, list[int]], pseudocount: float = 0.5) -> _PWM:
    counts = np.array([pfm["A"], pfm["C"], pfm["G"], pfm["T"]], dtype=float).T + float(pseudocount)
    probs = counts / counts.sum(axis=1, keepdims=True)
    log_odds = np.log2(probs / np.full((4,), 0.25, dtype=float))
    return _PWM(motif, matrix_id, log_odds)


def _scan_pwm(seq: str, pwm: _PWM, *, start_pos_1based: int) -> list[_MotifHit]:
    seq = seq.upper()
    hits: list[_MotifHit] = []
    length = pwm.length
    for i in range(0, len(seq) - length + 1):
        window = seq[i : i + length]
        score_fwd = pwm.score(window)
        rc = _reverse_complement(window)
        score_rev = pwm.score(rc)
        strand, score, matched = ("+", score_fwd, window) if score_fwd >= score_rev else ("-", score_rev, rc)
        hits.append(_MotifHit(pwm.name, pwm.matrix_id, start_pos_1based + i, start_pos_1based + i + length - 1, strand, float(score), matched))
    return hits


def _fig3c_score_range(pwm: _PWM) -> tuple[float, float]:
    return float(np.min(pwm.log_odds, axis=1).sum()), float(np.max(pwm.log_odds, axis=1).sum())


def _fig3c_score_fraction(pwm: _PWM, sequence: str) -> float:
    low, high = _fig3c_score_range(pwm)
    return float((pwm.score(sequence) - low) / max(high - low, 1e-9))


def _motif_family(motif: str) -> str:
    name = str(motif).upper()
    if name.startswith("CEBP"): return "C/EBP/bZIP"
    if name.startswith(("JUN","FOS","ATF","CREB","MAF","BACH","NFE2","BATF","HLF","DBP","TEF","NFIL3")): return "AP-1/bZIP"
    if name.startswith("FOX"): return "FOX/forkhead"
    if name.startswith("HNF4"): return "HNF4/nuclear receptor"
    if name.startswith(("NR1","NR2","NR3","RXR","PPAR","RARA","RARB","RARG")): return "nuclear receptor"
    if name.startswith(("MEIS","HOX","BARHL","BARX","LHX","OTX","VAX","BSX","PAX","NANOG","EMX","NKX","RAX","GSX","RHOXF","HMX","DLX","MSX","PITX","PRRX","PROP","PHOX","SIX")): return "homeobox"
    if name.startswith(("ARNT","HIF","MXI","MAX","CLOCK","TFAP4","TFE","USF","MYC","MLX","NPAS","MNT")): return "bHLH"
    if name.startswith(("KLF","SP")): return "KLF/SP"
    if name.startswith(("ETS","ELF","ELK","ERG","FLI","GABP","SPI")): return "ETS"
    if name.startswith("GATA"): return "GATA"
    if name.startswith(("RFX",)): return "RFX"
    if name.startswith(("ZNF","ZKSCAN","ZBTB","PATZ","HIC","BNC")): return "zinc finger"
    if name.startswith(("NFIA","NFIB","NFIC","NFIX")): return "NFI"
    if name.startswith(("TBX",)): return "T-box"
    if name.startswith(("RUNX",)): return "RUNX"
    if name.startswith(("SOX","SRY","TCF7","LEF1")): return "SOX/HMG"
    if name.startswith(("SMAD",)): return "SMAD"
    return str(motif)


def _fig3c_windows(hotspots: pd.DataFrame) -> pd.DataFrame:
    control = pd.DataFrame([{"hotspot": "C/EBP", "start_offset": CEBP_START, "end_offset": CEBP_END}])
    out = pd.concat([hotspots[["hotspot", "start_offset", "end_offset"]], control], ignore_index=True)
    return out.sort_values("start_offset").reset_index(drop=True)


def _fig3c_max_loss_choices(all_gene_scores: pd.DataFrame) -> pd.DataFrame:
    scores = all_gene_scores[all_gene_scores.gene.isin(GENES)].copy()
    grouped = scores.groupby(["edit_offset", "native_base", "alt_base"], as_index=False).agg(
        mean_three_gene_loss=("loss_vs_intact_T", "mean"), n_genes=("gene", "nunique")
    )
    if not grouped["n_genes"].eq(len(GENES)).all():
        raise ValueError("Incomplete three-gene predictions in native-locus ISM table")
    index = grouped.groupby("edit_offset")["mean_three_gene_loss"].idxmax()
    choices = grouped.loc[index, ["edit_offset", "native_base", "alt_base", "mean_three_gene_loss"]].copy()
    return choices.rename(columns={"alt_base": "most_disruptive_alt", "mean_three_gene_loss": "max_mean_three_gene_loss"}).set_index("edit_offset")


def _fig3c_mutated_motif_sequence(native_501: str, *, start_offset: int, end_offset: int, strand: str, choices: pd.DataFrame) -> tuple[str, str]:
    native_genomic = native_501[start_offset - _FIG3C_START_OFFSET : end_offset - _FIG3C_START_OFFSET + 1]
    chars = list(native_genomic)
    for offset in range(start_offset, end_offset + 1):
        if offset not in choices.index:
            continue
        chars[offset - start_offset] = str(choices.loc[offset, "most_disruptive_alt"])
    edited_genomic = "".join(chars)
    if strand == "-":
        return _reverse_complement(native_genomic), _reverse_complement(edited_genomic)
    return native_genomic, edited_genomic


def _fig3c_scan_and_score(native: str, jaspar_path: Path, all_gene_scores: pd.DataFrame, hotspots: pd.DataFrame) -> pd.DataFrame:
    choices = _fig3c_max_loss_choices(all_gene_scores)
    candidates = _fig3c_windows(hotspots)
    records = _parse_jaspar_pfm_file(jaspar_path)
    rows: list[dict[str, object]] = []
    for matrix_id, record in records.items():
        pwm = _pfm_to_pwm(motif=str(record["name"]), matrix_id=str(matrix_id), pfm=record["pfm"])  # type: ignore[arg-type]
        for hit in _scan_pwm(native, pwm, start_pos_1based=RS_POS + _FIG3C_START_OFFSET):
            start_offset, end_offset = hit.start - RS_POS, hit.end - RS_POS
            for candidate in candidates.itertuples(index=False):
                window_start, window_end = int(candidate.start_offset), int(candidate.end_offset)
                if start_offset > window_end or end_offset < window_start:
                    continue
                native_motif, edited_motif = _fig3c_mutated_motif_sequence(
                    native, start_offset=start_offset, end_offset=end_offset, strand=hit.strand, choices=choices
                )
                native_fraction = _fig3c_score_fraction(pwm, native_motif)
                if native_fraction < _FIG3C_MIN_NATIVE_SCORE_FRACTION:
                    continue
                edited_fraction = _fig3c_score_fraction(pwm, edited_motif)
                overlap = max(0, min(end_offset, window_end) - max(start_offset, window_start) + 1)
                overlap_fraction = overlap / max(1, end_offset - start_offset + 1)
                disruption = max(0.0, native_fraction - edited_fraction)
                rows.append({
                    "hotspot": str(candidate.hotspot), "window_start_offset": window_start, "window_end_offset": window_end,
                    "motif_family": _motif_family(hit.motif), "motif": hit.motif, "matrix_id": hit.matrix_id,
                    "motif_start_offset": start_offset, "motif_end_offset": end_offset, "strand": hit.strand,
                    "native_motif_sequence": native_motif, "max_loss_motif_sequence": edited_motif,
                    "native_pwm_score_fraction": native_fraction, "edited_pwm_score_fraction": edited_fraction,
                    "pwm_score_loss": disruption, "overlap_fraction": overlap_fraction,
                    "pwm_disruption_score": native_fraction * disruption * overlap_fraction,
                })
    exact = pd.DataFrame(rows)
    if exact.empty:
        raise RuntimeError("No qualifying JASPAR hits found")
    return exact.sort_values(["hotspot", "pwm_disruption_score"], ascending=[True, False])


def _fig3c_collapse(exact: pd.DataFrame, hotspots: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    best = exact.sort_values("pwm_disruption_score", ascending=False).groupby(["hotspot", "motif_family"], as_index=False).first()
    order = _fig3c_windows(hotspots)["hotspot"].tolist()
    heat = best.pivot(index="hotspot", columns="motif_family", values="pwm_disruption_score").fillna(0.0).reindex(order).fillna(0.0)
    family_rank = heat.max(axis=0).sort_values(ascending=False)
    selected = family_rank.head(_FIG3C_N_DISPLAY_FAMILIES).index.tolist()
    if "C/EBP/bZIP" in heat.columns and "C/EBP/bZIP" not in selected:
        selected = selected[: _FIG3C_N_DISPLAY_FAMILIES - 1] + ["C/EBP/bZIP"]
    selected = sorted(selected, key=lambda family: (family != "C/EBP/bZIP", -float(family_rank.get(family, 0.0))))
    return heat[selected], best


def _render3c(heat: pd.DataFrame, output: Path) -> None:
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica Neue", "DejaVu Sans"],
        "font.size": 9.3, "axes.labelsize": 9.3, "xtick.labelsize": 8.8, "ytick.labelsize": 10.2,
        "svg.fonttype": "none", "pdf.fonttype": 42,
    })
    fig, ax = plt.subplots(figsize=(3.20, 2.78))
    fig.subplots_adjust(left=0.17, right=0.81, bottom=0.43, top=0.96)
    values = heat.to_numpy(dtype=float)
    vmax = max(float(np.nanmax(values)), 1e-9)
    image = ax.pcolormesh(np.arange(values.shape[1] + 1) - 0.5, np.arange(values.shape[0] + 1) - 0.5, values, cmap="viridis", vmin=0, vmax=vmax, linewidth=0, edgecolors="none", antialiased=False, rasterized=False)
    ax.set_xlim(-0.5, values.shape[1] - 0.5); ax.set_ylim(values.shape[0] - 0.5, -0.5)
    ax.set_yticks(np.arange(len(heat.index))); ax.set_yticklabels(heat.index, fontweight="bold")
    ax.set_xticks(np.arange(len(heat.columns))); ax.set_xticklabels(heat.columns, rotation=90, ha="center", va="top")
    ax.tick_params(width=0.7, length=2.5, pad=2)
    for spine in ax.spines.values(): spine.set_linewidth(0.7)
    cbar = fig.colorbar(image, ax=ax, fraction=0.055, pad=0.025)
    if cbar.solids is not None: cbar.solids.set_rasterized(False)
    cbar.set_label("PWM disruption score", labelpad=4); cbar.ax.tick_params(width=0.7, length=2.5, pad=2); cbar.outline.set_linewidth(0.7)
    _save_svg(fig, output)


def run_fig3c(run_dir: Path, audit: Audit, fasta_path: Path) -> None:
    """Figure 3C: motif-family PWM-compatibility scan. Requires 3B's derived
    outputs (native_locus_501bp_all_gene_scores.tsv, ..._SORT1_hotspots.tsv)."""
    genome, _, _, _ = _ag()
    scores_path = run_dir / "derived/Figure3B_native_501bp_ism/native_locus_501bp_all_gene_scores.tsv"
    hotspots_path = run_dir / "derived/Figure3B_native_501bp_ism/native_locus_501bp_SORT1_hotspots.tsv"
    if not scores_path.exists() or not hotspots_path.exists():
        raise RuntimeError("Figure 3C requires Figure 3B's derived outputs; run panel 3B first (same --run-dir).")
    jaspar_path = run_dir / "raw" / "JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
    with audit.step("3C: download JASPAR 2024 CORE PFM matrices"):
        info = download(JASPAR_URL, jaspar_path)
        if not any(entry.get("path") == info["path"] for entry in audit.downloads):
            audit.downloads.append(info); audit.save()
    with audit.step("3C: scan JASPAR motifs and score PWM compatibility vs Figure 3B ISM"):
        interval = genome.Interval(CHROM, RS_POS, RS_POS).resize(501)
        native = _native_t(fasta_path, interval)
        all_gene_scores = pd.read_csv(scores_path, sep="\t")
        hotspots = pd.read_csv(hotspots_path, sep="\t")
        exact = _fig3c_scan_and_score(native, jaspar_path, all_gene_scores, hotspots)
        heat, best = _fig3c_collapse(exact, hotspots)
        out = run_dir / "derived/Figure3C_pwm_compatibility"; out.mkdir(parents=True, exist_ok=True)
        exact.to_csv(out / "Figure3C_PWM_disruption_exact_hits.tsv", sep="\t", index=False)
        best.to_csv(out / "Figure3C_PWM_disruption_best_family_hits.tsv", sep="\t", index=False)
        heat.to_csv(out / "Figure3C_PWM_disruption_values.tsv", sep="\t")
        _render3c(heat, run_dir / "figures/Figure3C.svg")


# --- Figures 3E/3F/3G: shared scramble-design / REF-ALT retention machinery
#
# All three panels build many REF/ALT sequence pairs around rs12740374 under
# a composition-preserving scramble background, score each with AlphaGenome,
# and compute a "retention vs an unscrambled/native baseline" ratio.  Ported
# from the working archive's report/panel_asymmetric_scramble/
# run_asymmetric_native_scramble.py (the shared "asym" library imported by
# both run_single_arm_recovery.py [3E] and run_component_necessity.py [3G])
# and the "wide_main_panel" 1bp-grid driver [3F] -- not part of this
# repository. See REPRODUCIBILITY_NEXT_STEPS.md R019.

_MINOR_MOTIF = "GTTGCTCAAT"
_MOTIF_START_1BASED = RS_POS - 1


def _motif_at_native_position(seq: str, seq_start0: int) -> str:
    start_index = _MOTIF_START_1BASED - 1 - seq_start0
    return seq[start_index : start_index + len(_MINOR_MOTIF)]


def _set_base(seq: str, seq_start0: int, pos1: int, base_value: str) -> str:
    index = int(pos1) - 1 - int(seq_start0)
    chars = list(seq)
    chars[index] = str(base_value).upper()
    return "".join(chars)


def _gene_tss_table() -> dict[str, int]:
    """SORT1/PSRC1/CELSR2 TSS positions -- identical values to what GENCODE
    v46 lookup produces (verified against sort1_hypothesis_panels.py's
    _load_gene_tss), sourced from the already-available GENE_TABLE so 3E/3F/
    3G need no GENCODE download at all."""
    return dict(zip(GENE_TABLE.gene_symbol, GENE_TABLE.tss_hg38.astype(int)))


def _score_design(
    run_dir: Path, audit: Audit, panel: str, interval: Any, sequences: list[str],
    design: pd.DataFrame, *, batch_size: int, max_workers: int,
) -> pd.DataFrame:
    """Score every unique sequence referenced by `design` for all three genes'
    liver RNA(TSS) signal, checkpointed per sequence hash (same pattern as
    run_fig3b). `design` must have `sequence_index` and `sequence_sha256`
    columns; returns one row per (design row, gene)."""
    genome, dna_client, dna_model, _ = _ag()
    cache = run_dir / "predictions" / f"Figure{panel}_sequence_cache"
    cache.mkdir(parents=True, exist_ok=True)
    hash_to_seq = dict(zip(design.sequence_sha256, [sequences[int(i)] for i in design.sequence_index]))
    values: dict[str, dict[str, float]] = {}
    for digest in hash_to_seq:
        path = cache / f"{digest}.tsv"
        if path.exists():
            values[digest] = pd.read_csv(path, sep="\t").set_index("gene").liver_rna_signal.to_dict()
    missing = [digest for digest in hash_to_seq if digest not in values]
    client = None
    for start in range(0, len(missing), batch_size):
        batch = missing[start : start + batch_size]
        if client is None:
            client = dna_client.create(api_key(), model_version=dna_model.ModelVersion.ALL_FOLDS, timeout=300)
        with audit.step(f"{panel}: score sequences {start + 1}-{start + len(batch)}"):
            outputs = client.predict_sequences(
                sequences=[hash_to_seq[digest] for digest in batch],
                requested_outputs={dna_client.OutputType.RNA_SEQ}, ontology_terms=[LIVER],
                intervals=[interval] * len(batch), progress_bar=False, max_workers=max_workers,
            )
            for digest, output in zip(batch, outputs, strict=True):
                gene_values = _summarize(output, interval)
                pd.DataFrame([{"gene": g, "liver_rna_signal": gene_values[g]} for g in GENES]).to_csv(
                    cache / f"{digest}.tsv", sep="\t", index=False
                )
                values[digest] = gene_values
            audit.add_api_calls(panel, len(batch))
            audit.add_api_requests(panel, 1)
    rows: list[dict[str, object]] = []
    for row in design.itertuples(index=False):
        gene_values = values[row.sequence_sha256]
        for gene in GENES:
            rows.append({**row._asdict(), "gene": gene, "liver_rna_signal": gene_values[gene]})
    return pd.DataFrame(rows)


# --- Figure 3G: expanded component-necessity audit ------------------------

_FIG3G_SEEDS = list(range(1740374, 1740382))
_FIG3G_COMPONENTS: list[tuple[str, list[tuple[int, int]], str]] = [
    ("native", [], "control"), ("full_315", [(-179, 135)], "control"),
    ("H1", [(-178, -167)], "hotspot"), ("H2", [(-152, -141)], "hotspot"),
    ("H3", [(-87, -76)], "hotspot"), ("H4", [(-61, -50)], "hotspot"),
    ("H5", [(-37, -26)], "hotspot"), ("H6", [(50, 61)], "hotspot"),
    ("C/EBP", [(-1, 8)], "positive_control"),
    ("upstream_arm", [(-179, -2)], "arm"), ("downstream_arm", [(9, 135)], "arm"),
    ("both_arms", [(-179, -2), (9, 135)], "arm"),
]


def _shuffle_groups(sequence: str, *, seq_start0: int, groups: list[list[int]], seed: int) -> tuple[str, int]:
    chars = list(sequence)
    changed: set[int] = set()
    for group_index, positions in enumerate(groups):
        indices = [position - 1 - seq_start0 for position in positions]
        original = [chars[index] for index in indices]
        shuffled = original.copy()
        for attempt in range(100):
            shuffled = original.copy()
            rng = np.random.default_rng(seed + 10_007 * (group_index + 1) + 1_000_000_007 * attempt)
            rng.shuffle(shuffled)
            if any(old != new for old, new in zip(original, shuffled, strict=True)):
                break
        else:
            raise ValueError("Component cannot be changed by composition-preserving shuffling.")
        for index, old, new in zip(indices, original, shuffled, strict=True):
            chars[index] = new
            if new != old:
                changed.add(index)
    return "".join(chars), len(changed)


def _fig3g_positions_for_intervals(intervals: list[tuple[int, int]]) -> list[list[int]]:
    return [[RS_POS + offset for offset in range(start, end + 1) if offset != 0] for start, end in intervals]


def _fig3g_build_design(ref_seq: str, seq_start0: int) -> tuple[list[str], pd.DataFrame]:
    sequences: list[str] = []
    rows: list[dict[str, object]] = []
    for seed in _FIG3G_SEEDS:
        for component_index, (component, intervals, family) in enumerate(_FIG3G_COMPONENTS):
            if component == "native":
                background, changed_positions = ref_seq, 0
            else:
                background, changed_positions = _shuffle_groups(
                    ref_seq, seq_start0=seq_start0,
                    groups=_fig3g_positions_for_intervals(intervals),
                    seed=seed + 1_000_003 * (component_index + 1),
                )
            pair: dict[str, str] = {}
            for allele, base_value in (("REF", "G"), ("ALT", "T")):
                sequence = _set_base(background, seq_start0, RS_POS, base_value)
                pair[allele] = sequence
                sequences.append(sequence)
                rows.append({
                    "sequence_index": len(sequences) - 1, "scramble_seed": seed, "component": component,
                    "component_family": family, "intervals_relative": ";".join(f"{a}:{b}" for a, b in intervals),
                    "allele": allele, "changed_positions_vs_native": changed_positions,
                    "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                    "motif_sequence_at_native_position": _motif_at_native_position(sequence, seq_start0),
                })
            differences = [i for i, (a, b) in enumerate(zip(pair["REF"], pair["ALT"], strict=True)) if a != b]
            rs_index = RS_POS - 1 - seq_start0
            if differences != [rs_index]:
                raise AssertionError(f"{component}/seed{seed}: allele pair is not matched.")
    return sequences, pd.DataFrame(rows)


def _fig3g_retention(scored: pd.DataFrame) -> pd.DataFrame:
    keys = ["scramble_seed", "component", "component_family", "gene"]
    ref = scored[scored.allele.eq("REF")].copy()
    alt = scored[scored.allele.eq("ALT")].copy()
    paired = ref.merge(alt, on=keys, suffixes=("_ref", "_alt"), validate="one_to_one")
    paired["delta_liver"] = paired.liver_rna_signal_alt - paired.liver_rna_signal_ref
    native = paired[paired.component.eq("native")][["scramble_seed", "gene", "delta_liver"]].rename(
        columns={"delta_liver": "native_delta_liver"}
    )
    paired = paired.merge(native, on=["scramble_seed", "gene"], validate="many_to_one")
    paired["retention"] = paired.delta_liver.abs() / paired.native_delta_liver.abs()
    paired["signed_retention"] = paired.delta_liver / paired.native_delta_liver
    return paired


def _fetch_reference(fasta_path: Path, interval: Any) -> str:
    """True (unsubstituted) reference sequence for `interval`, asserting the
    expected G allele at rs12740374 -- `_fig3g_build_design`/`_fig3e_build_design`
    substitute REF/ALT explicitly, unlike `_native_t` (used by 3B/3C)."""
    import pysam

    fasta = pysam.FastaFile(str(fasta_path))
    sequence = fasta.fetch(CHROM, int(interval.start), int(interval.end)).upper()
    index = RS_POS - 1 - int(interval.start)
    if sequence[index] != "G":
        raise ValueError(f"Reference allele mismatch at rs12740374: {sequence[index]} != G")
    return sequence


def run_fig3g(run_dir: Path, audit: Audit, fasta_path: Path, *, batch_size: int = 32, max_workers: int = 4) -> None:
    genome, _, _, _ = _ag()
    interval = genome.Interval(CHROM, RS_POS, RS_POS).resize(SEQ_LEN)
    ref_seq = _fetch_reference(fasta_path, interval)
    with audit.step("3G: build component-necessity scramble design"):
        sequences, design = _fig3g_build_design(ref_seq, int(interval.start))
        out = run_dir / "derived/Figure3G_component_necessity"; out.mkdir(parents=True, exist_ok=True)
        design.to_csv(out / "sequence_design.csv", index=False)
    scored = _score_design(run_dir, audit, "3G", interval, sequences, design, batch_size=batch_size, max_workers=max_workers)
    with audit.step("3G: compute retention and render"):
        retention = _fig3g_retention(scored)
        retention.to_csv(out / "retention_by_seed.csv", index=False)
        summary = retention.groupby(["component", "component_family", "gene"], as_index=False).agg(
            n_seeds=("scramble_seed", "nunique"), mean_retention=("retention", "mean"),
            median_retention=("retention", "median"), sd_retention=("retention", "std"),
            mean_signed_retention=("signed_retention", "mean"),
            sign_preserved_fraction=("signed_retention", lambda x: float(np.mean(np.asarray(x) > 0))),
        )
        summary.to_csv(out / "summary.csv", index=False)
        three_gene_mean = _fig3g_three_gene_mean(summary, retention)
        three_gene_mean.to_csv(out / "Figure3G_component_necessity_three_gene_mean_source.tsv", sep="\t", index=False)
        _render3g(three_gene_mean, run_dir / "figures/Figure3G.svg")


def _fig3g_three_gene_mean(summary: pd.DataFrame, retention: pd.DataFrame) -> pd.DataFrame:
    """Collapse the per-gene summary to the 3-gene-mean table Figure 3G plots
    and the release's outputs/source_data/Figure3G_component_necessity.tsv
    records. SEM is computed across the eight seed-level 3-gene means, not by
    averaging the three gene-specific SEMs (matches the working archive's
    make_panel.py)."""
    order = [name for name, _, _ in _FIG3G_COMPONENTS]
    frame = summary[summary.gene.isin(GENES)].copy()
    data = frame.groupby("component", as_index=True).agg(
        mean_retention=("mean_retention", "mean"), median_retention=("median_retention", "mean"), n_seeds=("n_seeds", "min")
    )
    per_seed = retention[retention.gene.isin(GENES)].groupby(["component", "scramble_seed"], as_index=False).agg(
        retention=("retention", "mean"), n_genes=("gene", "nunique")
    )
    sem_by_component = per_seed.groupby("component")["retention"].agg(lambda v: v.std(ddof=1) / np.sqrt(v.count()))
    # The true median-of-3-gene-means per seed, not the mean of the three
    # genes' individual medians (`data["median_retention"]` above) -- matches
    # the working archive's make_panel.py exactly.
    median_by_component = per_seed.groupby("component")["retention"].median()
    rows = [
        {"model": "ALL_FOLDS", "component": component, "mean_retention": data.loc[component, "mean_retention"],
         "sem_retention": sem_by_component.loc[component], "median_retention": median_by_component.loc[component],
         "n_seeds": int(data.loc[component, "n_seeds"])}
        for component in order
    ]
    return pd.DataFrame(rows)


def _render3g(three_gene_mean: pd.DataFrame, output: Path) -> None:
    order = [name for name, _, _ in _FIG3G_COMPONENTS]
    labels = {"native": "Intact", "full_315": "Full 315", "upstream_arm": "Up", "downstream_arm": "Down", "both_arms": "Both"}
    data = three_gene_mean.set_index("component")
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"],
        "font.size": 8.5, "axes.labelsize": 9, "xtick.labelsize": 7.2, "ytick.labelsize": 7.5, "axes.linewidth": 0.8,
    })
    x = np.arange(len(order), dtype=float)
    fig, ax = plt.subplots(figsize=(87 / 25.4, 54 / 25.4))
    fig.subplots_adjust(left=0.20, right=0.99, bottom=0.43, top=0.74)
    means = np.array([data.loc[item, "mean_retention"] for item in order])
    sems = np.array([data.loc[item, "sem_retention"] for item in order])
    ax.errorbar(x, means, yerr=sems, fmt="o", color="#2478b5", markersize=4.4, linewidth=0, elinewidth=0.9, capsize=1.8, zorder=4)
    ax.axhline(1, color="#777777", linestyle=":", linewidth=0.9); ax.axhline(0, color="#bbbbbb", linewidth=0.65)
    ax.set_xticks(x, [labels.get(item, item) for item in order])
    ax.tick_params(axis="x", labelsize=6.5, pad=7)
    plt.setp(ax.get_xticklabels(), rotation=90, ha="center", va="top", rotation_mode="anchor")
    ax.set_ylabel("Mean RNA retention\n(3 genes)", fontsize=8.5, labelpad=3)
    ax.set_ylim(-0.04, 1.10); ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#dddddd", linewidth=0.5); ax.set_axisbelow(True)
    _save_svg(fig, output)


# --- Figure 3E: directional single-arm motif-protected recovery -----------
#
# Ported from the working archive's report/panel_asymmetric_scramble/
# run_single_arm_recovery.py (not part of this repository).

_FIG3E_UPSTREAM_EXTENTS = [0, 2, 5, 10, 20, 40, 80, 120, 160, 220, 300, 400, 600, 800, 1000]
_FIG3E_DOWNSTREAM_EXTENTS = [0, 10, 20, 40, 80, 120, 160, 220, 300, 400, 600, 800, 1000]
_FIG3E_SEEDS = list(range(1740374, 1740382))
_FIG3E_LOCAL_HALF_BP = 1000


def _scramble_range(chars: list[str], *, seq_start0: int, start1: int, end1: int, seed: int) -> None:
    if int(end1) < int(start1):
        return
    indices = [pos1 - 1 - int(seq_start0) for pos1 in range(int(start1), int(end1) + 1)]
    usable = [i for i in indices if 0 <= i < len(chars) and chars[i].upper() in {"A", "C", "G", "T"}]
    if len(usable) <= 1:
        return
    values = [chars[i] for i in usable]
    rng = np.random.default_rng(int(seed))
    rng.shuffle(values)
    for index, value in zip(usable, values, strict=True):
        chars[index] = value


def _fig3e_scramble_template(ref_seq: str, *, seq_start0: int, local_half_bp: int, scramble_seed: int) -> str:
    """Shuffle each non-motif arm independently, never touching the
    protected -1..+8 C/EBP motif interval."""
    chars = list(ref_seq)
    seed_base = int(scramble_seed) + 7_919
    _scramble_range(chars, seq_start0=seq_start0, start1=RS_POS - local_half_bp, end1=RS_POS - 2, seed=seed_base + 11)
    _scramble_range(chars, seq_start0=seq_start0, start1=RS_POS + 9, end1=RS_POS + local_half_bp, seed=seed_base + 29)
    scrambled = "".join(chars)
    for start1, end1 in ((RS_POS - local_half_bp, RS_POS - 2), (RS_POS + 9, RS_POS + local_half_bp)):
        native_arm = ref_seq[start1 - 1 - seq_start0 : end1 - seq_start0]
        scrambled_arm = scrambled[start1 - 1 - seq_start0 : end1 - seq_start0]
        if sorted(native_arm) != sorted(scrambled_arm):
            raise AssertionError(f"Arm composition changed for {start1}:{end1}, seed {scramble_seed}.")
    native_motif = ref_seq[RS_POS - 2 - seq_start0 : RS_POS + 8 - seq_start0]
    scrambled_motif = scrambled[RS_POS - 2 - seq_start0 : RS_POS + 8 - seq_start0]
    if native_motif != scrambled_motif:
        raise AssertionError("The protected motif changed during non-motif shuffling.")
    return scrambled


def _fig3e_copy_native_positions(background: str, native: str, *, seq_start0: int, intervals: list[tuple[int, int]]) -> str:
    chars = list(background)
    for start1, end1 in intervals:
        if end1 < start1:
            continue
        for pos1 in range(int(start1), int(end1) + 1):
            chars[pos1 - 1 - seq_start0] = native[pos1 - 1 - seq_start0]
    return "".join(chars)


def _fig3e_build_design(ref_seq: str, seq_start0: int) -> tuple[list[str], pd.DataFrame]:
    sequences: list[str] = []
    rows: list[dict[str, object]] = []

    def add_pair(background: str, *, seed: int, condition: str, arm: str, extent_bp: int) -> None:
        pair: dict[str, str] = {}
        for allele, base_value in (("REF", "G"), ("ALT", "T")):
            sequence = _set_base(background, seq_start0, RS_POS, base_value)
            motif = _motif_at_native_position(sequence, seq_start0)
            pair[allele] = sequence
            sequences.append(sequence)
            rows.append({
                "sequence_index": len(sequences) - 1, "scramble_seed": int(seed), "condition": condition, "arm": arm,
                "extent_bp": int(extent_bp), "allele": allele, "motif_sequence_at_native_position": motif,
                "exact_minor_motif_at_native_position": motif == _MINOR_MOTIF,
                "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            })
        differences = [i for i, (a, b) in enumerate(zip(pair["REF"], pair["ALT"], strict=True)) if a != b]
        rs_index = RS_POS - 1 - seq_start0
        if differences != [rs_index]:
            raise AssertionError(f"{condition}/{arm}/{extent_bp}/seed{seed}: REF and ALT differ beyond rs12740374.")

    for seed in _FIG3E_SEEDS:
        motif_protected = _fig3e_scramble_template(ref_seq, seq_start0=seq_start0, local_half_bp=_FIG3E_LOCAL_HALF_BP, scramble_seed=seed)
        add_pair(ref_seq, seed=seed, condition="native", arm="none", extent_bp=-1)
        add_pair(motif_protected, seed=seed, condition="motif_only_baseline", arm="shared", extent_bp=0)
        for upstream in _FIG3E_UPSTREAM_EXTENTS:
            if upstream == 0:
                continue
            background = _fig3e_copy_native_positions(motif_protected, ref_seq, seq_start0=seq_start0, intervals=[(RS_POS - upstream, RS_POS - 2)])
            add_pair(background, seed=seed, condition="upstream_recovery", arm="upstream", extent_bp=upstream)
        for downstream in _FIG3E_DOWNSTREAM_EXTENTS:
            if downstream == 0:
                continue
            background = _fig3e_copy_native_positions(motif_protected, ref_seq, seq_start0=seq_start0, intervals=[(RS_POS + 9, RS_POS + downstream)])
            add_pair(background, seed=seed, condition="downstream_recovery", arm="downstream", extent_bp=downstream)

    design = pd.DataFrame(rows)
    alt_non_native = design[design.allele.eq("ALT") & ~design.condition.eq("native")]
    if not alt_non_native.exact_minor_motif_at_native_position.all():
        failed = alt_non_native[~alt_non_native.exact_minor_motif_at_native_position]
        raise AssertionError(f"The complete minor-allele C/EBP motif was not retained in {len(failed)} ALT constructs.")
    return sequences, design


def _fig3e_retention(scored: pd.DataFrame) -> pd.DataFrame:
    keys = ["scramble_seed", "condition", "arm", "extent_bp", "gene"]
    ref = scored[scored.allele.eq("REF")].copy()
    alt = scored[scored.allele.eq("ALT")].copy()
    paired = ref.merge(alt, on=keys, suffixes=("_ref", "_alt"), validate="one_to_one")
    paired["delta_liver"] = paired.liver_rna_signal_alt - paired.liver_rna_signal_ref
    native = paired[paired.condition.eq("native")][["scramble_seed", "gene", "delta_liver"]].rename(
        columns={"delta_liver": "native_delta_liver"}
    )
    paired = paired.merge(native, on=["scramble_seed", "gene"], validate="many_to_one")
    denominator = paired.native_delta_liver.abs().where(paired.native_delta_liver.abs() > 1e-12, np.nan)
    paired["retention_vs_native"] = paired.delta_liver.abs() / denominator
    paired["signed_retention_vs_native"] = paired.delta_liver / paired.native_delta_liver.where(paired.native_delta_liver.abs() > 1e-12, np.nan)
    return paired


def _fig3e_three_gene_source(retention: pd.DataFrame) -> pd.DataFrame:
    """The Figure 3E plotted/compact table: per-arm, per-series (three genes
    plus their mean) retention curves vs recovered extent. Matches the
    working archive's make_panel_E.py exactly, including its zero-extent
    "motif only" baseline point shared by both arms."""
    frame = retention[retention.gene.isin(GENES)].copy()
    arm_specs = {"upstream_recovery": "Upstream", "downstream_recovery": "Downstream"}
    rows: list[dict[str, object]] = []
    for condition, arm_label in arm_specs.items():
        curve = frame[frame.condition.eq(condition)].copy()
        baseline = frame[frame.condition.eq("motif_only_baseline")].copy()
        baseline = baseline.assign(extent_bp=0)
        curve = pd.concat([baseline, curve], ignore_index=True)
        gene_summary = curve.groupby(["gene", "extent_bp"], as_index=False).agg(mean=("retention_vs_native", "mean"))
        for gene in GENES:
            for row in gene_summary[gene_summary.gene.eq(gene)].sort_values("extent_bp").itertuples(index=False):
                rows.append({"arm": arm_label, "series": gene, "extent_bp": row.extent_bp, "mean_retention": row.mean})
        per_seed = curve.groupby(["scramble_seed", "extent_bp"], as_index=False).agg(
            retention=("retention_vs_native", "mean"), n_genes=("gene", "nunique")
        )
        if not per_seed.n_genes.eq(len(GENES)).all():
            raise ValueError("Incomplete gene triplet in directional-recovery source data")
        mean_summary = per_seed.groupby("extent_bp", as_index=False).agg(mean=("retention", "mean"), sd=("retention", "std"), n=("retention", "count"))
        mean_summary["sem"] = mean_summary["sd"] / np.sqrt(mean_summary["n"])
        for row in mean_summary.sort_values("extent_bp").itertuples(index=False):
            rows.append({"arm": arm_label, "series": "3-gene mean", "extent_bp": row.extent_bp, "mean_retention": row.mean, "sem_retention": row.sem})
    return pd.DataFrame(rows)


def _render3e(three_gene_source: pd.DataFrame, output: Path) -> None:
    specs = {"Upstream": ("#2478b5", "-"), "Downstream": ("#d62728", "--")}
    fig, ax = plt.subplots(figsize=(65 / 25.4, 45 / 25.4))
    fig.subplots_adjust(left=0.19, right=0.98, bottom=0.22, top=0.97)
    for arm, (color, linestyle) in specs.items():
        arm_rows = three_gene_source[three_gene_source.arm.eq(arm)]
        for gene in GENES:
            gene_rows = arm_rows[arm_rows.series.eq(gene)].sort_values("extent_bp")
            ax.plot(gene_rows.extent_bp, gene_rows.mean_retention, color=color, linestyle=linestyle, linewidth=0.7, alpha=0.35, zorder=2)
        mean_rows = arm_rows[arm_rows.series.eq("3-gene mean")].sort_values("extent_bp")
        ax.plot(mean_rows.extent_bp, mean_rows.mean_retention, color=color, linestyle=linestyle, linewidth=2.2, label=f"{arm}, 3-gene mean", zorder=4)
        ax.fill_between(mean_rows.extent_bp, mean_rows.mean_retention - mean_rows.sem_retention, mean_rows.mean_retention + mean_rows.sem_retention, color=color, alpha=0.11, linewidth=0)
    ax.axhline(1.0, color="#777777", linewidth=0.85, linestyle=":"); ax.axhline(0.0, color="#aaaaaa", linewidth=0.7)
    ax.set_xlim(0, 1000); ax.set_ylim(-0.04, 1.16)
    ax.set_xticks([0, 200, 400, 600, 800, 1000]); ax.set_xticklabels(["0\nmotif only", "200", "400", "600", "800", "1000"])
    ax.set_xlabel("Native sequence restored (bp)"); ax.set_ylabel("Mean liver RNA retention\nvs intact locus")
    ax.grid(axis="y", color="#dddddd", linewidth=0.5); ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    _save_svg(fig, output)


def run_fig3e(run_dir: Path, audit: Audit, fasta_path: Path, *, batch_size: int = 32, max_workers: int = 4) -> None:
    genome, _, _, _ = _ag()
    interval = genome.Interval(CHROM, RS_POS, RS_POS).resize(SEQ_LEN)
    ref_seq = _fetch_reference(fasta_path, interval)
    with audit.step("3E: build directional single-arm recovery design"):
        sequences, design = _fig3e_build_design(ref_seq, int(interval.start))
        out = run_dir / "derived/Figure3E_directional_recovery"; out.mkdir(parents=True, exist_ok=True)
        design.to_csv(out / "sequence_design.csv", index=False)
    scored = _score_design(run_dir, audit, "3E", interval, sequences, design, batch_size=batch_size, max_workers=max_workers)
    with audit.step("3E: compute retention and render"):
        retention = _fig3e_retention(scored)
        retention.to_csv(out / "retention_by_seed.csv", index=False)
        three_gene_source = _fig3e_three_gene_source(retention)
        three_gene_source.to_csv(out / "Figure3E_directional_scramble_recovery_three_gene_source.tsv", sep="\t", index=False)
        _render3e(three_gene_source, run_dir / "figures/Figure3E.svg")
