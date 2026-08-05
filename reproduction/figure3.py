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

from .common import Audit, api_key
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
    out=run_dir/"derived/Figure3B_native_501bp_ism"; out.mkdir(parents=True,exist_ok=True); mean.to_csv(out/"native_locus_501bp_three_gene_mean_scores.tsv",sep="\t",index=False,float_format="%.12g"); position.to_csv(out/"native_locus_501bp_three_gene_mean_position_summary.tsv",sep="\t",index=False,float_format="%.12g"); _hotspots(position).to_csv(out/"native_locus_501bp_SORT1_hotspots.tsv",sep="\t",index=False,float_format="%.12g"); _render3b(position,run_dir/"figures/Figure3B.svg")


def _render3b(position: pd.DataFrame, output: Path) -> None:
    fig,ax=plt.subplots(figsize=(6.2,2.2)); ax.plot(position.edit_offset,position.positive_max_loss,color="#2478b5",lw=.8); ax.axvspan(-1,8,color="#d62728",alpha=.15); ax.set(xlabel="Position from rs12740374 (bp)",ylabel="Maximum 3-gene RNA loss"); _save_svg(fig,output)
