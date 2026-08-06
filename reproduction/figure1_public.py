"""Public-data extensions for Figure 1C, 1D, and 1F."""

from __future__ import annotations

import gzip
import shutil
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .common import Audit, download, sha256_file
from .figure1 import GENES, VARIANT_POS, VARIANT_RSID, _save_svg

GTEX_URL = "https://storage.googleapis.com/adult-gtex/bulk-qtl/v7/single-tissue-cis-qtl/all_snp_gene_associations/Liver.allpairs.txt.gz"
GTEX_SHA256 = "372afa081939868407afef0322638b8a191e9f05b1451b6a873fdd04e5d81f67"
VCF_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
PANEL_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
PANEL_SHA256 = "b4023dc6ee2d62ee89c8d4d347db4d348e65518d66d346574cdae7a4bbd76858"
HIC_URL = "https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/25104375-a588-46e6-a382-663cee6c332f/4DNFICSTCJQZ.hic"
GENE_IDS = {"SORT1": "ENSG00000134243.7", "CELSR2": "ENSG00000143126.7", "PSRC1": "ENSG00000134222.12"}
TSS = {"CELSR2": 109249538, "PSRC1": 109283186, "SORT1": 109397918}


def _supplied(source: Path, destination: Path, expected: str | None, audit: Audit) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    reused = destination.exists()
    if not reused:
        shutil.copy2(source, destination)
    digest = sha256_file(destination)
    if expected and digest != expected:
        raise ValueError(f"Checksum mismatch for {source}: {digest}")
    audit.downloads.append({"url": f"supplied:{source}", "path": str(destination), "bytes": destination.stat().st_size, "sha256": digest, "reused": reused})
    audit.save()
    return destination


def prepare_figure1_public_inputs(run_dir: Path, audit: Audit, *, gtex_file: Path | None, vcf_file: Path | None, panel_file: Path | None) -> dict[str, Path]:
    raw = run_dir / "raw"
    gtex = _supplied(gtex_file, raw / "Liver.allpairs.txt.gz", GTEX_SHA256, audit) if gtex_file else Path(download(GTEX_URL, raw / "Liver.allpairs.txt.gz")["path"])
    if not gtex_file:
        if sha256_file(gtex) != GTEX_SHA256: raise ValueError("GTEx checksum mismatch")
        audit.downloads.append({"url": GTEX_URL, "path": str(gtex), "bytes": gtex.stat().st_size, "sha256": GTEX_SHA256, "reused": False}); audit.save()
    panel = _supplied(panel_file, raw / "integrated_call_samples_v3.20130502.ALL.panel", PANEL_SHA256, audit) if panel_file else Path(download(PANEL_URL, raw / "integrated_call_samples_v3.20130502.ALL.panel")["path"])
    if not panel_file:
        if sha256_file(panel) != PANEL_SHA256: raise ValueError("1000 Genomes panel checksum mismatch")
        audit.downloads.append({"url": PANEL_URL, "path": str(panel), "bytes": panel.stat().st_size, "sha256": PANEL_SHA256, "reused": False}); audit.save()
    vcf = raw / "1kgp_chr1_sort1.vcf.gz"
    if vcf_file:
        _supplied(vcf_file, vcf, None, audit)
        source_index = Path(str(vcf_file) + ".tbi")
        if source_index.exists(): _supplied(source_index, Path(str(vcf) + ".tbi"), None, audit)
    elif not vcf.exists():
        subprocess.run(["bcftools", "view", "-r", "chr1:109209432-109340504", "-Oz", "-o", str(vcf), VCF_URL], check=True)
        subprocess.run(["bcftools", "index", "-t", str(vcf)], check=True)
        audit.downloads.append({"url": VCF_URL + "#chr1:109209432-109340504", "path": str(vcf), "bytes": vcf.stat().st_size, "sha256": sha256_file(vcf), "reused": False}); audit.save()
    return {"gtex": gtex, "vcf": vcf, "panel": panel}


def _orient(row: pd.Series, value: float, ref: str, alt: str) -> float:
    if row.liftover_strand == "-":
        tr = str.maketrans("ACGT", "TGCA"); ref, alt = ref.translate(tr), alt.translate(tr)
    if alt == row.ldl_lowering_allele: return value
    if ref == row.ldl_lowering_allele: return -value
    raise ValueError(f"GTEx alleles do not match GLGC alleles for {row.rsid}")


def attach_eqtl_and_tagging(run_dir: Path, inputs: dict[str, Path], base: pd.DataFrame) -> pd.DataFrame:
    """Merge GTEx liver eQTL effects and the rs12740374 LD-tagging covariate
    onto an AlphaGenome gene-mask score table (any model regime). Shared by
    Figure 1C (ALL_FOLDS) and Figure S2B (FOLD_0) -- everything below this
    point depends only on the variant table and the external GTEx/1000G
    inputs, not on which model produced `ag_rna_liver_{gene}`."""
    base = base.copy()
    wanted_pos = set(base.pos_hg19.astype(int)); by_gene: dict[str, dict[int, list[tuple[float,float,float,str,str]]]] = {g:{} for g in GENES}
    id_gene = {v:k for k,v in GENE_IDS.items()}
    with gzip.open(inputs["gtex"], "rt") as handle:
        header = handle.readline().rstrip().split("\t"); ix = {v:i for i,v in enumerate(header)}
        for line in handle:
            fields=line.rstrip().split("\t"); gene=id_gene.get(fields[ix["gene_id"]])
            if not gene: continue
            variant=fields[ix["variant_id"]].split("_"); pos=int(variant[1])
            if pos in wanted_pos: by_gene[gene].setdefault(pos,[]).append((float(fields[ix["slope"]]),float(fields[ix["slope_se"]]),float(fields[ix["pval_nominal"]]),variant[2],variant[3]))
    for gene in GENES:
        vals=[]
        for _, row in base.iterrows():
            candidates=by_gene[gene].get(int(row.pos_hg19),[]); x=None
            for candidate in candidates:
                ref,alt=candidate[3],candidate[4]
                if row.liftover_strand == "-":
                    tr=str.maketrans("ACGT","TGCA"); ref,alt=ref.translate(tr),alt.translate(tr)
                if {ref,alt} == {row.ldl_raising_allele,row.ldl_lowering_allele}: x=candidate; break
            if x is not None:
                vals.append((_orient(row,x[0],x[3],x[4]),x[1],x[2]))
            elif candidates:
                # Preserve uncertainty/P for a same-position GTEx record, but do
                # not assign a signed effect when its alleles differ from GLGC.
                vals.append((np.nan,candidates[0][1],candidates[0][2]))
            else:
                vals.append((np.nan,np.nan,np.nan))
        base[f"eqtl_liver_{gene}"],base[f"eqtl_liver_{gene}_se"],base[f"eqtl_liver_{gene}_p"]=zip(*vals)
    panel=pd.read_csv(inputs["panel"],sep="\t"); eur=set(panel.loc[panel.super_pop.eq("EUR"),"sample"].astype(str))
    query=["bcftools","query","-f","%POS\t%ID\t%REF\t%ALT[\t%SAMPLE=%GT]\n",str(inputs["vcf"])]
    records: dict[int, list[list[str]]] = {}
    for line in subprocess.run(query,check=True,text=True,capture_output=True).stdout.splitlines():
        f=line.split("\t"); records.setdefault(int(f[0]), []).append(f)
    tag=records[VARIANT_POS][0]; sample_gt=lambda f:{x.split("=",1)[0]:x.split("=",1)[1] for x in f[4:]}
    taggt=sample_gt(tag); tag_alt=str(tag[3]).split(",").index("T")+1
    cov=[]; ptag=[]; status=[]
    for _,row in base.iterrows():
        candidates=records.get(int(row.pos), [])
        f=next((x for x in candidates if {row.ldl_lowering_allele,row.ldl_raising_allele} == set([x[2],*x[3].split(",")])),None)
        if f is None: cov.append(np.nan); ptag.append(np.nan); status.append("missing_from_vcf"); continue
        alleles=[f[2],*f[3].split(",")]
        try: low=alleles.index(row.ldl_lowering_allele)
        except ValueError: cov.append(np.nan); ptag.append(np.nan); status.append("allele_mismatch"); continue
        vg=sample_gt(f); n=lo=ta=both=0
        for s in eur & vg.keys() & taggt.keys():
            a=vg[s].replace("/","|").split("|"); b=taggt[s].replace("/","|").split("|")
            if len(a)!=2 or len(b)!=2 or "." in a+b: continue
            for x,y in zip(a,b): n+=1; lo+=int(int(x)==low); ta+=int(int(y)==tag_alt); both+=int(int(x)==low and int(y)==tag_alt)
        ps,pt,pb=lo/n,ta/n,both/n
        cov.append((pb-ps*pt)/(ps*(1-ps)) if 0<ps<1 else np.nan); ptag.append(pb/ps if ps else np.nan); status.append("matched")
    base["tagging_ptag_EUR"],base["tagging_covvar_EUR"],base["tagging_match_status"]=ptag,cov,status
    causal=base.loc[base.rsid.eq(VARIANT_RSID)].iloc[0]
    for gene in GENES:
        direct=base[f"ag_rna_liver_{gene}"]; model=direct+base.tagging_covvar_EUR*float(causal[f"ag_rna_liver_{gene}"])
        model.loc[base.rsid.eq(VARIANT_RSID)]=float(causal[f"ag_rna_liver_{gene}"])
        base[f"ag_model_snp_plus_covvar_rs127_for_plot_{gene}"]=model.fillna(direct)
    return base


def build_full_figure1c(run_dir: Path, inputs: dict[str, Path]) -> None:
    base = pd.read_csv(run_dir / "derived/Figure1C_middle_ag_scores.tsv", sep="\t")
    result = attach_eqtl_and_tagging(run_dir, inputs, base)
    out=run_dir/"derived/Figure1C_eqtl_direct_tagging.tsv"; result.to_csv(out,sep="\t",index=False,float_format="%.10g"); render_figure1c(out,run_dir/"figures/Figure1C.svg")


def render_figure1c(source: Path, output: Path) -> None:
    d=pd.read_csv(source,sep="\t"); fig,axes=plt.subplots(3,3,figsize=(6.5,5),sharex=True)
    cols=[lambda g:f"eqtl_liver_{g}",lambda g:f"ag_rna_liver_{g}",lambda g:f"ag_model_snp_plus_covvar_rs127_for_plot_{g}"]
    for i,g in enumerate(GENES):
        for j,fn in enumerate(cols):
            c=fn(g); causal=d.rsid.eq(VARIANT_RSID); axes[i,j].scatter(d.loc[~causal,"pos_mb"],d.loc[~causal,c],s=7,c="#777"); axes[i,j].scatter(d.loc[causal,"pos_mb"],d.loc[causal,c],s=20,c="#d62728"); axes[i,j].axhline(0,c="black",lw=.4)
    for j,t in enumerate(["GTEx liver","AG direct","AG + EUR tagging"]): axes[0,j].set_title(t,fontsize=8)
    _save_svg(fig,output)


def run_observed_contact(run_dir: Path, audit: Audit, hic_file: Path | None) -> None:
    import hicstraw
    source=str(hic_file or HIC_URL); start,end,res=108750000,109800000,2000
    hic=hicstraw.HiCFile(source)
    zoom=hic.getMatrixZoomData("1","1","oe","KR","BP",res)
    matrix=np.asarray(zoom.getRecordsAsMatrix(start,end,start,end),dtype=float)
    pred=run_dir/"predictions/Figure1D_observed_hic"; pred.mkdir(parents=True,exist_ok=True); np.savez_compressed(pred/"observed.npz",matrix=matrix,start=start,resolution=res)
    centers=start+(np.arange(matrix.shape[0])+.5)*res; offsets=(centers-(VARIANT_POS-1))/1000; idx=np.flatnonzero((offsets>=-60)&(offsets<=165)); pos=offsets[idx]
    out=run_dir/"derived/Figure1D_observed_hic.tsv"; pd.DataFrame(matrix[np.ix_(idx,idx)],index=pos,columns=pos).rename_axis("row_position_from_rs12740374_kb").to_csv(out,sep="\t")
    _render_heatmap(out,run_dir/"figures/Figure1D.svg","Observed 4DN O/E")
    audit.downloads.append({"url":source,"path":"remote byte-range access" if not hic_file else str(hic_file),"bytes":hic_file.stat().st_size if hic_file else 14732267572,"sha256":sha256_file(hic_file) if hic_file else "S3 ETag ca98fe976d7321969696347d167ba35c-1757","reused":False}); audit.save()
    _build_figure1f(run_dir,matrix,start,res)


def _bin(pos,start,res): return (pos-1-start)//res
def _window(m,a,b): return float(np.nanmean(m[max(0,a-1):a+2,max(0,b-1):b+2]))
def _stats(scale,m,start,res):
    a=_bin(VARIANT_POS,start,res); rows=[]
    for gene,tss in TSS.items():
        b=_bin(tss,start,res); dist=abs(a-b); null=np.array([_window(m,i,i+dist) for i in range(1,m.shape[0]-dist-1)])
        val=_window(m,a,b); pct=100*(np.sum(null<val)+.5*np.sum(null==val))/len(null)
        rows.append(dict(scale=scale,gene=gene,resolution_bp=res,flank_bins=1,window_width_bp_each_axis=3*res,anchor_bin=a,target_bin=b,variant_to_tss_bp=tss-VARIANT_POS,contact_value=val,same_distance_n=len(null),same_distance_mean=np.mean(null),same_distance_sd=np.std(null,ddof=1),same_distance_z=(val-np.mean(null))/np.std(null,ddof=1),same_distance_percentile=pct,same_distance_high_p_plus1=(np.sum(null>=val)+1)/(len(null)+1)))
    return rows
def _build_figure1f(run_dir,obs,ostart,ores):
    p=np.load(run_dir/"predictions/Figure1E_fold0_contact/ref_alt.npz"); rows=_stats("Observed 4DN 2kb 3x3",obs,ostart,ores)+_stats("AlphaGenome FOLD_0 3x3",p["ref"],int(p["start"]),int(p["resolution"])); out=run_dir/"derived/Figure1F_promoter_contact_percentiles.tsv"; pd.DataFrame(rows).to_csv(out,sep="\t",index=False); d=pd.DataFrame(rows); fig,ax=plt.subplots(figsize=(4,2.5));
    for scale,g in d.groupby("scale"): ax.plot(g.gene,g.same_distance_percentile,"o-",label=scale)
    ax.set_ylabel("Same-distance percentile"); ax.legend(fontsize=6); _save_svg(fig,run_dir/"figures/Figure1F.svg")
def _render_heatmap(source,output,label):
    d=pd.read_csv(source,sep="\t"); p=d.iloc[:,0]; fig,ax=plt.subplots(figsize=(2.5,2.5)); im=ax.imshow(d.iloc[:,1:],origin="lower",cmap="inferno",extent=(p.iloc[0],p.iloc[-1],p.iloc[0],p.iloc[-1])); fig.colorbar(im,ax=ax,label=label); _save_svg(fig,output)
