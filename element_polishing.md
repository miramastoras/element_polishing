# Polishing with element data

This document has experiments for polishing with element deepvariant calls

## 1. Workflow outline

1. Align element data to each haplotype, remove reads with too much divergence
2. Call variants with deepvariant on element bams and hifi bam
3. Filter element variants:
    - GQ
    - same reads aligning to same place on both haplotypes
4. Phase element variants
5.


## 2. Just polishing with homozygous element calls to correct "false hets":

Plan here: https://docs.google.com/presentation/d/1diAuHi9iJHjgWwmgzR0yjtPeLw7yESthG5F5cr9aB1w/edit#slide=id.g25e8d33c0ac_0_72



### 2.1 Filter alignments by de tag

subset to one chromosome
```
#!/bin/bash
#SBATCH --job-name=subbam
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=main
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=64
#SBATCH --output=%x.%j.log
#SBATCH --time=4:00:00

cd /private/groups/patenlab/mira/hprc_polishing/data/element_HG002/hprc_y2/alignments/minimap2
samtools view -@ 64 -b -h HG002_element_50x_all2pat.mm2.srt.bam h1tg000010l > HG002_element_50x_all2pat.mm2.srt.h1tg000010l.bam

samtools index HG002_element_50x_all2pat.mm2.srt.h1tg000010l.bam
```

```
time python3 /private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/de_dist/plot_de_distribution.py -b /private/groups/patenlab/mira/hprc_polishing/data/element_HG002/hprc_y2/alignments/minimap2/HG002_element_50x_all2pat.mm2.srt.h1tg000010l.bam -p /private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/de_dist/HG002_element_50x_all2pat.h1tg000010l.mm2
```

### 2.2 Polishing with different filters

#### 2.2.1 All homozygous "PASS" variants

Subset vcf to all 1/1 genotype calls
```
cd /private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/filter_experiments/all_PASS_hom_calls

# Maternal
grep "^#" HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.element_50X.deepvariant_1.5.WGS.PASS.vcf > HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf
grep -v "^#" HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.element_50X.deepvariant_1.5.WGS.PASS.vcf | grep "1/1" >> HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf

bcftools stats HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf | head -n 30
# SN	0	number of records:	5273

# Paternal
grep "^#" HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.deepvariant_1.5.WGS.PASS.vcf > HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf
grep -v "^#" HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.deepvariant_1.5.WGS.PASS.vcf | grep "1/1" >> HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf
# SN	0	number of records:	3927
```
Apply polish + dipcall + bedtools intersect
```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/workflows/applyPolish_dipcall.wdl
```
```
{
  "applyPolish_dipcall.hap2Fasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa",
  "applyPolish_dipcall.confidenceBedFile": "/private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed",
  "applyPolish_dipcall.hap1PolishingVcf": "/private/groups/patenlab/mira/hprc_polishing/data/element_HG002/hprc_y2/deepvariant/bwa-mem_bams/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf",
  "applyPolish_dipcall.hap2PolishingVcf": "/private/groups/patenlab/mira/hprc_polishing/data/element_HG002/hprc_y2/deepvariant/bwa-mem_bams/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf",
  "applyPolish_dipcall.sampleID": "HG002.trio_hifiasm_0.19.5.DC_1.2_40x.element.PASS.homalt.",
  "applyPolish_dipcall.referenceFasta": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
  "applyPolish_dipcall.dipcall_t.referenceIsHS38": true,
  "applyPolish_dipcall.referenceFastaFai": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai",
  "applyPolish_dipcall.dipcall_t.isMaleSample": true,
  "applyPolish_dipcall.hap1Fasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa"
}
```
```
cd /private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/filter_experiments/all_PASS_hom_calls

mkdir -p logs && time SINGULARITY_CACHEDIR=`pwd`/outputs/cache/.singularity/cache MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/outputs/cache/.cache/miniwdl toil-wdl-runner --logDebug --jobStore ./big_store --batchSystem slurm --batchLogsDir ./logs /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/workflows/applyPolish_dipcall.wdl applyPolish_dipcall_inputs.json -o applyPolish_dipcall_outputs -m applyPolish_dipcall_outputs.json  2>&1 | tee applyPolish_dipcall_log.txt
```

Happy through slurm
```
#!/bin/bash
#SBATCH --job-name=happy_elementhomalt
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=main
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=2:00:00

docker run -it --rm -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz /private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/filter_experiments/all_PASS_hom_calls/applyPolish_dipcall_outputs/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.element.PASS.homalt._hap1.polished.dipcall.vcf.gz -r /private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta -f /private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/filter_experiments/all_PASS_hom_calls/applyPolish_dipcall_outputs/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent_intersect_HG002.trio_hifiasm_0.19.5.DC_1.2_40x.element.PASS.homalt._hap1.polished.dipcall.bed -o /private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/filter_experiments/all_PASS_hom_calls/happy_out --pass-only --no-roc --no-json --engine=vcfeval --threads=32
```

Intersection with DeepPolisher edits
```
bcftools merge --force-samples HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf.gz HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf.gz -o HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.mat.combined.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf.gz

bcftools isec /private/groups/patenlab/mira/hprc_polishing/data/element_HG002/hprc_y2/deepvariant/bwa-mem_bams/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.mat.combined.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf.gz /private/groups/patenlab/mira/hprc_polishing/hprc_deepPolisher_wf_runs/HG002_y2_DCv1.2_PHv5_DPmodel5/HG002_y2_DCv1.2_PHv5_DPmodel5_polisher_output.vcf.gz -p ./
```

#### 2.2.1 All homozygous "PASS" variants, different GQ filters

```
bcftools view -f "PASS" -e 'FORMAT/GQ<=5' /private/groups/patenlab/mira/hprc_polishing/data/element_HG002/hprc_y2/deepvariant/bwa-mem_bams/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf.gz -Oz -o /private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/filter_experiments/homalt_calls_PASS_GQ_5/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.deepvariant_1.5.WGS.PASS.homalt.GQ5.vcf.gz

bcftools view -f "PASS" -e 'FORMAT/GQ<=5' /private/groups/patenlab/mira/hprc_polishing/data/element_HG002/hprc_y2/deepvariant/bwa-mem_bams/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf.gz -Oz -o /private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/filter_experiments/homalt_calls_PASS_GQ_5/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.element_50X.deepvariant_1.5.WGS.PASS.homalt.GQ5.vcf.gz


# 5915 calls / 9100 calls

bcftools view -f "PASS" -e 'FORMAT/GQ<=10' /private/groups/patenlab/mira/hprc_polishing/data/element_HG002/hprc_y2/deepvariant/bwa-mem_bams/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf.gz -Oz -o /private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/filter_experiments/homalt_calls_PASS_GQ_10/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.deepvariant_1.5.WGS.PASS.homalt.GQ10.vcf.gz

bcftools view -f "PASS" -e 'FORMAT/GQ<=10' /private/groups/patenlab/mira/hprc_polishing/data/element_HG002/hprc_y2/deepvariant/bwa-mem_bams/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf.gz -Oz -o /private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/filter_experiments/homalt_calls_PASS_GQ_10/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.element_50X.deepvariant_1.5.WGS.PASS.homalt.GQ10.vcf.gz

# 4263 / 9100 calls

bcftools view -f "PASS" -e 'FORMAT/GQ<=20' /private/groups/patenlab/mira/hprc_polishing/data/element_HG002/hprc_y2/deepvariant/bwa-mem_bams/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf.gz -Oz -o /private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/filter_experiments/homalt_calls_PASS_GQ_20/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.deepvariant_1.5.WGS.PASS.homalt.GQ20.vcf.gz

bcftools view -f "PASS" -e 'FORMAT/GQ<=20' /private/groups/patenlab/mira/hprc_polishing/data/element_HG002/hprc_y2/deepvariant/bwa-mem_bams/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf.gz -Oz -o /private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/filter_experiments/homalt_calls_PASS_GQ_20/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.element_50X.deepvariant_1.5.WGS.PASS.homalt.GQ20.vcf.gz


# 825 / 9100 calls

bcftools view -f "PASS" -e 'FORMAT/GQ<=30' /private/groups/patenlab/mira/hprc_polishing/data/element_HG002/hprc_y2/deepvariant/bwa-mem_bams/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.mat.combined.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf.gz -Oz -o /private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/filter_experiments/homalt_calls_PASS_GQ_30/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.mat.combined.element_50X.deepvariant_1.5.WGS.PASS.homalt.GQ30.vcf.gz
# 10
```
```
{
  "applyPolish_dipcall.hap2Fasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa",
  "applyPolish_dipcall.confidenceBedFile": "/private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed",
  "applyPolish_dipcall.hap1PolishingVcf": "/private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/filter_experiments/homalt_calls_PASS_GQ_20/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.deepvariant_1.5.WGS.PASS.homalt.GQ20.vcf.gz",
  "applyPolish_dipcall.hap2PolishingVcf": "/private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/filter_experiments/homalt_calls_PASS_GQ_20/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.element_50X.deepvariant_1.5.WGS.PASS.homalt.GQ20.vcf.gz",
  "applyPolish_dipcall.sampleID": "HG002.trio_hifiasm_0.19.5.DC_1.2_40x.element.PASS.homalt.GQ20",
  "applyPolish_dipcall.referenceFasta": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
  "applyPolish_dipcall.dipcall_t.referenceIsHS38": true,
  "applyPolish_dipcall.referenceFastaFai": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai",
  "applyPolish_dipcall.dipcall_t.isMaleSample": true,
  "applyPolish_dipcall.hap1Fasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa"
}
```
```
mkdir -p logs && time SINGULARITY_CACHEDIR=`pwd`/outputs/cache/.singularity/cache MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/outputs/cache/.cache/miniwdl toil-wdl-runner --logDebug --jobStore ./big_store --batchSystem slurm --batchLogsDir ./logs /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/workflows/applyPolish_dipcall.wdl applyPolish_dipcall_inputs.json -o applyPolish_dipcall_outputs -m applyPolish_dipcall_outputs.json  2>&1 | tee applyPolish_dipcall_log.txt
```

Running hap.py on all three GQ filters
```
bash /private/home/mmastora/progs/scripts/HG002_happy.sh /private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/filter_experiments/homalt_calls_PASS_GQ_5/applyPolish_dipcall_outputs/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.element.PASS.homalt.GQ5_hap1.polished.dipcall.vcf.gz /private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/filter_experiments/homalt_calls_PASS_GQ_5/applyPolish_dipcall_outputs/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent_intersect_HG002.trio_hifiasm_0.19.5.DC_1.2_40x.element.PASS.homalt.GQ5_hap1.polished.dipcall.bed /private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/filter_experiments/homalt_calls_PASS_GQ_5/homalt_calls_PASS_GQ_5_happy_out

bash /private/home/mmastora/progs/scripts/HG002_happy.sh /private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/filter_experiments/homalt_calls_PASS_GQ_10/applyPolish_dipcall_outputs/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.element.PASS.homalt.GQ10_hap1.polished.dipcall.vcf.gz /private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/filter_experiments/homalt_calls_PASS_GQ_10/applyPolish_dipcall_outputs/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent_intersect_HG002.trio_hifiasm_0.19.5.DC_1.2_40x.element.PASS.homalt.GQ10_hap1.polished.dipcall.bed /private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/filter_experiments/homalt_calls_PASS_GQ_10/homalt_calls_PASS_GQ_10_happy_out

bash /private/home/mmastora/progs/scripts/HG002_happy.sh /private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/filter_experiments/homalt_calls_PASS_GQ_20/applyPolish_dipcall_outputs/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.element.PASS.homalt.GQ20_hap1.polished.dipcall.vcf.gz /private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/filter_experiments/homalt_calls_PASS_GQ_20/applyPolish_dipcall_outputs/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent_intersect_HG002.trio_hifiasm_0.19.5.DC_1.2_40x.element.PASS.homalt.GQ20_hap1.polished.dipcall.bed /private/groups/patenlab/mira/hprc_polishing/element_polishing/homozygous_alleles/filter_experiments/homalt_calls_PASS_GQ_20/homalt_calls_PASS_GQ_20_happy_out
```
### 2.3 Implementing method to filter out variants whose reads aren't aligned to the same haplotype

#### Step 1: Project variants from one haplotype to the other

Align haplotypes together in paf format
```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/tasks/Asm2AsmAlignerPaf.wdl
```

```
{
  "asm2asmAlignerPaf.alignmentPaf.aligner": "minimap2",
  "asm2asmAlignerPaf.alignmentPaf.refAssembly": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa",
  "asm2asmAlignerPaf.alignmentPaf.options": "--cs -c --eqx",
  "asm2asmAlignerPaf.alignmentPaf.readFastq_or_queryAssembly": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa",
  "asm2asmAlignerPaf.alignmentPaf.suffix": "mat2pat"
}
```
```
mkdir -p logs && time SINGULARITY_CACHEDIR=`pwd`/outputs/cache/.singularity/cache MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/outputs/cache/.cache/miniwdl toil-wdl-runner --jobStore ./big_store2 --batchSystem slurm --batchLogsDir ./logs /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/tasks/Asm2AsmAlignerPaf.wdl asm2asmAlignerPaf_inputs.json -o toil_asm2asm_paf_out -m asm2asmAlignerPaf_outputs.json  2>&1 | tee toil_asm2asm_paf_log.txt
```

Convert both element vcf files to bed files without adding 10 bp to the indels
```
/private/groups/patenlab/mira/hprc_polishing/data/element_HG002/hprc_y2/deepvariant/bwa-mem_bams/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.element_50X.deepvariant_1.5.WGS.vcf.gz

/private/groups/patenlab/mira/hprc_polishing/data/element_HG002/hprc_y2/deepvariant/bwa-mem_bams/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.deepvariant_1.5.WGS.vcf.gz

export PATH=$PATH:/private/home/mmastora/progs/bin/
/private/home/mmastora/progs/bin/vcf2bed --do-not-split < /private/groups/patenlab/mira/hprc_polishing/data/element_HG002/hprc_y2/deepvariant/bwa-mem_bams/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.element_50X.deepvariant_1.5.WGS.PASS.vcf > /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.element_50X.deepvariant_1.5.WGS.vcf.bed

/private/home/mmastora/progs/bin/vcf2bed --do-not-split < /private/groups/patenlab/mira/hprc_polishing/data/element_HG002/hprc_y2/deepvariant/bwa-mem_bams/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.deepvariant_1.5.WGS.PASS.vcf > /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.deepvariant_1.5.WGS.vcf.bed

```
Project both to each haplotype
```
# project pat to mat

docker run --rm -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira mobinasri/flagger:latest python3 /home/programs/src/project_blocks_multi_thread.py --threads 10 --mode 'ref2asm' --paf /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/toil_asm2asm_paf_out/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat2pat.paf --blocks /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.deepvariant_1.5.WGS.vcf.bed --outputProjectable /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/pat_variants_projectable_to_mat.bed --outputProjection /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/pat_variants_projection_to_mat.bed

# project mat to pat
docker run --rm -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira mobinasri/flagger:latest python3 /home/programs/src/project_blocks_multi_thread.py --threads 10 --mode 'asm2ref' --paf /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/toil_asm2asm_paf_out/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat2pat.paf --blocks /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.element_50X.deepvariant_1.5.WGS.vcf.bed --outputProjectable /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/mat_variants_projectable_to_pat.bed --outputProjection /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/mat_variants_projection_to_pat.bed
```
Check how the projection looks in examples from https://docs.google.com/presentation/d/1fWFw_SnReOX2Um_wrX1iU37Bo-RJJNdf21YcHOyK0J0/edit#slide=id.g28b3561553a_0_85

Run projection on just homalt PASS variants to see how many project over

#### Step 2

input: projection
```
/private/home/mmastora/progs/bin/vcf2bed --do-not-split < /private/groups/patenlab/mira/hprc_polishing/data/element_HG002/hprc_y2/deepvariant/bwa-mem_bams/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.mat.combined.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf > /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.mat.combined.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf.bed

awk '{print $1"\t"$2-10"\t"$3+10"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.mat.combined.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf.bed > /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.mat.combined.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf.10bp.bed

# asking mobin why some variants are unprojectable
# the bed window size to project will need to be a parameter to test for this script

grep "h1" /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.mat.combined.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf.10bp.bed > test.bed

docker run --rm -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira mobinasri/flagger:latest python3 /home/programs/src/project_blocks_multi_thread.py --threads 10 --mode 'ref2asm' --paf /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/toil_asm2asm_paf_out/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat2pat.paf --blocks /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.mat.combined.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf.bed --outputProjectable /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/pat_variants_projectable_to_mat.bed --outputProjection /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/pat_variants_projection_to_mat.bed

# project mat to pat
docker run --rm -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira mobinasri/flagger:latest python3 /home/programs/src/project_blocks_multi_thread.py --threads 10 --mode 'asm2ref' --paf /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/toil_asm2asm_paf_out/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat2pat.paf --blocks /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.mat.combined.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf.bed --outputProjectable /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/mat_variants_projectable_to_pat.bed --outputProjection /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/mat_variants_projection_to_pat.bed
```

#### Step 3: Come up with a test set for writing the code

Choose one contig, one bam file
```
#!/bin/bash
#SBATCH --job-name=subbam
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=main
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=16
#SBATCH --output=%x.%j.log
#SBATCH --time=3:00:00

cd /private/groups/patenlab/mira/hprc_polishing/data/element_HG002/hprc_y2/alignments/bwa-mem

samtools view -@ 16 -b -h HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.srt.bam h1tg000001l > HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.h1tg000001l.srt.bam

samtools index HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.h1tg000001l.srt.bam
```
- starting with no bp expansion on either side
- paste projection and projectables side by side, cat mat and pats which is input to script,
```
paste -d "\t" /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/pat_variants_projectable_to_mat.bed /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/pat_variants_projection_to_mat.bed > /private/groups/patenlab/mira/hprc_polishing/element_polishing/check_reads_hap_aln/pat_variants_projected_to_mat.bed

paste -d "\t" /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/pat_variants_projectable_to_mat.bed /private/groups/patenlab/mira/hprc_polishing/element_polishing/link_haplotypes/mat_variants_projection_to_pat.bed > /private/groups/patenlab/mira/hprc_polishing/element_polishing/check_reads_hap_aln/mat_variants_projected_to_pat.bed

cat /private/groups/patenlab/mira/hprc_polishing/element_polishing/check_reads_hap_aln/pat_variants_projected_to_mat.bed /private/groups/patenlab/mira/hprc_polishing/element_polishing/check_reads_hap_aln/mat_variants_projected_to_pat.bed > variants_projected.bed

grep "h1tg000001l" variants_projected.bed
```

- give script same bamfile for both haplotypes, "projections" that have the same variant are the same and ones with different variants will be shifted upstream or completely different locations


Test bamfile: `/private/groups/patenlab/mira/hprc_polishing/data/element_HG002/hprc_y2/alignments/bwa-mem/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.h1tg000001l.srt.bam`
Test set:
```
h1tg000001l	95040761	95040762	.	5.6	T	C	PASS	.	GT:GQ:DP:AD:VAF:PL	./.:.:.:.:.:.	1/1:3:7:0,7:1:2,4,0	h1tg000001l	95040761	95040762	.	5.6	T	C	PASS	.	GT:GQ:DP:AD:VAF:PL	./.:.:.:.:.:.	1/1:3:7:0,7:1:2,4,0
h1tg000001l	4641590	4641591	.	22.4	C	CTT	PASS	.	GT:GQ:DP:AD:VAF:PL	./.:.:.:.:.:.	1/1:5:28:5,23:0.821429:20,2,0	h1tg000001l	4641590	4641591	.	22.4	C	CTT	PASS	.	GT:GQ:DP:AD:VAF:PL	./.:.:.:.:.:.	1/1:5:28:5,23:0.821429:20,2,0
h1tg000001l	76358599	76358600	.	4.1	C	T	PASS	.	GT:GQ:DP:AD:VAF:PL	./.:.:.:.:.:.	1/1:4:6:0,6:1:1,11,0	h1tg000001l	76358599	76358600	.	4.1	C	T	PASS	.	GT:GQ:DP:AD:VAF:PL	./.:.:.:.:.:.	1/1:4:6:0,6:1:1,11,0
h1tg000001l	7127121	7127122	.	19.4	A	AAT	PASS	.	GT:GQ:DP:AD:VAF:PL	./.:.:.:.:.:.	1/1:15:34:2,31:0.911765:19,17,0	h1tg000001l	156227861	156227862	.	19.4	A	AAT	PASS	.	GT:GQ:DP:AD:VAF:PL	./.:.:.:.:.:.	1/1:15:34:2,31:0.911765:19,17,0
h1tg000001l	8382273	8382274	.	4.4	G	C	PASS	.	GT:GQ:DP:AD:VAF:PL	./.:.:.:.:.:.	1/1:2:5:1,4:0.8:0,3,0	h1tg000001l	8392273	8392274	.	4.4	G	C	PASS	.	GT:GQ:DP:AD:VAF:PL	./.:.:.:.:.:.	1/1:2:5:1,4:0.8:0,3,0
h1tg000001l	112169074	112169075	.	19.4	A	AT	PASS	.	GT:GQ:DP:AD:VAF:PL	./.:.:.:.:.:.	1/1:17:125:2,123:0.984:19,21,0	h1tg000001l	112169100	112169101	.	19.4	A	AT	PASS	.	GT:GQ:DP:AD:VAF:PL	./.:.:.:.:.:.	1/1:17:125:2,123:0.984:19,21,0
```

- keeping original vcf entry by pasting it back in after vcf2bed step
- paste projection and projectables side by side,cat mat and pats which is input to script,
- for every bed block, extract read names from bamfile on each haplotype
- compare two lists, only print variant to new vcf file if they are different

```
python3 check_haplotype_read_alignments.py --hap1Bam /Users/miramastoras/Desktop/element_polishing_files/check_reads_aligned/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.h1tg000001l.srt.bam --hap2Bam /Users/miramastoras/Desktop/element_polishing_files/check_reads_aligned/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.h1tg000001l.srt.bam --hap1Blocks /Users/miramastoras/Desktop/element_polishing_files/check_reads_aligned/test_projections.bed --hap2Blocks /Users/miramastoras/Desktop/element_polishing_files/check_reads_aligned/test_projections.bed --inVcf /Users/miramastoras/Desktop/element_polishing_files/check_reads_aligned/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.deepvariant_1.5.WGS.PASS.homalt.vcf --outVcf /Users/miramastoras/Desktop/element_polishing_files/check_reads_aligned/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.element_50X.deepvariant_1.5.WGS.PASS.homalt.filt.vcf
```

Confirmed it worked on the test set. Now running it on the whole bam file and set of projections to see how many homozygous variants are kept.

Using the exact projections, without expanding by 10bp:


## 3. Integrating heterozygous element calls: correcting false homozygous regions

**Linking haplotypes**

- align two haplotypes together, make paf file
- for all variant calls on both haplotypes, find their location on the other haplotype with project blocks
    look in IGV at a few examples of these to understand the projection and how many bp to expand by
- parse output of both projections:
- build some type of data structure to hold the different categories of alleles and flag what the polishing edit would do
- one haplotype has no variant call, one haplotype has a change back to that allele

Write a function to loop through each variant call, find its projection on the other haplotype, and catalog the alleles before and after polishing with element. Find if there is a variant call on the other haplotype in that exact location, and remove it from the list.
