configfile: "config/config.yaml"
SAMPLES = config["samples"]

print(expand("{sample}", sample=SAMPLES))

rule all:
    input:
        # expand("quality-check/fastqc/raw/{sample}-R1_fastqc.html", sample=SAMPLES),
        # expand("quality-check/fastqc/raw/{sample}-R1_fastqc.zip", sample=SAMPLES),
        # expand("quality-check/fastqc/raw/{sample}-R2_fastqc.html", sample=SAMPLES),
        # expand("quality-check/fastqc/raw/{sample}-R2_fastqc.zip", sample=SAMPLES),
        # "quality-check/multiqc/trim/multiqc_report.html",
        # expand("trimmed/{sample}-R1.fastq.gz", sample=SAMPLES),
        # expand("trimmed/{sample}-R2.fastq.gz", sample=SAMPLES),
        # expand("quality-check/fastqc/{mode}/{sample}-R1_fastqc.zip", sample=SAMPLES, mode="trimmed"),
        # expand("quality-check/fastqc/{mode}/{sample}-R2_fastqc.zip", sample=SAMPLES, mode="trimmed")
        # expand("alignment/{sample}.bam", sample=SAMPLES),
        # expand("alignment/{sample}.sam", sample=SAMPLES),
        # expand("alignment/{sample}.sorted.bam", sample=SAMPLES)
        "counts/counts.txt"


rule fastqc_raw:
    input:
        # R1 = lambda wildcards: config["samples"[wildcards.sample]["R1"]],
        R1 = "raw_data/{sample}-R1.fastq.gz",
        R2 = "raw_data/{sample}-R2.fastq.gz"
    output:
        R1_html = "quality-check/fastqc/raw/{sample}-R1_fastqc.html",
        R1_zip = "quality-check/fastqc/raw/{sample}-R1_fastqc.zip",
        R2_html = "quality-check/fastqc/raw/{sample}-R2_fastqc.html",
        R2_zip = "quality-check/fastqc/raw/{sample}-R2_fastqc.zip"
    benchmark:
        "benchmarks/{sample}-R1.fq.benchmark.txt"
    threads:
        4
    shell:
        """
        fastqc -t {threads} {input.R1} {input.R2} -o quality-check/fastqc/raw/
        """

rule multiqc_raw:
    input:
        expand("quality-check/fastqc/raw/{sample}-R1_fastqc.zip", sample=SAMPLES),
        expand("quality-check/fastqc/raw/{sample}-R2_fastqc.zip", sample=SAMPLES)
    output:
        directory("quality-check/multiqc/fq/multiqc_data"),
        "quality-check/multiqc/fq/multiqc_report.html"
    log:
        "logs/multiqc.log"
    shell:
        """
        (multiqc {input} -o quality-check/multiqc/fq/) 2> {log}
        """
        
rule fastp:
    input:
        "quality-check/multiqc/fq/multiqc_report.html",
        R1 = "raw_data/{sample}-R1.fastq.gz",
        R2 = "raw_data/{sample}-R2.fastq.gz"
    output:
        O1 = "trimmed/{sample}-R1.fastq.gz",
        O2 = "trimmed/{sample}-R2.fastq.gz",
        json = "reports/fastp/{sample}-fastp.json",
        html = "reports/fastp/{sample}-fastp.html"
    threads:
        4
    params:
        polyX = 10
    log:
        "logs/{sample}-fastp.log"
    shell:
        """
        (fastp -w {threads} -i {input.R1} -I {input.R2} -o {output.O1} -O {output.O2} -j {output.json} \
        -h {output.html} -x --poly_x_min_len {params.polyX}) 2> {log}
        """

rule fastqc_trimmed:
    input:
        # R1 = lambda wc: if wc.mode == 'trimmed': "trimmed/{wc.sample}-R1.fastq.gz" else "raw-data/{wc.sample}-R1.fastq.gz"
        # R2 = lambda wc: if wc.mode == 'trimmed': "trimmed/{wc.sample}-R2.fastq.gz" else "raw-data/{wc.sample}-R2.fastq.gz"
        R1 = "trimmed/{sample}-R1.fastq.gz",
        R2 = "trimmed/{sample}-R2.fastq.gz"
    output:
        R1_html = "quality-check/fastqc/trimmed/{sample}-R1_fastqc.html",
        R1_zip = "quality-check/fastqc/trimmed/{sample}-R1_fastqc.zip",
        R2_html = "quality-check/fastqc/trimmed/{sample}-R2_fastqc.html",
        R2_zip = "quality-check/fastqc/trimmed/{sample}-R2_fastqc.zip"
    threads:
        4
    shell:
        """
        fastqc -t {threads} {input.R1} {input.R2} -o quality-check/fastqc/trimmed/
        """

rule multiqc_trimmed:
    input:
        expand("quality-check/fastqc/trimmed/{sample}-R1_fastqc.zip", sample=SAMPLES),
        expand("quality-check/fastqc/trimmed/{sample}-R2_fastqc.zip", sample=SAMPLES)
    output:
        directory("quality-check/multiqc/trim/multiqc_data"),
        "quality-check/multiqc/trim/multiqc_report.html"
    log:
        "logs/multiqc.log"
    shell:
        """
        (multiqc {input} -o quality-check/multiqc/trim/) 2> {log}
        """

rule hisat2_index:
    input:
        "quality-check/multiqc/trim/multiqc_report.html",
        fa = config["ref_genome"]
    output:
        expand("{prefix}.{i}.ht2", prefix=config["genome_index"], i=range(1,9))
    log:
        "logs/hisat2_index.log"
    threads:
        4
    shell:
        """
        (hisat2-build -p {threads} {input.fa} {config[genome_index]}) 2> {log}
        """

rule hisat2_alignment:
    input:
        idx = expand("reference/index/genome.{i}.ht2", i=range(1,9)),
        R1 = "trimmed/{sample}-R1.fastq.gz",
        R2 = "trimmed/{sample}-R2.fastq.gz"
    output:
        sam = "alignment/{sample}.sam"
    log:
        "logs/hisat2_alignment_{sample}.log"
    params:
        idx_prefix = config["genome_index"]
    threads:
        4
    shell:
        """
        (hisat2 -x {params.idx_prefix} -1 {input.R1} -2 {input.R2} \
         -p {threads} -S {output.sam}) 2> {log}
        """

rule sam_to_bam:
    input:
        "alignment/{sample}.sam"
    output:
        "alignment/{sample}.bam"
    log:
        "logs/samtobam_{sample}.log"
    shell:
        """
        (samtools view -b {input} -o {output}) 2> {log}
        """

rule sort_index:
    input:
        "alignment/{sample}.bam"
    output:
        sorted = "alignment/{sample}.sorted.bam"
        # indexed = "alignment/{sample}.sorted.bam.bai"
    log:
        "logs/samtools_sort_{sample}.log"
    shell:
        """
        (samtools sort {input} -o {output.sorted} --write-index) 2> {log}
        """


## For populating read group information through picard:
# picard AddOrReplaceReadGroups \
#   -I 672.sorted.bam \
#   -O 672.sorted.rg.bam \
#   --RGID "H7M7WDSXF.1" \          # FlowcellID.Lane
#   --RGLB "TACTTAGC+TGGTACCT" \    # Dual index sequences (library identifier)
#   --RGPL "ILLUMINA" \             # Always "ILLUMINA" for NovaSeq/HiSeq
#   --RGPU "H7M7WDSXF.1.672" \      # FlowcellID.Lane.Sample (unique per sample)
#   --RGSM "672" \                  # Sample name (must match your sample ID)
#   --CREATE_INDEX true

rule feature_counts:
    input:
        expand("alignment/{sample}.sorted.bam", sample=SAMPLES)
    output:
        counts = "counts/counts.txt",
        summary = "counts/counts.txt.summary"
    threads:
        4
    params:
        strand = 2
    shell:
        """
        featureCounts -a {config[gtf]} -o {output.counts} \
        -p  -B -C -s {params.strand} -T {threads} {input}
        """
    
