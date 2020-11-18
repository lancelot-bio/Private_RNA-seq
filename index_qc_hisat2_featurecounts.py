REP_INDEX = {"sample1","sample2","sample3","sample4","sample5","sample6","sample7","sample8"}
FA = {"/home/junzhang/dataset/genome/7.Ananas_comosus"}
FQ = {"/home/junzhang/0.cleandata"}
HI2_INDEX = {"/home/junzhang/sample/0.index/Acomosus_321_v3.fa.HISAT2.index"}
ANNO_PATH = {"/home/junzhang/dataset/genome/7.Ananas_comosus/Acomosus_321_v3.gene_exons.gtf"}

rule all:
    input:
        expand("QC/{rep}_1.fastp",rep=REP_INDEX),
        expand("QC/{rep}_2.fastp",rep=REP_INDEX),
        expand("QC/{rep}.fastp.html",rep=REP_INDEX),
        expand("QC/{rep}.fastp.json",rep=REP_INDEX),
        expand("0.index/"),
        expand("1.align/{rep}.hisat2.sam.summary",rep=REP_INDEX),
        expand("1.align/BAM/{rep}.sorted.bam",rep=REP_INDEX),
        expand("2.counts/counts.featurecoundts")



rule index:
    input:
       fa = expand("{fa_path}/Acomosus_321_v3.fa",fa_path=FA)
    output:
        "0.index/Acomosus_321_v3.fa.HISAT2.index"
    log:
        "0.index/Acomosus_321_v3.fa.HISAT2.index.log"
    shell:
        "hisat2-build {input[0]} {output[0]} -p 24 > {log} 2>&1 "
#HISAT2-build

rule qc:
    input:
        "/home/junzhang/0.cleandata/{rep}_R1.fq",
        "/home/junzhang/0.cleandata/{rep}_R2.fq"
    output:
        "QC/{rep}_1.fastp",
        "QC/{rep}_2.fastp",
        "QC/{rep}.fastp.html",
        "QC/{rep}.fastp.json"
    log:
        "QC/{rep}.fastp.log"
    shell:
        "fastp -w 6 -i {input[10]}} -I {input[1]]} -o {output[0]} -O {output[1]} -h  {output[2]} -j {output[3]} > {log} 2>&1"
# fastp, 6 thread, output html and json

rule HISAT2:
    input:
        "QC/{rep}_1.fastp",
        "QC/{rep}_2.fastp"
    output:
        "1.align/{rep}.hisat2.sam",
        "1.align/{rep}.hisat2.sam.summary"
    log:
        "1.align/{rep}.hisat2.log"
    shell:
        "hisat2 -x {HI2_INDEX} -1 {input[0]} -2 {input[1]} -S {output[0]} --new-summary --summary-file {output[1]} -p 24 > {log} 2>&1"

rule Samtools:
    input:
        "1.align/{rep}.hisat2.sam"
    output:
        "1.align/BAM/{rep}.sorted.bam"
    log:
        "1.align/BAM/{rep}.bam.log"
    shell:
        "samtools sort -O BAM -o {output[0]} -@ 8 {input[0]} > {log} 2>&1"
#Samtools V1.7

rule featurecounts:
    input:
        bam = expand("1.align/BAM/{rep}.sorted.bam",rep=REP_INDEX)
    output:
        "2.counts/counts.featurecoundts"
    log:
        "2.counts/counts.featurecoundts.log"
    shell:
        "featureCounts -t exon -a {ANNO_PATH} -T 6 -o {output} {input.bam} > {log} 2>&1 "

