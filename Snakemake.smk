import snakemake.parser
import json
from os.path import join, basename, dirname

# globals ----------------------------

configfile: 'config.yml'
# Full path to an uncompressed FASTA file with all chromosome sequences.
GTF = config['GTF']
GENOME = config['GENOME']
# Full path to a folder where intermediate output files will be created.
OUT_DIR = config['OUT_DIR']

FILES = json.load(open(config['SAMPLES_JSON']))

SAMPLES = sorted(FILES.keys())


rule all:
        input:
                join(dirname(GENOME), 'STAR', basename(GENOME).rstrip(".fasta")),
                expand(OUT_DIR + "/STAR/"+"{sample}.Aligned.sortedByCoord.out.bam",sample = SAMPLES),
                OUT_DIR + "/gene_count_matrix.csv",


rule STAR_index:
        input:
                gtf = GTF,
                genome = GENOME
        output:
            index = directory(join(dirname(GENOME), 'STAR', basename(GENOME).rstrip(".fasta")))
        log:
            'logs/STAR_index.log'
        shell:
                """
                STAR --runThreadN 20 \
                --runMode genomeGenerate  \
                --genomeDir {output.index} \
                --genomeFastaFiles {input.genome}  \
                --sjdbGTFfile {input.gtf}
                """
rule STAR:
        input:
                r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
                r2 = lambda wildcards: FILES[wildcards.sample]['R2'],
                index = rules.STAR_index.output.index
        output:
                bam = OUT_DIR + "/STAR/{sample}.Aligned.sortedByCoord.out.bam",
                trans_bam = OUT_DIR + "/STAR/{sample}.Aligned.toTranscriptome.out.bam",
                log_file = OUT_DIR + "/STAR/{sample}.Log.final.out"
        params:
                outFileNamePrefix = OUT_DIR + "/STAR/{sample}."
        shell:
                """
                STAR --runThreadN 20 \
                    --genomeDir {input.index} \
                    --outSAMtype BAM SortedByCoordinate \
                    --outBAMcompression 9 \
                    --limitBAMsortRAM 60000000000 \
                    --readFilesCommand zcat \
                    --quantMode TranscriptomeSAM GeneCounts \
                    --readFilesIn {input.r1} {input.r2} \
                    --outFileNamePrefix {params.outFileNamePrefix}
                """

rule STRINGTIE:
        input:
                bam = rules.STAR.output.bam,
                gtf = GTF
        output:
                out = OUT_DIR +"/ballgown/{sample}/{sample}.gtf"
        shell: "stringtie -e -B -p 8 -G {input.gtf} -o {output} {input.bam}"

rule count_tables:
            input:
                expand(OUT_DIR +"/ballgown/{sample}/{sample}.gtf", sample=SAMPLES)
            output:
                output_g = OUT_DIR + "/gene_count_matrix.csv",
            params:
                StringtieScript = config["PrepDE"],
                output_result = config["OUT_DIR"] +"/ballgown"
            shell:
                """
                python {params.StringtieScript} -i {params.output_result} -g {output.output_g} -l 150
                """