BASE_DIR = "./"
REF_NAME = "glyma.Wm82.gnm1.FCtY.genome_main.fna"
REF = BASE_DIR + "ref/" + REF_NAME
THREAD=15
DEPTH=5

workdir: BASE_DIR

SAMPLES,  = glob_wildcards(BASE_DIR+"samples/{sample}_1.fastq")

rule all:
    input: "output/merge/merged.qc.vcf.gz"

rule alignment:
    input:
        fwd="samples/{sample}_1.fastq",
        rev="samples/{sample}_2.fastq",
    params: REF=REF, THREAD=THREAD
    output: "output/{sample}/{sample}.sam"
    shell: """
        mkdir -p output/{wildcards.sample}
        bwa mem -t {params.THREAD} {params.REF} {input.fwd} {input.rev} -o {output}

    """

rule fixmate:
    input: "output/{sample}/{sample}.sam"
    output: "output/{sample}/{sample}.bam"
    shell: "samtools fixmate -O bam {input} {output}"

rule sorting:
    input: "output/{sample}/{sample}.bam"
    params: THREAD=THREAD
    output: "output/{sample}/{sample}.sort.bam"
    shell: "sambamba sort -o {output} {input} --tmpdir output/{wildcards.sample}/tmp"

rule mpileup:
    input: "output/{sample}/{sample}.sort.bam"
    params: DEPTH=DEPTH, REF=REF
    output: "output/{sample}/{sample}.m.gvcf.gz"
    shell: "bcftools mpileup -g {params.DEPTH} -Oz -o {output} -f {params.REF} {input}"

rule call:
    input: "output/{sample}/{sample}.m.gvcf.gz"
    params: DEPTH=DEPTH
    output: "output/{sample}/{sample}.call.gvcf.gz"
    shell: "bcftools call -g {params.DEPTH} -Oz -o {output} {input}"

rule index:
    input: "output/{sample}/{sample}.call.gvcf.gz"
    output: "output/{sample}/{sample}.call.gvcf.gz"
    shell: "bcftools index {input}"

rule clean:
    shell: "rm output/{wildcards.sample}/{wildcards.sample}.sam output/{wildcards.sample}/{wildcards.sample}.bam output/{wildcards.sample}/{wildcards.sample}.m.gvcf.gz"

rule merge:
    input: expand("output/{sample}/{sample}.call.gvcf.gz", sample=SAMPLES)
    output: "output/merge/merged.gvcf.gz"
    params: REF=REF
    run:
        inputs = " ".join(input)
        shell("bcftools merge -g {params.REF} -m both -Oz -o {output} + inputs")

rule filtering:
    input: "output/merge/merged.gvcf.gz"
    params: FMISSING="F_MISSING < 0.05", MAF="MAF>0.02", QUAL="%QUAL>=30", MQ="MQ>=30", TYPE="TYPE='snp'"
    output: "output/merge/merged.qc.vcf.gz"
    shell: 'bcftools view -i "{params.FMISSING} && {params.MAF} && {params.QUAL} && {params.MQ} && {params.TYPE}" -Oz -o {output} {intput}'