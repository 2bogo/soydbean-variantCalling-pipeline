BASE_DIR = "/data/snakemake_tutorial/"
REF_NAME = "glyma.Wm82.gnm1.FCtY.genome_main.fna"
REF = BASE_DIR + "ref/" + REF_NAME

workdir: BASE_DIR

SAMPLES,  = glob_wildcards(BASE_DIR+"samples/raw/{smaple}_1.fastq")

rule all:
    input: "samples/gather/merged.qc.gvcf.gz"

rule alignment:
    input: fwd="samples/raw/{sample}_1.fastq", rev="samples/raw/{sample}_2.fastq"
    params: REF=REF, THREAD=15
    output: "samples/{sample}/{sample}.sam"
    shell: """
        mkdir -p samples/{wildcards.sample}
        mkdir -p samples/gather
        bwa mem -t {params.THREAD} {params.REF} {input.fwd} {input.rev} -o {output}
    """

rule fixmate:
    input: "samples/{sample}/{sample}.sam"
    output: "samples/{sample}/{sample}.bam"
    shell: "samtools fixmate -O bam {input} {output}"

rule sorting:
    input: "samples/{sample}/{sample}.bam"
    output: "samples/{sample}/{sample}.sort.bam"
    shell: "sambamba sort -o  {input} {output} --tmpdir samples/{wildcards.sample}/tmp"

rule mpileup:
    input: "samples/{sample}/{sample}.sort.bam"
    params: DEPTH=5, REF=REF
    output: "samples/{sample}/{sample}.gvcf.gz"
    shell: "bcftools mpileup -g {params.DEPTH} -Oz -o {output} -f {params.REF} {input}"

rule call:
    input: "samples/{sample}/{sample}.gvcf.gz"
    params: DEPTH=5
    output: "samples/{sample}/{sample}.call.gvcf.gz"
    shell: "bcftools call -g {params.DEPTH} -Oz -o {output} {input}"

rule index:
    input: "samples/{sample}/{sample}.call.gvcf.gz"
    output: "samples/{sample}/{sample}.call.gvcf.gz"
    shell: """
        rm samples/{wildcards.sample}/{wildcards.sample}.sam samples/{wildcards.sample}/{wildcards.sample}.bam samples/{wildcards.sample}/{wildcards.sample}.gvcf.gz
        bcftools index {input}
        """

rule merge:
    input: expand("samples/{sample}/{sample}.call.gvcf.gz", sample=SAMPLES)
    output: "samples/gather/merged.gvcf.gz"
    params: REF=REF
    run:
        inputs = " ".join(input)
        shell("bcftools merge -g {params.REF} -m both -Oz -o {output}" + inputs)

rule filtering:
    input: "samples/gather/merged.gvcf.gz"
    params: FMISSING="F_MISSING < 0.05", MAF="MAF>0.02", QUAL="%QUAL<30", MQ="MQ<30"
    output: "samples/gather/merged.qc.gvcf.gz"
    shell: 'bcftools view -i "{params.FMISSING} && {params.MAF}" -e "N_ALT!=0 && {params.QUAL} && {params.MQ}" -Oz -o {output} {intput}'