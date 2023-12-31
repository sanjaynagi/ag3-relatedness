rule ZarrToHaplotypesVCF:
    """
    Write out haplotypes VCF files from provided malariagen_data
    """
    output:
        haplotypeVCF = "results/vcfs/{sample_set}_{contig}.vcf"
    conda:
        "../envs/pythonGenomics.yaml"
    log:
        "logs/ZarrToVCF_haplotypes/{sample_set}_{contig}.log"
    resources:
        tot=1
    priority: 10
    params:
        basedir=workflow.basedir,
        dataset=dataset
    script:
        "../scripts/ZarrToVCF_haplotypes.py"


gzippedVCF = getVCFs(gz=True)
rule BGZip:
    """
    This is overwriting log files at the
    """
    input:
        calls = getVCFs(gz=False)
    output:
        calls_gz = gzippedVCF
    log:
        "logs/bgzip/{contig}.log"
    priority: 100
    shell:
        """
        bgzip {input.calls} 2> {log}
        """

rule BcftoolsIndex:
    input:
        calls = getVCFs(gz=True)
    output:
        calls_gz = "results/vcfs/{sample_set}_{contig}.vcf.gz.csi",
    log:
        "logs/bcftoolsIndex/{sample_set}_{contig}.log",
    shell:
        """
        bcftools index {input.calls} 2> {log}
        """

rule Tabix:
    input:
        calls = getVCFs(gz=True)
    output:
        calls_tbi = "results/vcfs/{sample_set}_{contig}.vcf.gz.tbi",
    log:
        "logs/tabix/{sample_set}_{contig}.log",
    shell:
        """
        tabix {input.calls} 2> {log}
        """


rule concatVCFs:
    input:
        calls = expand("results/vcfs/{{sample_set}}_{contig}.vcf.gz", contig=config["contigs"]),
        tbi = expand("results/vcfs/{{sample_set}}_{contig}.vcf.gz.tbi", contig=config["contigs"]),
        csi = expand("results/vcfs/{{sample_set}}_{contig}.vcf.gz.csi", contig=config["contigs"]),
    output:
        cattedVCF = "results/vcfs/wholegenome/{sample_set}.vcf.gz",
    log:
        "logs/bcftoolsConcat/{sample_set}.log",
    threads: 8
    shell:
        """
        bcftools concat -o {output.cattedVCF} -O z --threads {threads} {input.calls} 2> {log}
        """


rule BcftoolsIndex_cattedVCF:
    input:
        calls = lambda wildcards: getVCFs(gz=True, allcontigs=True)
    output:
        calls_gz = "results/vcfs/wholegenome/{sample_set}.vcf.gz.csi",
    log:
        "logs/bcftoolsIndex/{sample_set}.log",
    shell:
        """
        bcftools index {input.calls} 2> {log}
        """


rule Tabix_cattedVCF:
    input:
        calls = lambda wildcards: getVCFs(gz=True, allcontigs=True)
    output:
        calls_tbi = "results/vcfs/wholegenome/{sample_set}.vcf.gz.tbi",
    log:
        "logs/tabix/{sample_set}.log",
    shell:
        """
        tabix {input.calls} 2> {log}
        """