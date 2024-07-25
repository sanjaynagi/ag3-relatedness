rule malariagen_dataToVCF:
    """
    Write out haplotypes VCF files from provided malariagen_data
    """
    output:
        haplotypeVCF = config['release'] + "_results/vcfs/{sample_set}_{contig}.vcf"
    conda:
        "../envs/pythonGenomics.yaml"
    log:
        "logs/ZarrToVCF/{sample_set}_{contig}.log"
    resources:
         tot=1
    priority: 10
    params:
        basedir=workflow.basedir,
        dataset=dataset,
        release=config['release'],
        haplotypes=config['haplotypes'],
    script:
        "../scripts/ZarrToVCF.py"

rule BGZip:
    """
    This is overwriting log files at the
    """
    input:
        calls = config['release'] + "_results/vcfs/{sample_set}_{contig}.vcf"
    output:
        calls_gz = config['release'] + "_results/vcfs/{sample_set}_{contig}.vcf.gz"
    log:
        "logs/bgzip/{sample_set}/{contig}.log"
    shell:
        """
        bgzip {input.calls} 2> {log}
        """

rule BcftoolsIndex:
    input:
        calls = config['release'] + "_results/vcfs/{sample_set}_{contig}.vcf.gz"
    output:
        calls_gz = config['release'] + "_results/vcfs/{sample_set}_{contig}.vcf.gz.csi",
    log:
        "logs/bcftoolsIndex/{sample_set}_{contig}.log",
    shell:
        """
        bcftools index {input.calls} 2> {log}
        """

rule Tabix:
    input:
        calls = config['release'] + "_results/vcfs/{sample_set}_{contig}.vcf.gz"
    output:
        calls_tbi = config['release'] + "_results/vcfs/{sample_set}_{contig}.vcf.gz.tbi",
    log:
        "logs/tabix/{sample_set}_{contig}.log",
    shell:
        """
        tabix {input.calls} 2> {log}
        """


rule concatVCFs:
    input:
        calls = expand(config['release'] + "_results/vcfs/{{sample_set}}_{contig}.vcf.gz", contig=config["contigs"]),
        tbi = expand(config['release'] + "_results/vcfs/{{sample_set}}_{contig}.vcf.gz.tbi", contig=config["contigs"]),
        csi = expand(config['release'] + "_results/vcfs/{{sample_set}}_{contig}.vcf.gz.csi", contig=config["contigs"]),
    output:
        cattedVCF = config['release'] + "_results/vcfs/wholegenome/{sample_set}.vcf.gz",
    log:
        "logs/bcftoolsConcat/{sample_set}.log",
    threads: 8
    shell:
        """
        bcftools concat -o {output.cattedVCF} -O z --threads {threads} {input.calls} 2> {log}
        """


rule BcftoolsIndex_cattedVCF:
    input:
        calls = config['release'] + "_results/vcfs/wholegenome/{sample_set}.vcf.gz"
    output:
        calls_gz = config['release'] + "_results/vcfs/wholegenome/{sample_set}.vcf.gz.csi",
    log:
        "logs/bcftoolsIndex/{sample_set}.log",
    shell:
        """
        bcftools index {input.calls} 2> {log}
        """


rule Tabix_cattedVCF:
    input:
        calls = config['release'] + "_results/vcfs/wholegenome/{sample_set}.vcf.gz"
    output:
        calls_tbi = config['release'] + "_results/vcfs/wholegenome/{sample_set}.vcf.gz.tbi",
    log:
        "logs/tabix/{sample_set}.log",
    shell:
        """
        tabix {input.calls} 2> {log}
        """