

rule GenomeIndex:
    """
    Index the reference genome with samtools
    """
    input:
        ref=config["reference"]["genome"],
    output:
        idx=config["reference"]["genome"] + ".fai",
    log:
        "logs/GenomeIndex.log",
    wrapper:
        "v0.69.0/bio/samtools/faidx"

rule ZarrToHaplotypesVCF:
    """
    Write out haplotypes VCF files from provided malariagen_data
    """
    output:
        haplotypeVCF = expand("resources/vcfs/{dataset}_{{contig}}.haplotypes.vcf", dataset=dataset)
    conda:
        "../envs/pythonGenomics.yaml"
    log:
        "logs/ZarrToVCF_haplotypes/{contig}.log"
    params:
        ag3_sample_sets = ag3_sample_sets,
        basedir=workflow.basedir,
        dataset=dataset
    script:
        "../scripts/ZarrToVCF_haplotypes.py"


gzippedVCF = getVCFs(gz=True, bothAllelisms=True)
rule BGZip:
    """
    This is overwriting log files at the
    """
    input:
        calls = getVCFs(gz=False, bothAllelisms=True)
    output:
        calls_gz = gzippedVCF
    log:
        "logs/bgzip/{contig}_{allelism}.log" if config['VCF']['activate'] is False else "logs/bgzip/{contig}.log"
    shell:
        """
        bgzip {input.calls} 2> {log}
        """

rule BcftoolsIndex:
    input:
        calls = getVCFs(gz=True)
    output:
        calls_gz = "resources/vcfs/{dataset}_{contig}.{allelism}.vcf.gz.csi",
    log:
        "logs/bcftoolsIndex/{dataset}_{contig}.{allelism}.log",
    shell:
        """
        bcftools index {input.calls} 2> {log}
        """

rule Tabix:
    input:
        calls = getVCFs(gz=True)
    output:
        calls_tbi = "resources/vcfs/{dataset}_{contig}.{allelism}.vcf.gz.tbi",
    log:
        "logs/tabix/{dataset}_{contig}_{allelism}.log",
    shell:
        """
        tabix {input.calls} 2> {log}
        """


rule concatVCFs:
    input:
        calls = lambda wildcards: getVCFs(gz=True, allelism=wildcards.allelism, allcontigs=False, allcontigsseparately=True),
        tbi = lambda wildcards: [vcf+".tbi" for vcf in getVCFs(gz=True, allelism=wildcards.allelism, allcontigs=False, allcontigsseparately=True)],
        csi = lambda wildcards: [vcf+".csi" for vcf in getVCFs(gz=True, allelism=wildcards.allelism, allcontigs=False, allcontigsseparately=True)],
    output:
        cattedVCF = "resources/vcfs/wholegenome/{dataset}.{allelism}.vcf.gz",
    log:
        "logs/bcftoolsConcat/{dataset}.{allelism}.log",
    threads: 8
    shell:
        """
        bcftools concat -o {output.cattedVCF} -O z --threads {threads} {input.calls} 2> {log}
        """


rule BcftoolsIndex_cattedVCF:
    input:
        calls = lambda wildcards: getVCFs(gz=True, allelism=wildcards.allelism, allcontigs=True)
    output:
        calls_gz = "resources/vcfs/wholegenome/{dataset}.{allelism}.vcf.gz.csi",
    log:
        "logs/bcftoolsIndex/{dataset}.{allelism}.log",
    shell:
        """
        bcftools index {input.calls} 2> {log}
        """

rule Tabix_cattedVCF:
    input:
        calls = lambda wildcards: getVCFs(gz=True, allelism=wildcards.allelism, allcontigs=True)
    output:
        calls_tbi = "resources/vcfs/wholegenome/{dataset}.{allelism}.vcf.gz.tbi",
    log:
        "logs/tabix/{dataset}_{allelism}.log",
    shell:
        """
        tabix {input.calls} 2> {log}
        """