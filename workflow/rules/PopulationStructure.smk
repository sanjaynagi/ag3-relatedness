rule mask_inversions:
    """
    Mask inversions in vcf files
    """
    input:
        vcf = "results/vcfs/wholegenome/{sample_set}.vcf.gz",
        csi = "results/vcfs/wholegenome/{sample_set}.vcf.gz.csi",
        tbi = "results/vcfs/wholegenome/{sample_set}.vcf.gz.tbi",
    output:
        vcf = "results/vcfs/wholegenome.masked/{sample_set}.vcf.gz"
    log:
        log = "logs/mask_inversions/{sample_set}.log"
    params:
        mask_regions = "^" + "^,".join(config['mask_region'])
    shell:
        """
        bcftools view -t {params.mask_regions} -O z -o {output.vcf} {input.vcf} 2> {log}
        """

rule BcftoolsIndex_masked:
    input:
        vcf = "results/vcfs/wholegenome.masked/{sample_set}.vcf.gz"
    output:
        vcf_gz = "results/vcfs/wholegenome.masked/{sample_set}.vcf.gz.csi"
    log:
        "logs/bcftoolsIndex/{sample_set}_masked.log",
    shell:
        """
        bcftools index {input.vcf} 2> {log}
        """

rule Tabix_masked:
    input:
        vcf = "results/vcfs/wholegenome.masked/{sample_set}.vcf.gz"
    output:
        vcf_tbi = "results/vcfs/wholegenome.masked/{sample_set}.vcf.gz.tbi",
    log:
        "logs/tabix/{sample_set}.masked.log"
    shell:
        """
        tabix {input.vcf} 2> {log}
        """

rule ngsRelate:
    """
    Mask inversions in vcf files
    """
    input:
        vcf = "results/vcfs/wholegenome.masked/{sample_set}.vcf.gz",
        csi = "results/vcfs/wholegenome.masked/{sample_set}.vcf.gz.csi",
        tbi = "results/vcfs/wholegenome.masked/{sample_set}.vcf.gz.tbi",
    output:
        "results/relatedness/ngsRelate.{sample_set}"
    log:
        log = "logs/ngsRelate/{sample_set}.log"
    params:
        tag = 'GT',
        basedir=workflow.basedir,
    threads: 24
    shell:
        """
        {params.basedir}/scripts/NgsRelate/ngsRelate -h {input.vcf} -O {output} -c 1 -T {params.tag} -p {threads} 2> {log}
        """