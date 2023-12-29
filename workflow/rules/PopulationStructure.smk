
rule ngsRelate:
    """
    Run NGSRelate on VCF files
    """
    input:
        vcf = "resources/vcfs/wholegenome/{sample_set}.vcf.gz",
        csi = "resources/vcfs/wholegenome/{sample_set}.vcf.gz.csi",
        tbi = "resources/vcfs/wholegenome/{sample_set}.vcf.gz.tbi",
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