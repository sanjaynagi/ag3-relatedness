rule set_kernel:
    input:
        f'{workflow.basedir}/envs/pythonGenomics.yaml'
    output:
        touch("resources/.kernel.set")
    conda: f'{workflow.basedir}/envs/pythonGenomics.yaml'
    log:
        "logs/set_kernel.log"
    shell: 
        """
        python -m ipykernel install --user --name probe 2> log
        """


def getVCFs(gz=True, allcontigs=False, allcontigsseparately=False):

    if allcontigs == False:
        if gz == True:
            genotypes = expand("resources/vcfs/{sample_set}.{{contig}}.vcf.gz", sample_set=sample_sets)
        elif gz == False:
            genotypes = expand("resources/vcfs/{sample_set}.{{contig}}.vcf", sample_set=sample_sets)
    elif allcontigs == True:
        if gz == True:
            genotypes =  expand("resources/vcfs/wholegenome/{sample_set}.vcf.gz", sample_set=sample_sets)
        elif gz == False:
            genotypes = expand("resources/vcfs/wholegenome/{sample_set}.vcf", sample_set=sample_sets)
    
    if allcontigsseparately:
        genotypes = expand(genotypes, contig=contigs)

    return(genotypes)



# def getSelectedOutputs(wildcards):

#     """
#     Function that returns a list of the desired outputs for the rule all, depending on the config.yaml
#     configuration file.
#     """

#     selected_input = []
   
#     if config['Relatedness']['activate']:
#         selected_input.extend(
#             expand(
#                 "results/relatedness/ngsRelate.{sample_set}",
#             sample_set=sample_sets,
#             )
#         )

#     return(selected_input)



def singleTrue(iterable):
    iterator = iter(iterable)
#    # consume from "i" until first true or it's exhausted
    has_true = any(iterator) 
#    # carry on consuming until another true value / exhausted
    has_another_true = any(iterator) 
#    # True if exactly one true found
    return has_true and not has_another_true