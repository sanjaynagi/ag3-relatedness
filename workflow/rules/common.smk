rule set_kernel:
    input:
        f'{workflow.basedir}/envs/pythonGenomics.yaml'
    output:
        touch("results/.kernel.set")
    conda: f'{workflow.basedir}/envs/pythonGenomics.yaml'
    log:
        "logs/set_kernel.log"
    shell: 
        """
        python -m ipykernel install --user --name probe 2> log
        """




def singleTrue(iterable):
    iterator = iter(iterable)
#    # consume from "i" until first true or it's exhausted
    has_true = any(iterator) 
#    # carry on consuming until another true value / exhausted
    has_another_true = any(iterator) 
#    # True if exactly one true found
    return has_true and not has_another_true