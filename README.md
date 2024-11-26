
# ----- WELCOME TO THE COMET PIPELINE ----- #


The coMet pipeline (standing for COMmunity METrics) is a pipeline developped to create communities following assembly rules and compute phylogenetic diversity metrics based on these communties. 
It is written mainly in R (the community simulator is written in C) and operated by command line. It is, by nature, designed to be launched on remote machines in slurm environement to create numerous replicates.
However, it is possible to operate it completely on local machines. 

# ---------------------- #
TO USE IT: 

The first step is the modification on the defaut configuration file "Foo.R" following the user needs. This configuration file will be used to pilot the integrity of the pipeline. 
Parameters details are described directly in this file. 

The second step is to launch the pipeline using this command line: ./Bash_coMet_Pipeline.sh Foo.R "Foo" 
The first argument is the name of the configuration file, the second argument is the name of the Scenario.

# ---------------------- #
NOTES: 
- Script 1 and 6 are designed to be used on a -slurm- remote machine to allow the creation of multiple replicates. However, you can still launch these two script for one replicate at a time on your local machine.
- This pipeline will not be further upgraded. 
