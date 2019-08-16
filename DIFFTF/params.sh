###############################
# TEMPLATE FOR PARAMS.SH FILE #
#      ALWAYS UP TO DATE      #
###############################

# GENERAL NOTE: All paths can be relative as well as absolute. Usually, relative paths should be prefered

#################
# INPUT OPTIONS #
#################

# Config file to use
configFile="config_hg19.json"

# Snakefile to use
snakefile="/g/scb2/zaugg/carnold/Projects/diffTF/src/Snakefile"

#######################
# PERFORMANCE OPTIONS #
#######################

# Use a dry run for testing purposes?
dryRun=false

# Number of cores to use. For local execution, the minimum of this the real number of cores and this number is taken.
# When in cluster mode, the maximum number of CPUs per rule. See the separate cluster specification
nCores=16

# Disable locking of the .snakemake directory. Make to to NEVER run more than 1 Snakemake analysis for the same folder in parallel. Usually, there is a check for that, this switch disables it for performance reasons
nolock=true

# Shadow directory. See Snakemake help
shadowDir=""

# Enable various additional output of Snakemake such as verbose messages and printing shell commands
useVerbose=false

#################################
# CONDA AND SINGULARITY OPTIONS #
#################################
useConda=false

# Directory in which all conda environments are stored. As we use them repeatedly, a designated directory to store them in is useful. Defaults to "/g/scb2/zaugg/zaugg_shared/Programs/Snakemake/conda" unless specified otherwise
condaDir="/g/scb2/zaugg/zaugg_shared/Programs/Snakemake/conda"

# Should singularity be used?
useSingularity=true

# Singularity prefix. Only relevant if useSingularity=true. Defaults to "/g/scb2/zaugg/zaugg_shared/Programs/Snakemake/singularity"
singularityPrefix="/g/scb2/zaugg/zaugg_shared/Programs/Snakemake/singularity"

# Additional arguments for singularity. Only relevant if useSingularity=true
singularityArgs="--bind /g/scb,/g/scb2,/scratch"


###################
# CLUSTER OPTIONS #
###################

submitToCluster=true

# Use SLURM mode? This can always be set to TRUE currently and is redundant, as only SLURM is supported.
useSLURM=true

# Cluster configuration file
clusterConfig="cluster.json"

# Maximum number of simultaenous jobs
maxJobsCluster=300

# Maximum number of times a job should be reexecuted when failing. For SLURM, 0 is usually fine for SLURM, as it is reliable enough
maxRestartsPerJob=0

################
# RULE OPTIONS #
################

# Default: "", which means use all rules defined in the Snakefile.
# If you want to use only a particular set of rules, specify them here by their rule name and separate them by spaces
allowedRules=""

# Start from only a specific rule? If yes, name the rule here, which has to correspond to the rule name in the Snakefile
runSpecificRule=""

# Run also all downstream rules when runSpecificRule is specified? If set to false, ONLY the specified rule will be run; otherwise, all downstream rules are also triggered
runAlsoDownstreamRules=true

# Should rules marked as incomplete be rerun automatically? Usually, this is a good idea.
rerunIncomplete=true

# Should the WHOLE pipeline be rerun? use with care, setting this to true will rerun everything
forceRerunAll=false

# Should Snakemake abort after the first error or continue running?
abortAfterFirstError=false


#################
# FILES OPTIONS #
#################

# Use --notemp for developing purposes: Ignore temp() declarations. This is useful when running only a part of the workflow,
# since temp() would lead to deletion of probably needed files by other parts of the workflow.
ignoreTemp=false

# Only touches output files. useful when you manually copied data and you want to ensure that rules are not run unnecessarily because of timestamp reasons.
touchOutputFiles=false


#######################
# DEVELOPMENT OPTIONS #
#######################

#Ask Christian
ignoreZeroSizedFiles=true

runCustomCommand=""

# Only execute the last part: Calling Snakemake
skipSummaryAndDAG=true

# You can usually leave this untouched unless you also suffer from the "bug" that dot has a problem producing PDF files. Then change to "svg" pr "png"
workflowGraphFileType="pdf"
