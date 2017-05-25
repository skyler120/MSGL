
#!/bin/bash
# Telling how many nodes and processors should be used.
#PBS -l nodes=1:ppn=8
# Naming the file
#PBS -N run-sims-sgl
# Outputting error
#PBS -j oe
# Not sure what the two next lines do
#PBS -q default
#PBS -S /bin/bash
#PBS -m abe
#PBS -M ss3349@cornell.edu

######## I WONDER IF I NEED THIS LINE ######
cd $PBS_O_WORKDIR

# Telling cluster that you are using R
R --vanilla > run-sims-sgl.out <<EOF

# Looking for what machines are available to use.
setwd("/home/fs01/ss3349/MSGL")
#install all packages needed



#################################################################################
#  All of the above 'R --vanilla...' is for the cluster
#  All of the below 'R --vanilla...' is an R file
#  This is the beginning of a 'regular' R file
#################################################################################

######## HERE BEGINS THE SIMULATION ########




EOF
