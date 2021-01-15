################################################################################################
# PROTAX for aligned sequences (both train and test data must follow the same global alignment)
# About the statistical model associated with each taxon node:
# for simplicity, here using best and 2nd best similarity as predictors
# (for other choices, do the required changes in c functions)
################################################################################################

###########################
# Preliminaries
###########################

# 1) compile c functions and copy binaries to scripts directory

cd c
make
cp create_xdata_best2 trainclassify_best2 classify_best2 ../scripts
cd ..

# 2) set environmental variable PROTAX to point to scripts directory
#    e.g. replace path information ../ based on your computer system

PROTAX=`pwd`/scripts

###########################
# Training the model
###########################

# prepare following 3 files:
#
# - taxonomy file (tree structure for all taxa in all levels in text format): nid pid level tname
# - seqid2tax file (taxonomy information of each sequence in text format): seqid maxlevel tname
# - refs.aln file (aligned sequences in FASTA format, MUST use capital letters A,C,G,T and - for gap)
#                 (one sequence must be in one line, not in multiple lines, i.e. total number of lines in file is 2*number_of_seqs)

# 0) create directory for model parameters and go there, also copy taxonomy, seqid2tax, and refs.aln there

mkdir mymodel
cd mymodel
cp ../taxonomy .
cp ../seqid2tax .
cp ../refs.aln .

# 1) set variable NUM_TAXLEVELS based on your taxonomy, e.g. if there are 4 levels then:

NUM_TAXLEVELS=4

# 2) add priors for unknown taxa in taxonomy, one value for each level
#    this relates to how much you think there are taxa not included in your taxonomy,
#    larger values add more uncertainty to all predictions
#    unk prior is level-specific where units correspond to leaf nodes of the taxonomy (value * prior(known_species))
#    NOTE: the number of priors in the ,,, list must equal to $NUM_TAXLEVELS

perl $PROTAX/taxonomy_priors.pl 10,1,0.1,0.01 taxonomy > taxonomy.priors

for ((LEVEL=1; LEVEL<=$NUM_TAXLEVELS; LEVEL++))
do
 perl $PROTAX/thintaxonomy.pl $LEVEL taxonomy.priors > tax$LEVEL
done

# 3) generate training data for each level, here 10000 training samples per level

for ((LEVEL=1; LEVEL<=$NUM_TAXLEVELS; LEVEL++))
do
 echo "LEVEL $LEVEL"
 perl $PROTAX/seqid2taxlevel.pl $LEVEL seqid2tax > ref.tax$LEVEL
 perl $PROTAX/get_all_reference_sequences.pl $LEVEL tax$LEVEL ref.tax$LEVEL rseqs$LEVEL
 perl $PROTAX/taxrseq2numeric.pl $LEVEL tax$LEVEL refs.aln > rseqs${LEVEL}.numeric
 perl $PROTAX/generate_training_data.pl $LEVEL tax$LEVEL ref.tax$LEVEL rseqs$LEVEL 10000 1 no train$LEVEL 
 perl $PROTAX/traindat2numeric.pl refs.aln train$LEVEL > train${LEVEL}.numeric
done

# to check what kind of training data there is for each level:

for ((LEVEL=1; LEVEL<=$NUM_TAXLEVELS; LEVEL++))
do
 echo "LEVEL $LEVEL:"
 cut -f4 -d" " train$LEVEL | cut -f1 -d"," | sort | uniq -c
done

# 4) calculate xdat file (sequence similarity predictors), scale the values and save the scaling parameters for later use
#    ...this can take a while if large training data...

for ((LEVEL=1; LEVEL<=$NUM_TAXLEVELS; LEVEL++))
do
 echo "LEVEL $LEVEL"
 $PROTAX/create_xdata_best2 tax$LEVEL refs.aln rseqs${LEVEL}.numeric train${LEVEL}.numeric > train${LEVEL}.xdat
done

for ((LEVEL=1; LEVEL<=$NUM_TAXLEVELS; LEVEL++))
do
 perl $PROTAX/scale_xdat.pl sc$LEVEL train${LEVEL}.xdat > train${LEVEL}.scxdat
done

# 5) parameter estimation in R
#    MCMC for parameters separately in each taxonomy level
#    you need to check the convergence and continue iterations or re-initialize adaptive proposal if needed

R
# use correct path to R scripts based on your system
source("../scripts/amcmc.rcode_noweight.txt")
library(compiler)
logprior=cmpfun(logprior)
loglikelihood=cmpfun(loglikelihood)
adaptiveMCMC=cmpfun(adaptiveMCMC)

num.params=1+4
ind=1001:2000

dev.new(width=10,height=5)

dat=read.xdata("train1.scxdat")
pp1=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1)
traceplot.all(pp1,ind,num.levels=1, title="L1")

# if problems with mixing, usually one or two re-adaptations are enough, run commented lines
# initstate=initialize.adaptation(pp1$params[2000,])
# pp1=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1,prev.state=initstate)
# traceplot.all(pp1,ind,num.levels=1, title="L1 readapted")

k=which.max(pp1$postli[ind])
write.postparams(pp1,"mcmc1",ind[k])

dat=read.xdata("train2.scxdat")
pp1=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1)
traceplot.all(pp1,ind,num.levels=1, title="L2")

# if problems with mixing, usually one or two re-adaptations are enough, run commented lines
# initstate=initialize.adaptation(pp1$params[2000,])
# pp1=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1,prev.state=initstate)
# traceplot.all(pp1,ind,num.levels=1, title="L2 readapted")

k=which.max(pp1$postli[ind])
write.postparams(pp1,"mcmc2",ind[k])

dat=read.xdata("train3.scxdat")
pp1=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1)
traceplot.all(pp1,ind,num.levels=1, title="L3")

# if problems with mixing, usually one or two re-adaptations are enough, run commented lines
# initstate=initialize.adaptation(pp1$params[2000,])
# pp1=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1,prev.state=initstate)
# traceplot.all(pp1,ind,num.levels=1, title="L3 readapted")

k=which.max(pp1$postli[ind])
write.postparams(pp1,"mcmc3",ind[k])

dat=read.xdata("train4.scxdat")
pp1=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1)
traceplot.all(pp1,ind,num.levels=1, title="L4")

# if problems with mixing, usually one or two re-adaptations are enough, run commented lines
# initstate=initialize.adaptation(pp1$params[2000,])
# pp1=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1,prev.state=initstate)
# traceplot.all(pp1,ind,num.levels=1, title="L4 readapted")

k=which.max(pp1$postli[ind])
write.postparams(pp1,"mcmc4",ind[k])

# continue to all levels you have L5, L6, ..., L$NUM_TAXLEVELS

##################################
# Check model with training data
##################################

#### check classification with training samples (in this example using those which were used for training level 4)

perl $PROTAX/init_train_taxprob2numeric.pl train4 > query0.prob.numeric

cut -f3-6 -d" " mcmc1 > par1
cut -f3-6 -d" " mcmc2 > par2
cut -f3-6 -d" " mcmc3 > par3
cut -f3-6 -d" " mcmc4 > par4

# parent probs from previous level classification
for LEVEL in 1 2 3 4
do
 echo "LEVEL $LEVEL" 
 PREVLEVEL=$((LEVEL-1))
 IFILE=query${PREVLEVEL}.prob.numeric
 OFILE=query${LEVEL}.prob
 # add scaling
 $PROTAX/trainclassify_best2 tax$LEVEL refs.aln rseqs${LEVEL}.numeric par$LEVEL sc$LEVEL 0.01 train4.numeric $IFILE > $OFILE
 perl $PROTAX/train_taxprob2numeric.pl tax$LEVEL $OFILE > ${OFILE}.numeric
done

##### calculate correctness (note: correct label is the one which training sample mimicked, it can be e.g. unknown species)

for LEVEL in 1 2 3 4
do
 perl $PROTAX/trainsample2correct.pl $LEVEL taxonomy train4 > query${LEVEL}.tax
 perl $PROTAX/trainsample2addcor.pl query${LEVEL}.prob query${LEVEL}.tax > query${LEVEL}.cor
done

##### bias accuracy plots
R
source("../scripts/amcmc.rcode_noweight.txt")

nimi=c("Phylum","Class","Order","Family")

par(mfrow=c(2,2))
for (i in 1:4) {
 file=sprintf("query%d.cor",i)
 a=read.table(file,header=F)
 accuracy.plot(a[,3],a[,4],name=sprintf("Model performance: %s",nimi[i]))
}

#####################################
# PROTAX classification for test data
#####################################

# put level-specific parameters in single file in correct order

echo -n "" > model.pars
echo -n "" > model.scs
echo -n "" > model.rseqs.numeric

for ((LEVEL=1; LEVEL<=$NUM_TAXLEVELS; LEVEL++))
do
 cut -f3-6 -d" " mcmc$LEVEL >> model.pars
 cat sc$LEVEL >> model.scs
 cat rseqs${LEVEL}.numeric >> model.rseqs.numeric
done

# test sequences in FASTA file 'test.fa'
# NOTE: test.fa must follow the same alignment as refs.aln
# for each test sequence, output list of taxon,probability pairs
# NOTE: output includes all taxonomy levels NOT sorted by probability

$PROTAX/classify_best2 taxonomy.priors refs.aln model.rseqs.numeric model.pars model.scs 0.01 test.fa > test.out

