# Going to try and make the tips align for this tree
setwd("/Users/izzynovick/Desktop/Tineoid_Moth_Project/UCE Analysis/Ancestral state recon")
install.packages("ape", dependencies = TRUE)
library("ape")
library(phytools)


# Going to try and reload a new tree formatted in figtree
# Read in the new ladderized left tree
test_tree<-readNexus(file="test_ladderized_tree_for_asr",
                     format="raxml")
plotTree(test_tree)

# Read in new data, edited to fit labels
moth_data2<-read.csv("synanthropy_char_matrix_csv_2.csv")

# Trying to reorer data so in the same order as the tree:
# 1. Get the species names from the phylogenetic tree
tree_species2 <- test_tree$tip.label
tree_species2

print(test_tree$tip.label)

# Print species names from the CSV data
print(moth_data2$Species)

# Match the species in the csv to the order of the species in the tree
matched_species <- intersect(test_tree$tip.label, moth_data2$Species)
matched_species

# Reorder rows of CSV file based on order of species names in tree
data_reordered <- moth_data2[match(matched_species, moth_data2$Species), ]
data_reordered
class(data_reordered)

# Write reordered data back to a CSV file
write.csv(data_reordered, "reordered_data.csv", row.names = FALSE)

# Making new dataframes for keratinophagy and synanthropy
keratinophagy_df2 <- data_reordered[, c("Species","Keratinophagy")]
synanthropy_df2 <- data_reordered[, c("Species","Synanthropy")]
keratinophagy_df2
synanthropy_df2

# Set the rownames to be the species names
rownames(keratinophagy_df2) <- keratinophagy_df2$Species
rownames(synanthropy_df2) <- synanthropy_df2$Species

# Make the new data frames to be used
# Habitat
habitat_data2<-setNames(synanthropy_df2$Synanthropy,
                        rownames(synanthropy_df2))
class(habitat_data2)
habitat_data2

# Diet
diet_data2<-setNames(keratinophagy_df2$Keratinophagy,
                     rownames(keratinophagy_df2))
class(diet_data2)

diet_data2

############## Now going to rerun all the previous models with the new tree.

# Fitting an equal rates model to the keratinophagy data and fitzjohn prior
diet_er2<-fitMk(test_tree,diet_data2,model="ER",
               pi="fitzjohn")
diet_er2

# Fitted (or set) value of Q:
# 0         1         2
# 0 -0.041393  0.020696  0.020696
# 1  0.020696 -0.041393  0.020696
# 2  0.020696  0.020696 -0.041393

# Fitted (or set) value of pi:
#  0        1        2 
# 0.944438 0.041194 0.014368 

# Log-likelihood: -24.526861 

# Fitting the ARD model with fitzjohn prior for diet

diet_ard2<-fitMk(test_tree,diet_data2,model="ARD",
                pi="fitzjohn")
diet_ard2

# Fitted (or set) value of Q:
# 0         1          2
# 0 -0.042909  0.032085  0.010824
# 1  0.000000 -0.032853  0.032853
# 2  0.054735  0.000000 -0.054735

# Fitted (or set) value of pi:
# 0 1 2 
# 0.900084 0.001028 0.098889 

# Log-likelihood: -23.302427 

# This one has a little better log-likelihood

# Going to fit a flat prior and compare it to Fitzjohn 
flatprior.fit_diet_ard2<-fitMk(test_tree,diet_data2,model="ARD")
flatprior.fit_diet_ard2

# Fitted (or set) value of Q:
# 0         1          2
# 0 -0.058942  0.040116  0.018826
# 1  0.000000 -0.042340  0.042340
# 2  0.198795  0.000000 -0.198795

# Fitted (or set) value of pi:
# 0        1        2 
# 0.333333 0.333333 0.333333 

# Log-likelihood: -23.555095


# Fitting an equal rates model to the synanthropy data
habitat_er2<-fitMk(test_tree,habitat_data2,model="ER",
                  pi="fitzjohn")
habitat_er2

# Fitted (or set) value of Q:
#  0          1          2
# 0 -0.161259  0.080630  0.080630
# 1  0.080630 -0.161259  0.080630
# 2  0.080630  0.080630 -0.161259

# Fitted (or set) value of pi:
# 0        1        2 
# 0.366189 0.378035 0.255776 

# Log-likelihood: -35.664653 

# Fitting the ARD model with habitat and fitzjohn prior
habitat_ard2<-fitMk(test_tree,habitat_data2,model="ARD",
                   pi="fitzjohn")
habitat_ard2

# Fitted (or set) value of Q:
#  0          1          2
# 0 -0.190385  0.190385  0.000000
# 1  0.187255 -0.197337  0.010082
# 2  0.002426  0.000000 -0.002426

# Log-likelihood: -28.835819 

# ARD had a better log-likelihood than ER for habitat/synanthropy

# Let’s compare the diet models using the normal machinery of likelihood – in this case, AIC.

anova(diet_er2,diet_ard2)

# log(L) d.f.      AIC    weight
# diet_er2  -24.52686    1 51.05372 0.97758962
# diet_ard2 -23.30243    6 58.60485 0.02241038

# For diet, most of the weight (97%) falls within the ER model, and AIC is lower (better) for ER model

# Let’s compare the habitat models using the normal machinery of likelihood – in this case, AIC.
anova(habitat_er2,habitat_ard2)

# log(L) d.f.      AIC    weight
# habitat_er2  -35.66465    1 73.32931 0.1383771
# habitat_ard2 -28.83582    6 69.67164 0.8616229

# For habitat, most of the weight (86%) falls within the ARD model, and has lower (better) AIC in ARD

# Let’s compute marginal ancestral states for diet under each of these two models as follows.

# ER model:
diet_er2.anc<-ancr(diet_er2)
diet_er2.anc

# Log-likelihood = -24.526861

# ARD model:
diet_ard2.anc<-ancr(diet_ard2)
diet_ard2.anc

# Log-likelihood = -23.302427 

# Log-likelihood for ARD is higher than ER for diet

# Let’s compute marginal ancestral states for habitat under each of these two models as follows.

# ER model:
habitat_er2.anc<-ancr(habitat_er2)
habitat_er2.anc

# Log-likelihood = -35.664653 

# ARD model:
habitat_ard2.anc<-ancr(habitat_ard2)
habitat_ard2.anc

# Log-likelihood = -28.835819 

# Trying to use fithrm to make an ordered model for diet
ordered_Mk_diet2<-fitHRM(test_tree,diet_data2,niter=10,
                        parallel=TRUE,ncores=10,umbral=TRUE,ordered=TRUE,
                        order=sort(unique(diet_data2)),pi="fitzjohn",ncat=1,
                        logscale=TRUE)
ordered_Mk_diet2

# Fitted (or set) value of Q:
#  0         1          2
# 0 -0.038655  0.038655 0.000000
# 1  0.000000 -0.042175 0.042175
# 2  0.000000  0.000000 0.000000

# Fitted (or set) value of pi:
# 0 1 2 
# 1 0 0 

# Log-likelihood: -24.660978 

# Diet: Comparing ordered MK model with fitzjohn prior with ER model for diet and ARD model for diet

anova(ordered_Mk_diet2,diet_er2,diet_ard2)

#            log(L) d.f.      AIC    weight
# object    -24.66098    4 57.32196 0.04082486
# diet_er2  -24.52686    1 51.05372 0.93767965
# diet_ard2 -23.30243    6 58.60485 0.02149548

# AIC is lowest in ER, weight is mostly carried in the ER model, log likelihood is slightly better in ARD

# Look at directionality of transitions as well as transition rate
plot(as.Qmatrix(ordered_Mk_diet2),width=TRUE,color=TRUE,xlim=c(-1.5,1),
     text=FALSE,max.lwd=4)

# Let’s compute marginal ancestral states for diet under the ordered mk model as follows:
diet_orderedmk2.anc<-ancr(ordered_Mk_diet2)
diet_orderedmk2.anc


# Let's plot the diet ordered mk model now
cols<-viridisLite::viridis(n=3)

er_cex<-apply(diet_orderedmk2.anc$ace,1,function(x) if(any(x>0.95)) 0.2 else 0.6)
plot(diet_orderedmk2.anc,args.nodelabels=list(cex=er_cex,piecol=cols),
     mar=c(0.1,0.1,1.1,0.1))
mtext("a) diet ordered mk model marginal ancestral states",line=0,adj=0)



# Trying to use fithrm to make an ordered model for habitat
ordered_Mk_hab2<-fitHRM(test_tree,habitat_data2,niter=10,
                       parallel=TRUE,ncores=10,umbral=TRUE,ordered=TRUE,
                       order=sort(unique(habitat_data2)),pi="fitzjohn",ncat=1,
                       logscale=TRUE)
ordered_Mk_hab2

# Fitted (or set) value of Q:
# 0         1          2
# 0 -0.190151  0.190151 0.000000
# 1  0.187264 -0.197246 0.009982
# 2  0.000000  0.000000 0.000000

# Fitted (or set) value of pi:
# 0        1        2 
# 0.5235 0.4765 0.0000 

# Log-likelihood: -28.835929 

# Habitat: Comparing ordered MK model with fitzjohn prior with ER model for habitat and ARD model for habitat
anova(ordered_Mk_hab2,habitat_er2,habitat_ard2)

# log(L) d.f.      AIC     weight
# object       -28.83593    4 65.67186 0.86423879
# habitat_er2  -35.66465    1 73.32931 0.01878625
# habitat_ard2 -28.83582    6 69.67164 0.11697496

# Most of the weight carried in ordered MK, ordered MK and ARD have comparable log likelihoods and ordered MK has the best AIC

# See the directions of the transitions
plot(as.Qmatrix(ordered_Mk_hab2),width=TRUE,color=TRUE,xlim=c(-1.5,1),
     text=FALSE,max.lwd=4)

# Let’s compute marginal ancestral states for habitat under the ordered mk model as follows:
hab_orderedmk2.anc<-ancr(ordered_Mk_hab2)
hab_orderedmk2.anc

# Let's plot the habitat ordered mk model now
cols<-viridisLite::viridis(n=3)

er_cex<-apply(hab_orderedmk2.anc$ace,1,function(x) if(any(x>0.95)) 0.2 else 0.6)
plot(hab_orderedmk2.anc,args.nodelabels=list(cex=er_cex,piecol=cols),
     mar=c(0.1,0.1,1.1,0.1))
mtext("a) habitat ordered mk model marginal ancestral states",line=0,adj=0)


# Trying ordered mk with habitat but going to use pi=estimated to see if it changes the root and to compare to other models
ordered_Mk_hab_est2<-fitHRM(test_tree,habitat_data2,niter=10,
                           parallel=TRUE,ncores=10,umbral=TRUE,ordered=TRUE,
                           order=sort(unique(habitat_data2)),pi="estimated",ncat=1,
                           logscale=TRUE)
ordered_Mk_hab_est2

# Fitted (or set) value of Q:
#  0         1          2
# 0 -0.174700  0.174700  0.000000
# 1  0.189341 -0.206459  0.017118
# 2  0.000000  0.281659 -0.281659

# Fitted (or set) value of pi:
# 0 1 2 
# 0 0 1

# Log-likelihood: -29.338525 


# Going to compare this to the other habitat models
anova(ordered_Mk_hab_est2,ordered_Mk_hab2,habitat_er2,habitat_ard2)

# log(L) d.f.      AIC      weight
# object          -29.33853    4 66.67705 0.34332715
# ordered_Mk_hab2 -28.83593    4 65.67186 0.56752215
# habitat_er2     -35.66465    1 73.32931 0.01233642
# habitat_ard2    -28.83582    6 69.67164 0.07681428

# Weight is carried mostly in ordered MK with fitzjohn prior, AIC is lowest in ordered MK with fitzjohn prior

# Going to do ordered MK with estimated prior for diet now
ordered_Mk_diet_est2<-fitHRM(test_tree,diet_data2,niter=10,
                            parallel=TRUE,ncores=10,umbral=TRUE,ordered=TRUE,
                            order=sort(unique(diet_data2)),pi="estimated",ncat=1,
                            logscale=TRUE)
ordered_Mk_diet_est2

# Fitted (or set) value of Q:
# 0         1          2
# 0 -0.006257  0.006257  0.000000
# 1  0.065829 -0.118682  0.052852
# 2  0.000000  0.335199 -0.335199

# Fitted (or set) value of pi:
# 0 1 2 
# 0 0 1

# Log-likelihood: -26.436962 

# Going to now comapre this to the other diet models
anova(ordered_Mk_diet_est2,ordered_Mk_diet2,diet_er2,diet_ard2)

# log(L) d.f.      AIC     weight
# object           -26.43696    4 60.87392 0.006864887
# ordered_Mk_diet2 -24.66098    4 57.32196 0.040544606
# diet_er2         -24.52686    1 51.05372 0.931242589
# diet_ard2        -23.30243    6 58.60485 0.021347918

# Log-likelihood was best under ARD, AIC was best for MK with Fitzjohn, weight mostly carried in ER

# See the directions of the transitions for habitat with estimated prior
plot(as.Qmatrix(ordered_Mk_hab_est2),width=TRUE,color=TRUE,xlim=c(-1.5,1),
     text=FALSE,max.lwd=4)

# Let’s compute marginal ancestral states for habitat under estimated prior under the ordered mk model as follows:

hab_est_orderedmk2.anc<-ancr(ordered_Mk_hab_est2)
hab_est_orderedmk2.anc

# Let's plot the diet ordered mk model with fitzjohn prior now
cols<-viridisLite::viridis(n=3)

er_cex<-apply(hab_est_orderedmk2.anc$ace,1,function(x) if(any(x>0.95)) 0.2 else 0.6)
plot(hab_est_orderedmk2.anc,args.nodelabels=list(cex=er_cex,piecol=cols),
     mar=c(0.1,0.1,1.1,0.1))
mtext("a) habitat ordered mk model marginal ancestral states",line=0,adj=0)

# Trying fithrm ordered MK with habitat but going to create my own pi to see if it changes the root

#### This is how you specify your own prior state for pi!!
# Define a custom prior distribution function for pi
custom_pi_prior <- function(num_states) {
  # Set one state (e.g., state 0) to have a known value (e.g., 1)
  # Distribute the remaining probability mass among other states
  prior_probs <- rep(0, num_states)
  prior_probs[1] <- 1  # Set state 0 to have a known value
  if(num_states > 1) {
    remaining_prob_mass <- 1 - prior_probs[1]
    remaining_prob_per_state <- remaining_prob_mass / (num_states - 1)
    prior_probs[-1] <- remaining_prob_per_state
  }
  return(prior_probs)
}


# Specify the custom prior distribution for pi
pi_prior <- custom_pi_prior(num_states = 3)

# Doing the ordered mk for habitat with custom root prior 
ordered_Mk_hab_prior2<-fitHRM(test_tree,habitat_data2,niter=10,
                             parallel=TRUE,ncores=10,umbral=TRUE,ordered=TRUE,
                             order=sort(unique(habitat_data2)),pi=pi_prior,ncat=1,
                             logscale=TRUE)
ordered_Mk_hab_prior2

# Fitted (or set) value of Q:
#  0         1          2
# 0 -0.187695  0.187695 0.000000
# 1  0.185302 -0.195790 0.010488
# 2  0.000000  0.000000 0.000000

# Fitted (or set) value of pi:
# 0 1 2 
# 1 0 0 

# Log-likelihood: -28.790819 

# Going to compare all habitat models done so far
anova(ordered_Mk_hab_prior2,ordered_Mk_hab_est2,ordered_Mk_hab2,habitat_er2,habitat_ard2)

# log(L) d.f.      AIC      weight
# object              -28.79082    4 65.58164 0.372533097
# ordered_Mk_hab_est2 -29.33853    4 66.67705 0.215426423
# ordered_Mk_hab2     -28.83593    4 65.67186 0.356101365
# habitat_er2         -35.66465    1 73.32931 0.007740695
# habitat_ard2        -28.83582    6 69.67164 0.048198420

# Log-likelihood comparable for all except ER, AIC comparable for all ordered mk models, and weight mostly split between all ordered mk model with custom prior and ordered mk with fitzjohn

# See the directions of the transitions
dev.off()
plot(as.Qmatrix(ordered_Mk_hab_prior2),width=TRUE,color=TRUE,xlim=c(-1.5,1),
     text=FALSE,max.lwd=4)

# Let’s compute marginal ancestral states for habitat under custom prior under the ordered mk model as follows:

hab_prior_orderedmk2.anc<-ancr(ordered_Mk_hab_prior2)
hab_prior_orderedmk2.anc

# Let's plot the habitat ordered mk model with custom prior now
cols<-viridisLite::viridis(n=3)

er_cex<-apply(hab_prior_orderedmk2.anc$ace,1,function(x) if(any(x>0.95)) 0.2 else 0.6)
plot(hab_prior_orderedmk2.anc,args.nodelabels=list(cex=er_cex,piecol=cols),
     mar=c(0.1,0.1,1.1,0.1))
mtext("a) habitat ordered mk model with custom prior marginal ancestral states",line=0,adj=0)

# Doing the ordered mk for diet with custom root prior 
ordered_Mk_diet_prior2<-fitHRM(test_tree,diet_data2,niter=10,
                              parallel=TRUE,ncores=10,umbral=TRUE,ordered=TRUE,
                              order=sort(unique(diet_data2)),pi=pi_prior,ncat=1,
                              logscale=TRUE)
ordered_Mk_diet_prior2

# Fitted (or set) value of Q:
# 0         1          2
# 0 -0.038655  0.038655 0.000000
# 1  0.000000 -0.042175 0.042175
# 2  0.000000  0.000000 0.000000

# Fitted (or set) value of pi:
# 0 1 2 
# 1 0 0 

# Log-likelihood: -24.660978 

# Going to compare all diet models done so far
anova(ordered_Mk_diet_prior2,ordered_Mk_diet_est2,ordered_Mk_diet2,diet_er2,diet_ard2)

#                       log(L) d.f.      AIC     weight
# object               -24.66098    4 57.32196 0.038964793
# ordered_Mk_diet_est2 -26.43696    4 60.87392 0.006597398
# ordered_Mk_diet2     -24.66098    4 57.32196 0.038964793
# diet_er2             -24.52686    1 51.05372 0.894956914
# diet_ard2            -23.30243    6 58.60485 0.020516101

# Log-likelihood best in ARD with fitzjohn. AIC best in ER with fitzjohn. Weight carried mostly in ER with fitzjohn

# See the directions of the transitions
plot(as.Qmatrix(ordered_Mk_diet_prior2),width=TRUE,color=TRUE,xlim=c(-1.5,1),
     text=FALSE,max.lwd=4)

# Let’s compute marginal ancestral states for diet under custom prior under the ordered mk model as follows:

diet_prior_orderedmk2.anc<-ancr(ordered_Mk_diet_prior2)
diet_prior_orderedmk2.anc

# Let's plot the diet ordered mk model with custom prior now
cols<-viridisLite::viridis(n=3)

er_cex<-apply(diet_prior_orderedmk2.anc$ace,1,function(x) if(any(x>0.95)) 0.2 else 0.6)
plot(diet_prior_orderedmk2.anc,args.nodelabels=list(cex=er_cex,piecol=cols),
     mar=c(0.1,0.1,1.1,0.1))
mtext("a) diet ordered mk model marginal ancestral states",line=0,adj=0)

### Now going to do ARD and Er models with custom prior and compare everything/plot everything
# First doing habitat ARD model
habitat_ard_custom2<-fitMk(test_tree,habitat_data2,model="ARD",
                          pi=pi_prior)
habitat_ard_custom2

# Fitted (or set) value of Q:
#  0          1          2
#  0 -0.187695  0.187695 0.000000
# 1  0.185302 -0.195790 0.010488
# 2  0.000000  0.000000 0.000000

# Fitted (or set) value of pi:
# 0 1 2 
# 1 0 0 

# Log-likelihood: -28.790819 

# See the directions of the transitions
plot(as.Qmatrix(habitat_ard_custom2),width=TRUE,color=TRUE,xlim=c(-1.5,1),
     text=FALSE,max.lwd=4)

# Let’s compute marginal ancestral states for habitat under custom prior under the ARD model as follows:

habitat_ard_custom2.anc<-ancr(habitat_ard_custom2)
habitat_ard_custom2.anc

# Let's plot the habitat ARD model with custom prior now
cols<-viridisLite::viridis(n=3)

er_cex<-apply(habitat_ard_custom2.anc$ace,1,function(x) if(any(x>0.95)) 0.2 else 0.6)
plot(habitat_ard_custom2.anc,args.nodelabels=list(cex=er_cex,piecol=cols),
     mar=c(0.1,0.1,1.1,0.1))
mtext("a) diet ard model with custom prior marginal ancestral states",line=0,adj=0)

# Now going to do habitat with ER model
habitat_er_custom2<-fitMk(test_tree,habitat_data2,model="ER",
                         pi=pi_prior)
habitat_er_custom2
# Fitted (or set) value of Q:
# 0          1          2
# 0 -0.154942  0.077471  0.077471
# 1  0.077471 -0.154942  0.077471
# 2  0.077471  0.077471 -0.154942

# Fitted (or set) value of pi:
# 0 1 2 
# 1 0 0 

# Log-likelihood: -35.59382 

# See the directions of the transitions
plot(as.Qmatrix(habitat_er_custom2),width=TRUE,color=TRUE,xlim=c(-1.5,1),
     text=FALSE,max.lwd=4)

# Let’s compute marginal ancestral states for habitat under custom prior under the ER model as follows:

habitat_er_custom2.anc<-ancr(habitat_er_custom2)
habitat_er_custom2.anc

# Let's plot the habitat ER model with custom prior now
dev.off()
cols<-viridisLite::viridis(n=3)

er_cex<-apply(habitat_er_custom2.anc$ace,1,function(x) if(any(x>0.95)) 0.2 else 0.6)
plot(habitat_er_custom2.anc,args.nodelabels=list(cex=er_cex,piecol=cols),
     mar=c(0.1,0.1,1.1,0.1))
mtext("a) habitat er model with custom prior marginal ancestral states",line=0,adj=0)

# Now going to compare stats from all the models for habitat
anova(habitat_ard_custom2,habitat_er_custom2,ordered_Mk_hab_prior2,ordered_Mk_hab_est2,ordered_Mk_hab2,habitat_er2,habitat_ard2)

# log(L) d.f.      AIC      weight
# habitat_ard_custom2   -28.79082    6 69.58164 0.047620332
# habitat_er_custom2    -35.59382    1 73.18764 0.007848001
# ordered_Mk_hab_prior2 -28.79082    4 65.58164 0.351869307
# ordered_Mk_hab_est2   -29.33853    4 66.67705 0.203477078
# ordered_Mk_hab2       -28.83593    4 65.67186 0.336349016
# habitat_er2           -35.66465    1 73.32931 0.007311331
# habitat_ard2          -28.83582    6 69.67164 0.045524934

# Best log-likelihood: all except any ERs. Best AICs: all ordered MKs. Best weight: ordered mk with custom prior, then ordered mk with fitzjohn

# Now doing diet ARD model with custom prior
diet_ard_custom2<-fitMk(test_tree,diet_data2,model="ARD",
                       pi=pi_prior)
diet_ard_custom2

# Fitted (or set) value of Q:
# 0         1          2
# 0 -0.050471  0.035145  0.015326
# 1  0.000000 -0.037591  0.037591
# 2  0.132919  0.000000 -0.132919

# Fitted (or set) value of pi:
# 0 1 2 
# 1 0 0 

# Log-likelihood: -23.096062

# See the directions of the transitions
plot(as.Qmatrix(diet_ard_custom2),width=TRUE,color=TRUE,xlim=c(-1.5,1),
     text=FALSE,max.lwd=4)

# Let’s compute marginal ancestral states for diet under custom prior under the ARD model as follows:

diet_ard_custom2.anc<-ancr(diet_ard_custom2)
diet_ard_custom2.anc

# Let's plot the diet ARD model with custom prior now
dev.off()
cols<-viridisLite::viridis(n=3)

er_cex<-apply(diet_ard_custom2.anc$ace,1,function(x) if(any(x>0.95)) 0.2 else 0.6)
plot(diet_ard_custom2.anc,args.nodelabels=list(cex=er_cex,piecol=cols),
     mar=c(0.1,0.1,1.1,0.1))
mtext("a) diet ard model with custom prior marginal ancestral states",line=0,adj=0)

########## Now doing diet ER model with custom prior. This is the best one
diet_er_custom2<-fitMk(test_tree,diet_data2,model="ER",
                      pi=pi_prior)
diet_er_custom2

# Fitted (or set) value of Q:
# 0         1         2
# 0 -0.043019  0.021510  0.021510
# 1  0.021510 -0.043019  0.021510
# 2  0.021510  0.021510 -0.043019

# Fitted (or set) value of pi:
# 0 1 2 
# 1 0 0 

# Log-likelihood: -24.468464 

# See the directions of the transitions
plot(as.Qmatrix(diet_er_custom2),width=TRUE,color=TRUE,xlim=c(-1.5,1),
     text=FALSE,max.lwd=4)

# Let’s compute marginal ancestral states for diet under custom prior under the ER model as follows:

diet_er_custom2.anc<-ancr(diet_er_custom2)
diet_er_custom2.anc

# Let's plot the diet ER model with custom prior now
dev.off()
cols<-viridisLite::viridis(n=3)

er_cex<-apply(diet_er_custom2.anc$ace,1,function(x) if(any(x>0.95)) 0.2 else 0.6)
plot(diet_er_custom2.anc,args.nodelabels=list(cex=er_cex,piecol=cols),
     mar=c(0.1,0.1,1.1,0.1))
mtext("a) diet er model with custom prior marginal ancestral states",line=0,adj=0)

# Now going to compare all models for diet
anova(diet_ard_custom2,diet_er_custom2,ordered_Mk_diet_prior2,ordered_Mk_diet_est2,ordered_Mk_diet2,diet_er2,diet_ard2)

#                          log(L) d.f.      AIC     weight
# diet_ard_custom2       -23.09606    6 58.19212 0.012775333
# diet_er_custom2        -24.46846    1 50.93693 0.480637588
# ordered_Mk_diet_prior2 -24.66098    4 57.32196 0.019739061
# ordered_Mk_diet_est2   -26.43696    4 60.87392 0.003342157
# ordered_Mk_diet2       -24.66098    4 57.32196 0.019739061
# diet_er2               -24.52686    1 51.05372 0.453373609
# diet_ard2              -23.30243    6 58.60485 0.010393192

# Best log-likelihoods: ARD with custom prior, ARD with fitzjohn.
# Best AICs: ER with custom prior, ER with fitzjohn
# Best weight: ER with custom prior, ER with fitzjohn

# Going to do symmetric model for diet with custom prior
diet_sym_custom2<-fitMk(test_tree,diet_data2,model="SYM",
                       pi=pi_prior)
diet_sym_custom2

# Fitted (or set) value of Q:
#  0         1         2
# 0 -0.033359  0.021703  0.011656
# 1  0.021703 -0.051370  0.029667
# 2  0.011656  0.029667 -0.041323

# Fitted (or set) value of pi:
# 0 1 2 
# 1 0 0 

# Log-likelihood: -24.214784 

# See the directions of the transitions
plot(as.Qmatrix(diet_sym_custom2),width=TRUE,color=TRUE,xlim=c(-1.5,1),
     text=FALSE,max.lwd=4)

# Let’s compute marginal ancestral states for diet under custom prior under the ER model as follows:

diet_sym_custom2.anc<-ancr(diet_sym_custom2)
diet_sym_custom2.anc

# Let's plot the diet ER model with custom prior now
dev.off()
cols<-viridisLite::viridis(n=3)

er_cex<-apply(diet_sym_custom2.anc$ace,1,function(x) if(any(x>0.95)) 0.2 else 0.6)
plot(diet_sym_custom2.anc,args.nodelabels=list(cex=er_cex,piecol=cols),
     mar=c(0.1,0.1,1.1,0.1))
mtext("a) diet symmetric model with custom prior marginal ancestral states",line=0,adj=0)

# Now going to compare all models for diet
anova(diet_sym_custom2,diet_ard_custom2,diet_er_custom2,ordered_Mk_diet_prior2,ordered_Mk_diet_est2,ordered_Mk_diet2,diet_er2,diet_ard2)

#                           log(L) d.f.      AIC      weight
# diet_sym_custom2       -24.21478    3 54.42957 0.077346289
# diet_ard_custom2       -23.09606    6 58.19212 0.011787209
# diet_er_custom2        -24.46846    1 50.93693 0.443462054
# ordered_Mk_diet_prior2 -24.66098    4 57.32196 0.018212318
# ordered_Mk_diet_est2   -26.43696    4 60.87392 0.003083653
# ordered_Mk_diet2       -24.66098    4 57.32196 0.018212318
# diet_er2               -24.52686    1 51.05372 0.418306843
# diet_ard2              -23.30243    6 58.60485 0.009589317

# Making transition matrices for diet and habitat:
diet_transition_matrix <- matrix(c(0, 1, 0,   # Transition from 0 to 1
                                   0, 0, 1,   # Transition from 1 to 2
                                   0, 0, 0),  # No transition from 2
                                 nrow = 3, byrow = TRUE)

# For habitat type
habitat_transition_matrix <- matrix(c(0, 1, 0,   # Transition from 0 to 1
                                      0, 0, 1,   # Transition from 1 to 2
                                      0, 0, 0),  # No transition from 2
                                    nrow = 3, byrow = TRUE)

# Now doing with model that is custom diet transition matrix and custom prior. This was the best with previous tree

diet_transmat_custom2<-fitMk(test_tree,diet_data2,model=diet_transition_matrix,
                            pi=pi_prior)
diet_transmat_custom2
diet_transition_matrix

# Fitted (or set) value of Q:
# 0         1        2
# 0 -0.040086  0.040086 0.000000
# 1  0.000000 -0.040086 0.040086
# 2  0.000000  0.000000 0.000000

# Fitted (or set) value of pi:
# 0 1 2 
# 1 0 0 

# Log-likelihood: -24.667471 

# See the directions of the transitions
plot(as.Qmatrix(diet_transmat_custom2),width=TRUE,color=TRUE,xlim=c(-1.5,1),
     text=FALSE,max.lwd=4)

# Let’s compute marginal ancestral states for diet under custom prior under custom model as follows:

diet_transmat_custom2.anc<-ancr(diet_transmat_custom2)
diet_transmat_custom2.anc

# Let's plot the diet custom model with custom prior now
dev.off()
cols<-viridisLite::viridis(n=3)

er_cex<-apply(diet_transmat_custom2.anc$ace,1,function(x) if(any(x>0.95)) 0.2 else 0.6)
plot(diet_transmat_custom2.anc,args.nodelabels=list(cex=er_cex,piecol=cols), fsize=0.6, ftype="i",
     mar=c(0.1,0.1,1.1,0.1))
mtext("diet custom Q model with custom prior marginal ancestral states", line=0,adj=0)


# Now compare this to all other diet models
anova(flatprior.fit_diet_ard2,diet_sym_custom2,diet_transmat_custom2,diet_ard_custom2,diet_er_custom2,ordered_Mk_diet_prior2,ordered_Mk_diet_est2,ordered_Mk_diet2,diet_er2,diet_ard2)

# log(L) d.f.      AIC      weight
# flatprior.fit_diet_ard2 -23.55509    6 59.11019 0.005433183
# diet_sym_custom2        -24.21478    3 54.42957 0.056420690
# diet_transmat_custom2   -24.66747    1 51.33494 0.265111141
# diet_ard_custom2        -23.09606    6 58.19212 0.008598246
# diet_er_custom2         -24.46846    1 50.93693 0.323485913
# ordered_Mk_diet_prior2  -24.66098    4 57.32196 0.013285079
# ordered_Mk_diet_est2    -26.43696    4 60.87392 0.002249388
# ordered_Mk_diet2        -24.66098    4 57.32196 0.013285079
# diet_er2                -24.52686    1 51.05372 0.305136301
# diet_ard2               -23.30243    6 58.60485 0.006994982

# Best log likelihood: ARD with custom prior and ARD with fitzjohn
# Best AIC: ER with custom prior
# Best weight: ER with custom prior

# Now doing custom habitat model with custom prior
hab_transmat_custom2<-fitMk(test_tree,habitat_data2,model=habitat_transition_matrix,
                           pi=pi_prior)
hab_transmat_custom2

# Fitted (or set) value of Q:
# 0         1        2
# 0 -0.0689  0.0689 0.0000
# 1  0.0000 -0.0689 0.0689
# 2  0.0000  0.0000 0.0000

# Fitted (or set) value of pi:
# 0 1 2 
# 1 0 0 

# Log-likelihood: -39.580177 

# See the directions of the transitions
plot(as.Qmatrix(hab_transmat_custom2),width=TRUE,color=TRUE,xlim=c(-1.5,1),
     text=FALSE,max.lwd=4)

# Let’s compute marginal ancestral states for habitat under custom prior under custom model as follows:

hab_transmat_custom2.anc<-ancr(hab_transmat_custom2)
hab_transmat_custom2.anc

# Let's plot the habitat custom model with custom prior now
dev.off()
cols<-viridisLite::viridis(n=3)

er_cex<-apply(hab_transmat_custom2.anc$ace,1,function(x) if(any(x>0.95)) 0.2 else 0.6)
plot(hab_transmat_custom2.anc,args.nodelabels=list(cex=er_cex,piecol=cols),
     mar=c(0.1,0.1,1.1,0.1))
mtext("a) habitat custom Q model with custom prior marginal ancestral states",line=0,adj=0)


# Now compare to all other habitat models
anova(hab_transmat_custom2,habitat_ard_custom2,habitat_er_custom2,ordered_Mk_hab_prior2,ordered_Mk_hab_est2,ordered_Mk_hab2,habitat_er2,habitat_ard2)

#                         log(L) d.f.      AIC       weight
# hab_transmat_custom2  -39.58018    1 81.16035 0.0001456944
# habitat_ard_custom2   -28.79082    6 69.58164 0.0476133944
# habitat_er_custom2    -35.59382    1 73.18764 0.0078468579
# ordered_Mk_hab_prior2 -28.79082    4 65.58164 0.3518180415
# ordered_Mk_hab_est2   -29.33853    4 66.67705 0.2034474326
# ordered_Mk_hab2       -28.83593    4 65.67186 0.3363000118
# habitat_er2           -35.66465    1 73.32931 0.0073102661
# habitat_ard2          -28.83582    6 69.67164 0.0455183013

# Habitat with custom Q matrix with unidirectional transitions.. is bad

# Going to change the Q so strong transition from 0 to 1, weaker transition from 1 to 2
hab_adj_q <- matrix(c(0, 2, 0,   # Strong Transition from 0 to 1
                      0, 0, 1,   # Weaker Transition from 1 to 2
                      0, 0, 0),  # No transition from 2
                    nrow = 3, byrow = TRUE)

# Now doing custom habitat model (strong transition from 0-1) with custom prior
hab_transmat_strong1_custom2<-fitMk(test_tree,habitat_data2,model=hab_adj_q,
                                   pi=pi_prior)
hab_transmat_strong1_custom2

# Fitted (or set) value of Q:
# 0         1        2
# 0 -0.078443  0.078443 0.000000
# 1  0.000000 -0.029424 0.029424
# 2  0.000000  0.000000 0.000000

# Fitted (or set) value of pi:
# 0 1 2 
# 1 0 0 

# Log-likelihood: -38.999199 

# See the directions of the transitions
plot(as.Qmatrix(hab_transmat_strong1_custom2),width=TRUE,color=TRUE,xlim=c(-1.5,1),
     text=FALSE,max.lwd=4)

# Let’s compute marginal ancestral states for habitat under custom prior under custom model as follows:

hab_transmat_strong1_custom2.anc<-ancr(hab_transmat_strong1_custom2)
hab_transmat_strong1_custom2.anc

# Let's plot the habitat custom model with custom prior now
dev.off()
cols<-viridisLite::viridis(n=3)

er_cex<-apply(hab_transmat_strong1_custom2.anc$ace,1,function(x) if(any(x>0.95)) 0.2 else 0.6)
plot(hab_transmat_strong1_custom2.anc,args.nodelabels=list(cex=er_cex,piecol=cols),
     mar=c(0.1,0.1,1.1,0.1))
mtext("a) habitat custom Q model with custom prior marginal ancestral states",line=0,adj=0)


# Now compare to all other habitat models
anova(hab_transmat_strong1_custom2,hab_transmat_custom2,habitat_ard_custom2,habitat_er_custom2,ordered_Mk_hab_prior2,ordered_Mk_hab_est2,ordered_Mk_hab2,habitat_er2,habitat_ard2)

#                                log(L) d.f.      AIC       weight
# hab_transmat_strong1_custom2 -38.99920    2 81.99840 9.581251e-05
# hab_transmat_custom2         -39.58018    1 81.16035 1.456805e-04
# habitat_ard_custom2          -28.79082    6 69.58164 4.760883e-02
# habitat_er_custom2           -35.59382    1 73.18764 7.846106e-03
# ordered_Mk_hab_prior2        -28.79082    4 65.58164 3.517843e-01
# ordered_Mk_hab_est2          -29.33853    4 66.67705 2.034279e-01
# ordered_Mk_hab2              -28.83593    4 65.67186 3.362678e-01
# habitat_er2                  -35.66465    1 73.32931 7.309566e-03
# habitat_ard2                 -28.83582    6 69.67164 4.551394e-02

# still pretty bad

# Going to change the Q to have strong transition from 0 to 1, regular transition from 1 to 2, and weak reverse transition from 2 to 1
hab_reverse_q <- matrix(c(1, 2, 0,   # strong Transition from 0 to 1, regular transition from 0 to 0
                          0, 1, 1,   # regular Transition from 1 to 2, regular transition form 1 to 1
                          0, 1, 0),  # No transition from 2, regular transition from 2 to 1
                        nrow = 3, byrow = TRUE)

# Now doing custom habitat model (strong transition from 0-1) with custom prior
hab_transmat_reverse_custom2<-fitMk(test_tree,habitat_data2,model=hab_reverse_q,
                                   pi=pi_prior)
hab_transmat_reverse_custom2

# Fitted (or set) value of Q:
# 0         1        2
# 0 -0.078467  0.078467  0.000000
# 1  0.000000 -0.027366  0.027366
# 2  0.000000  0.027366 -0.027366

# Fitted (or set) value of pi:
# 0 1 2 
# 1 0 0 

# Log-likelihood: -39.109056 

# See the directions of the transitions
plot(as.Qmatrix(hab_transmat_reverse_custom2),width=TRUE,color=TRUE,xlim=c(-1.5,1),
     text=FALSE,max.lwd=4)

# Let’s compute marginal ancestral states for habitat under custom prior under custom model as follows:

hab_transmat_reverse_custom2.anc<-ancr(hab_transmat_reverse_custom2)
hab_transmat_reverse_custom2.anc

# Let's plot the habitat custom model with custom prior now
dev.off()
cols<-viridisLite::viridis(n=3)

er_cex<-apply(hab_transmat_reverse_custom2.anc$ace,1,function(x) if(any(x>0.95)) 0.2 else 0.6)
plot(hab_transmat_reverse_custom2.anc,args.nodelabels=list(cex=er_cex,piecol=cols),
     mar=c(0.1,0.1,1.1,0.1))
mtext("a) habitat custom Q model with custom prior marginal ancestral states",line=0,adj=0)



# Now compare to all other habitat models
anova(hab_transmat_reverse_custom2,hab_transmat_strong1_custom2,hab_transmat_custom2,habitat_ard_custom2,habitat_er_custom2,ordered_Mk_hab_prior2,ordered_Mk_hab_est2,ordered_Mk_hab2,habitat_er2,habitat_ard2)

#                                log(L) d.f.      AIC       weight
# hab_transmat_reverse_custom2 -39.10906    2 82.21811 8.583703e-05
# hab_transmat_strong1_custom2 -38.99920    2 81.99840 9.580429e-05
# hab_transmat_custom2         -39.58018    1 81.16035 1.456680e-04
# habitat_ard_custom2          -28.79082    6 69.58164 4.760475e-02
# habitat_er_custom2           -35.59382    1 73.18764 7.845433e-03
# ordered_Mk_hab_prior2        -28.79082    4 65.58164 3.517541e-01
# ordered_Mk_hab_est2          -29.33853    4 66.67705 2.034105e-01
# ordered_Mk_hab2              -28.83593    4 65.67186 3.362389e-01
# habitat_er2                  -35.66465    1 73.32931 7.308938e-03
# habitat_ard2                 -28.83582    6 69.67164 4.551003e-02

# Still bad

########## Now going to do a symmetrical model with custom prior. This one was the best with previous tree and is the best with this one too
hab_sym_customprior2<-fitMk(test_tree,habitat_data2,model="SYM",
                           pi=pi_prior)
hab_sym_customprior2

# Fitted (or set) value of Q:
# 0          1         2
# 0 -0.185447  0.185447  0.000000
# 1  0.185447 -0.196165  0.010718
# 2  0.000000  0.010718 -0.010718

# Fitted (or set) value of pi:
# 0 1 2 
# 1 0 0 

# Log-likelihood: -28.799003 


library(phytools)
# See the directions of the transitions
plot(as.Qmatrix(hab_sym_customprior2),width=TRUE,color=TRUE,xlim=c(-1.5,1),
     text=FALSE,max.lwd=4)

# Let’s compute marginal ancestral states for habitat under custom prior under sym model as follows:

hab_sym_customprior2.anc<-ancr(hab_sym_customprior2)
hab_sym_customprior2.anc

# Let's plot the habitat sym model with custom prior now
dev.off()
cols<-viridisLite::viridis(n=3)


er_cex<-apply(hab_sym_customprior2.anc$ace,1,function(x) if(any(x>0.95)) 0.2 else 0.6)
plot(hab_sym_customprior2.anc,args.nodelabels=list(cex=er_cex,piecol=cols),
     mar=c(0.1,0.1,1.1,0.1), align.tip.label= TRUE, label.offset = 1, no.margin = T)
mtext("a) habitat sym model with custom prior marginal ancestral states",line=0,adj=0)


# Now compare to all other habitat models
anova(hab_sym_customprior2,hab_transmat_reverse_custom2,hab_transmat_strong1_custom2,hab_transmat_custom2,habitat_ard_custom2,habitat_er_custom2,ordered_Mk_hab_prior2,ordered_Mk_hab_est2,ordered_Mk_hab2,habitat_er2,habitat_ard2)

# log(L) d.f.      AIC       weight
# hab_sym_customprior2         -28.79900    3 63.59801 4.867513e-01
# hab_transmat_reverse_custom2 -39.10906    2 82.21811 4.405574e-05
# hab_transmat_strong1_custom2 -38.99920    2 81.99840 4.917142e-05
# hab_transmat_custom2         -39.58018    1 81.16035 7.476389e-05
# habitat_ard_custom2          -28.79082    6 69.58164 2.443307e-02
# habitat_er_custom2           -35.59382    1 73.18764 4.026658e-03
# ordered_Mk_hab_prior2        -28.79082    4 65.58164 1.805373e-01
# ordered_Mk_hab_est2          -29.33853    4 66.67705 1.044002e-01
# ordered_Mk_hab2              -28.83593    4 65.67186 1.725742e-01
# habitat_er2                  -35.66465    1 73.32931 3.751303e-03
# habitat_ard2                 -28.83582    6 69.67164 2.335796e-02

# This one has the best AIC and carries the most weight

# Going to try ordered single rate for habitat and custom prior
## ordered single rate

ordered1_model<-matrix(c(0,1,0,1,0,1,0,1,0),3,3,byrow=TRUE,
                       dimnames=list(c("a","b","c"),c("a","b","c")))
fitOrdered1<-fitMk(test_tree,habitat_data2,model=ordered1_model, pi=pi_prior)
fitOrdered1

# Fitted (or set) value of Q:
# 0         1         2
# 0 -4.063343  4.063343  0.000000
# 1  4.063343 -8.126687  4.063343
# 2  0.000000  4.063343 -4.063343

# Log-likelihood: -29.882054 

# See the directions of the transitions
plot(as.Qmatrix(fitOrdered1),width=TRUE,color=TRUE,xlim=c(-1.5,1),
     text=FALSE,max.lwd=4)

# Let’s compute marginal ancestral states for habitat under custom prior under sym model as follows:

fitOrdered1.anc<-ancr(fitOrdered1)
fitOrdered1.anc

# Let's plot the habitat sym model with custom prior now
dev.off()
cols<-viridisLite::viridis(n=3)

er_cex<-apply(fitOrdered1.anc$ace,1,function(x) if(any(x>0.95)) 0.2 else 0.6)
plot(fitOrdered1.anc,args.nodelabels=list(cex=er_cex,piecol=cols),
     mar=c(0.1,0.1,1.1,0.1))
mtext("a) habitat ordered single rate with custom prior marginal ancestral states",line=0,adj=0)


# Now compare to all other habitat models
anova(fitOrdered1,hab_sym_customprior2,hab_transmat_reverse_custom2,hab_transmat_strong1_custom2,hab_transmat_custom2,habitat_ard_custom2,habitat_er_custom2,ordered_Mk_hab_prior2,ordered_Mk_hab_est2,ordered_Mk_hab2,habitat_er2,habitat_ard2)

#                                log(L) d.f.      AIC       weight
# fitOrdered1                 -29.88205    1 61.76411 0.0390159892
# hab_sym_customprior         -25.94129    3 57.88257 0.2717113750
# hab_transmat_reverse_custom -33.88971    2 71.77942 0.0002608831
# hab_transmat_strong1_custom -33.91455    2 71.82911 0.0002544811
# hab_transmat_custom         -34.18938    1 70.37876 0.0005255267
# habitat_ard_custom          -25.16453    6 62.32906 0.0294147145
# habitat_er_custom           -32.35694    1 66.71387 0.0032840941
# ordered_Mk_hab_prior        -25.20797    4 58.41593 0.2081084213
# ordered_Mk_hab_est          -25.20968    4 58.41936 0.2077520209
# ordered_Mk_hab              -25.21321    4 58.42642 0.2070202848
# habitat_er                  -32.36079    1 66.72157 0.0032714727
# habitat_ard                 -25.16569    6 62.33137 0.0293807365

# Sym custom prior is still the best

# Now going to do an ordered symmetric model with custom prior
## ordered symmetric model
model_sym<-matrix(c(0,1,0,2,0,1,0,2,0),3,3,byrow=TRUE,
                  dimnames=list(c("a","b","c"),c("a","b","c")))
fitOrdered_sym<-fitMk(test_tree,habitat_data2,model=model_sym, pi=pi_prior)
fitOrdered_sym

# Fitted (or set) value of Q:
# 0          1          2
# 0 -6.190584   6.190584   0.000000
# 1 14.885980 -21.076564   6.190584
# 2  0.000000  14.885980 -14.885980

# Log-likelihood: -27.765121 

# See the directions of the transitions
plot(as.Qmatrix(fitOrdered_sym),width=TRUE,color=TRUE,xlim=c(-1.5,1),
     text=FALSE,max.lwd=4)

# Let’s compute marginal ancestral states for habitat under custom prior under ordered sym model as follows:

fitOrdered_sym.anc<-ancr(fitOrdered_sym)
fitOrdered_sym.anc

# Let's plot the habitat sym model with custom prior now
dev.off()
cols<-viridisLite::viridis(n=3)

er_cex<-apply(fitOrdered_sym.anc$ace,1,function(x) if(any(x>0.95)) 0.2 else 0.6)
plot(fitOrdered_sym.anc,args.nodelabels=list(cex=er_cex,piecol=cols),
     mar=c(0.1,0.1,1.1,0.1))
mtext("a) habitat ordered symmetrical rate with custom prior marginal ancestral states",line=0,adj=0)


# Now compare to all other habitat models
anova(fitOrdered_sym,fitOrdered1,hab_sym_customprior2,hab_transmat_reverse_custom2,hab_transmat_strong1_custom2,hab_transmat_custom2,habitat_ard_custom2,habitat_er_custom2,ordered_Mk_hab_prior2,ordered_Mk_hab_est2,ordered_Mk_hab2,habitat_er2,habitat_ard2)

# log(L) d.f.      AIC       weight
# fitOrdered_sym              -27.76512    2 59.53024 0.1065144133
# fitOrdered1                 -29.88205    1 61.76411 0.0348602240
# hab_sym_customprior         -25.94129    3 57.88257 0.2427701973
# hab_transmat_reverse_custom -33.88971    2 71.77942 0.0002330953
# hab_transmat_strong1_custom -33.91455    2 71.82911 0.0002273752
# hab_transmat_custom         -34.18938    1 70.37876 0.0004695505
# habitat_ard_custom          -25.16453    6 62.32906 0.0262816234
# habitat_er_custom           -32.35694    1 66.71387 0.0029342908
# ordered_Mk_hab_prior        -25.20797    4 58.41593 0.1859418749
# ordered_Mk_hab_est          -25.20968    4 58.41936 0.1856234363
# ordered_Mk_hab              -25.21321    4 58.42642 0.1849696406
# habitat_er                  -32.36079    1 66.72157 0.0029230137
# habitat_ard                 -25.16569    6 62.33137 0.0262512646


# Let's plot the habitat sym model with custom prior now

# Now going to do a symmetrical model with custom prior with new tree
hab_sym_customprior2<-fitMk(test_tree,habitat_data2,model="SYM",
                            pi=pi_prior)
hab_sym_customprior2

# Fitted (or set) value of Q:
#  0         1         2
# 0 -0.185447  0.185447  0.000000
# 1  0.185447 -0.196165  0.010718
# 2  0.000000  0.010718 -0.010718

# Fitted (or set) value of pi:
#  0 1 2 
# 1 0 0 

# Log-likelihood: -28.799003 


library(phytools)
# See the directions of the transitions
plot(as.Qmatrix(hab_sym_customprior2),width=TRUE,color=TRUE,xlim=c(-1.5,1),
     text=FALSE,max.lwd=4)

# Let’s compute marginal ancestral states for habitat under custom prior under sym model as follows:

hab_sym_customprior2.anc<-ancr(hab_sym_customprior2)
hab_sym_customprior2.anc

# Let's plot the habitat sym model with custom prior now
dev.off()
cols<-viridisLite::viridis(n=3)


er_cex<-apply(hab_sym_customprior2.anc$ace,1,function(x) if(any(x>0.95)) 0.2 else 0.6)
plot(hab_sym_customprior2.anc,args.nodelabels=list(cex=er_cex,piecol=cols),
     mar=c(0.1,0.1,1.1,0.1), align.tip.label= TRUE, label.offset = 1, no.margin = T)
mtext("a) habitat sym model with custom prior marginal ancestral states",line=0,adj=0)



# Now hab_sym_customprior.anc contains the aligned tree structure
dev.off()
cols<-viridisLite::viridis(n=3)


er_cex<-apply(hab_sym_customprior.anc$ace,1,function(x) if(any(x>0.95)) 0.2 else 0.6)
plot(hab_sym_customprior2.anc,args.nodelabels=list(cex=er_cex,piecol=cols),
     mar=c(0.1,0.1,1.1,0.1), align.tip.label= TRUE, label.offset = 1, no.margin = T)
mtext("Marginal ancestral states for symmetric model of synanthropic evolution with custom prior",line=0,adj=0)



##### Making the figure legend

# Best diet model: ER with custom prior
dev.off()
library(viridisLite)
cols <- plasma(n = 3)
#cols<-viridisLite::viridis(n=3)

er_cex<-apply(diet_er_custom2.anc$ace,1,function(x) if(any(x>0.95)) 0.2 else 0.6)
plot(diet_er_custom2.anc,args.nodelabels=list(cex=er_cex,piecol=cols),
     mar=c(0.1,0.1,1.1,0.1))
mtext("a) diet er model with custom prior marginal ancestral states",line=0,adj=0)

legend("bottomleft", legend = c("Non-keratinophagous", "Facultatively keratinophagous", "Obligately keratinophagous"),
       fill = cols, title = "Dietary Strategy",cex=0.8)

# Let's plot the habitat sym model with custom prior now, best habitat model
dev.off()
cols<-viridisLite::viridis(n=3)


er_cex<-apply(hab_sym_customprior2.anc$ace,1,function(x) if(any(x>0.95)) 0.2 else 0.6)
plot(hab_sym_customprior2.anc,args.nodelabels=list(cex=er_cex,piecol=cols),
     mar=c(0.1,0.1,1.1,0.1), align.tip.label= TRUE, label.offset = 1, no.margin = T)
mtext("a) habitat sym model with custom prior marginal ancestral states",line=0,adj=0)


legend("bottomleft", legend = c("Non-synanthropic", "Facultatively synanthropic", "Obligately synanthropic"),
       fill = cols, title = "Synanthropy Status",cex=0.8)





