################################################################################################################################################
##                                                                                                                      
##   USE THE HDP OUTPUT FILES TO ASSESS SIGNATURES AND COMBINE THE FILES
##                                                                                                                      
##  Date: 19 AUGUST 2021                                                                                                                    
##  
##  Author: Tim Coorens adapted by Moritz Przybilla                                                                                                                    
##           
##                                                                                                                      
################################################################################################################################################

# clear workspace beforehand
rm(list = ls())

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("BiocManager", "RColorBrewer", "lsa", "lattice", "reshape2", "ggrepel", "readr", "stringr", "tidyverse", "hdp", "sigfit", "BSgenome.Hsapiens.UCSC.hg38", "nrmisc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

#####################################################################################
# READ IN DATA
#####################################################################################

# create output directory
output.dir <- "/Users/mp34/sanger/team154pc/mp34/lung/NTCU/WGS/HDP/clones"
dir.create(output.dir)

# get mutations and patient sample sheet with clusters from NDP
mutations= all_counts = read.table(paste0(output.dir, "/trinuc_mut_mat_full.txt"))
key_table=read.table(paste0(output.dir, "/key_table_full.txt"), header = T)
freq=table(key_table$Patient)

# load cosmic data from sigfit too
data("cosmic_signatures_v3.2")
ref <- t(cosmic_signatures_v3.2)

# get the variant files which are used for clustering and the final tree
hdp.list <- list.files(output.dir, pattern = "\\.Rdata", full.names = T, all.files = T)
hdp.list <- hdp.list[-grep("COSMIC|multi|hdp2refsigs|sigs_per_patient", hdp.list)]
hdp.list <- lapply(hdp.list, readRDS)

# combine objets
mut_example_multi <- hdp_multi_chain(hdp.list)
mut_example_multi

par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
p1 <- lapply(chains(mut_example_multi), plot_lik, bty="L", start=1000)
p2 <- lapply(chains(mut_example_multi), plot_numcluster, bty="L")
p3 <- lapply(chains(mut_example_multi), plot_data_assigned, bty="L")

#####################################################################################
# EXTRACT THE MUTATIONAL SIGNATURES AKA COMPONENTS
#####################################################################################

mut_example_multi <- hdp_extract_components(mut_example_multi, cos.merge = 0.95)
mut_example_multi

par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
plot_comp_size(mut_example_multi, bty="L")
trinuc_context <- sapply(strsplit(colnames(mut_count), '\\.'), `[`, 4)
group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                              each=16))

# pick your colours
mut_colours <- c(RColorBrewer::brewer.pal(10, 'Paired')[seq(1,10,2)], 'grey70')

pdf(paste0(output.dir, "/extracted_components_hdp.pdf"))
plot_comp_distn(mut_example_multi, cat_names=trinuc_context,
                grouping=group_factor, col=mut_colours,
                col_nonsig="grey80", show_group_labels=TRUE)
dev.off()

# write to file
write_rds(mut_example_multi, paste0(output.dir, "/HDP_multi_chain.Rdata"))

#####################################################################################
# CHECK THE COMPONENTS ACROSS ALL SAMPLES
#####################################################################################

pdf(paste0(output.dir, "/muts_attributed.pdf"))
plot_comp_size(mut_example_multi, bty="L")
dev.off()

trinuc_context <- sapply(strsplit(colnames(mut_count), '\\.'), `[`, 4)
group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                              each=16))
mut_colours=c("dodgerblue","black","red","grey70","olivedrab3","plum2")

# visualize the raw components
for (i in 0:mut_example_multi@numcomp){
  pdf(paste0(output.dir,"/hdp_component_",i,".pdf"),width=12,height=4)
  
  plot_comp_distn(mut_example_multi, cat_names=trinuc_context,
                  grouping=group_factor, col=mut_colours,comp=i,
                  col_nonsig="grey80", show_group_labels=TRUE)
  dev.off()
}

# investigate the signature contribution across samples
plot_dp_comp_exposure(mut_example_multi,
                      dpindices=2:34, incl_numdata_plot=FALSE,
                      col=c(RColorBrewer::brewer.pal(12, "Paired"),"magenta","firebrick"),
                      incl_nonsig=TRUE, cex.names=0.8,
                      dpnames = row.names(mutations)[2:34],
                      ylab_exp = 'Signature exposure', leg.title = 'Signature')

pdf(paste0(output.dir, "/signature_attribution.pdf"),width=10,height=8)
plot_dp_comp_exposure(mut_example_multi, dpindices=5:length(mut_example_multi@comp_dp_counts), incl_nonsig = T,ylab_exp = 'Signature exposure', leg.title = 'Signature',
                      col=c(RColorBrewer::brewer.pal(12, "Set3"),"magenta","firebrick",RColorBrewer::brewer.pal(8, "Dark2")))
dev.off()

#####################################################################################
# PERFORM SIGNATURE DECONVOLUTION
#####################################################################################

options(stringsAsFactors = F)

## DECONVOLUTE HDP SIGNATURES INTO REFERENCE SIGNATURES AND ARRIVE AT FINAL SET
# extract the mean signature assignments from HDP 
dp_distn <- comp_dp_distn(mut_example_multi)
ndp <- nrow(dp_distn$mean)
ncomp <- ncol(dp_distn$mean)
mean_assignment <- t(dp_distn$mean[1:nrow(dp_distn$mean),,drop=FALSE])
write.table(mean_assignment,paste0(output.dir, "/mean_assignment_hdp.txt"))

# create a dataframe with the hdp sig estimates
mean_sigs=as.data.frame(t(comp_categ_distn(mut_example_multi)$mean))
lower=mean_sigs=as.data.frame(t(comp_categ_distn(mut_example_multi)$mean))
write.table(mean_sigs, paste0(output.dir,"/hdp_sigs.txt"))

mut.cols = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)

#Load HDP signatures
# hdp_sigs=read.table(paste0(output.dir,"/hdp_sigs.txt"))
hdp_sigs=mean_sigs

#Assess cosine similarities for all reference signatures
cosine_matrix=data.frame(matrix(nrow=ncol(hdp_sigs), ncol=ncol(ref)))
rownames(cosine_matrix)=colnames(hdp_sigs)
colnames(cosine_matrix)=colnames(ref)

# determine cosine similarity between signatures of reference and this
for (n in 1:nrow(cosine_matrix)) {
  for (m in 1:ncol(cosine_matrix)) {
    cosine_matrix[n,m] <- cosine(x=hdp_sigs[,rownames(cosine_matrix)[n]],
                                 y=ref[,colnames(cosine_matrix)[m]])
  }
}

write.table(cosine_matrix, paste0(output.dir, "/Cosine_similarities.txt"),sep="\t",quote=F)

#plot output
pdf(paste0(output.dir, "/cosine_similarities.pdf"), height=5, width=15)
color.palette = colorRampPalette(c("white", "orange", "purple"))
levelplot(t(cosine_matrix[dim(cosine_matrix)[1]:1,]),col.regions=color.palette, aspect="fill", scales=list(x=list(rot=90)))
dev.off()

colnames(hdp_sigs)=gsub("X","N",colnames(hdp_sigs))

#####################################################################################
# FIRST ITERATION - DECOMPOSING HDP SIGNATURES INTO ALL SUSPECTED SIGNATURES
#####################################################################################

#Make selection of priors sigs - alternatively, use all COSMIC sigs, but this leads to overfitting
gdsigs=c("SBS1", "SBS5", "SBS11", "SBS18", "SBS32", "SBS36")

add="" #add something to titles to differentiate multiple runs [OPTIONAL]

signatures=t(ref[,gdsigs])
sample_list=paste0("N",c(0:(ncol(hdp_sigs)-1))) 
colnames(hdp_sigs)=paste0("N",c(0:(ncol(hdp_sigs)-1))) 
profiles=hdp_sigs[,sample_list]

signature_fraction = matrix(NA,nrow=nrow(signatures),ncol=length(sample_list))
rownames(signature_fraction) = rownames(signatures)
colnames(signature_fraction) = sample_list
maxiter <- 1000

for (j in 1:length(sample_list)) {
  freqs = profiles[,j]
  freqs[is.na(freqs)] = 0
  # EM algowith to estimate the signature contribution
  alpha = runif(nrow(signatures)); alpha=alpha/sum(alpha) # Random start (seems to give ~identical results)
  for (iter in 1:maxiter) {
    contr = t(array(alpha,dim=c(nrow(signatures),96))) * t(signatures)
    probs = contr/array(rowSums(contr),dim=dim(contr))
    probs = probs * freqs
    old_alpha = alpha
    alpha = colSums(probs)/sum(probs)
    if (sum(abs(alpha-old_alpha))<1e-5) {
      break
    }
  }
  # Saving the signature contributions for the sample
  print(j/length(sample_list))
  signature_fraction[,j] = alpha
}

#Plot initial deconvolution and save results
pdf(paste0(output.dir, "/Deconvolution_hdp_sigs_R1_",add,"_priors.pdf"), height=5, width=10)
color.palette = colorRampPalette(c("white", "orange", "purple"))
levelplot((signature_fraction[nrow(signature_fraction):1,]),col.regions=color.palette, aspect="fill", scales=list(x=list(rot=90)))
dev.off()

write.table(signature_fraction, paste0(output.dir, "/hdp_known_sigs_broken_down_into_pcawg_gd_sigs_",add,"_priors.txt"), sep="\t", col.names=T, row.names = T, quote=F)

#####################################################################################
# SECOND ITERATION - DECOMPOSING THE HDP SIGNATURES WHICH DID NOT RESOLVE AGAIN 
#####################################################################################

sigs_deconv_R2=list()
for(n in 1:length(sample_list)){
  sigs_deconv_R2[[n]]=rownames(signature_fraction)[signature_fraction[,n]>0.05]
}
names(sigs_deconv_R2)=colnames(signature_fraction)

# check the signatures to go forward
sigs_to_deconv=names(sigs_deconv_R2)[unlist(lapply(sigs_deconv_R2,length))>1]

ref_sigs_R2=sort(unique(unlist(sigs_deconv_R2)))
signature_fractionR2=matrix(NA,ncol=length(sigs_to_deconv),nrow=length(ref_sigs_R2))
rownames(signature_fractionR2)=ref_sigs_R2
colnames(signature_fractionR2)=sigs_to_deconv
#repeat the deconvolution with the identified constitutive signatures
n=1
for(s in sigs_to_deconv){
  gdsigs <- sigs_deconv_R2[[s]]
  signatures <- t(ref[,gdsigs])
  
  signature_fraction = matrix(NA,nrow=nrow(signatures),ncol=length(sample_list))
  rownames(signature_fraction) = rownames(signatures)
  colnames(signature_fraction) = sample_list
  maxiter <- 1000
  
  freqs = profiles[,s]
  freqs[is.na(freqs)] = 0
  
  alpha = runif(nrow(signatures)); alpha=alpha/sum(alpha) # Random start (seems to give ~identical results)
  for (iter in 1:maxiter) {
    contr = t(array(alpha,dim=c(nrow(signatures),96))) * t(signatures)
    probs = contr/array(rowSums(contr),dim=dim(contr))
    probs = probs * freqs
    old_alpha = alpha
    alpha = colSums(probs)/sum(probs)
    if (sum(abs(alpha-old_alpha))<1e-5) {
      break
    }
  }
  # Saving the signature contributions for the sample
  signature_fractionR2[gdsigs,n] = alpha
  n=n+1
  reconsbs <- rep(0,96)
  for (g in gdsigs) {
    reconsbs=reconsbs+(ref[,g]*alpha[g])
  }
  cosine_reconst=cosine(x=reconsbs, y=hdp_sigs[,s])
  print(paste0(s,": ",cosine_reconst))
  pdf(paste0(output.dir, "/HDP_",s,"_reconstitution_",add,"_priors.pdf"))
  par(mfrow=c(length(alpha)+2,1))
  par(mar=c(1,2,4,1))
  barplot(hdp_sigs[,s], col=mut.cols, main=paste0("HDP ",s),names.arg="")
  barplot(reconsbs, col=mut.cols, main=paste0("Reconstituted ",s," cosine similarity to original: ", round(cosine_reconst, digits=2)))
  for (g in gdsigs) {
    barplot(ref[,g], col=mut.cols, main=paste0("PCAWG ", g, " accounts for ", round(alpha[g], digits=2)))
  }
  dev.off()
}

# # Component 5 only reaches a reconstructed cosine similarity of 0.6
# # thus, we will reset and include it without priors in the analysis (see below)
sigs_deconv_R2$N2="N2"
saveRDS(sigs_deconv_R2,paste0(output.dir, "/hdp2refsigs.Rdata"))

#Combine hdp signatures that did not get deconvolved and reference signatures into final table
final_sigs=cbind(ref[,ref_sigs_R2], hdp_sigs[,c("N2")])
colnames(final_sigs) <- c(colnames(ref[,ref_sigs_R2]), "N2")

#Rename the HDP components that didn't get deconvoluted
write.table(final_sigs,paste0(output.dir, "/final_sigs_2023_08_20.txt"))

#####################################################################################
# FIT SIGNATURES TO OBSERVED COUNTS 
#####################################################################################

# # read in the final signatures
# final_sigs=read.table(paste0(output.dir, "/final_sigs_2023_07_18.txt"))

# read in the trinucleotide context information for all samples
hdp_exposures=t(read.table(paste0(output.dir, "/mean_assignment_hdp.txt")))
rownames(hdp_exposures) <- rownames(mutations)[1:length(mut_example_multi@comp_dp_counts)]
hdp_sigs=read.table(paste0(output.dir, "/hdp_sigs.txt"))
colnames(hdp_sigs)=gsub("X","N",colnames(hdp_sigs))
colnames(hdp_exposures)=gsub("X","N",colnames(hdp_exposures))

# 
colnames(hdp_exposures)=paste0("N",0:(ncol(hdp_exposures)-1))
hdp2ref=readRDS(paste0(output.dir, "/hdp2refsigs.Rdata"))

# make names consistent
colnames(all_counts)=rownames(final_sigs)=colnames(cosmic_signatures_v3.2)

# AT THIS POINT WE CAN ADD ANY COMPONENT TO THE FINAL SIGS IF IT HAS NOT PROPERLY RESOLVED
# This enables us to do "de-novo" extraction
# Combine hdp signatures that did not get deconvolved and reference signatures into final table
# Here, we include component 5 given that it only has a cosine similarity of 0.6

# combine the reference signatures with the inferred signatures
colnames(all_counts)=rownames(final_sigs)=colnames(cosmic_signatures_v3.2)
# final_sigs=cbind(final_sigs,hdp_sigs)

# get the patients of interest
patients=unique(substr(rownames(all_counts),1,11))
# patients <- patients[-1]

# set up the parameters
sigs_per_patient_matrix=matrix(0,ncol=length(unique(unlist(hdp2ref))),nrow=length(patients))
colnames(sigs_per_patient_matrix)=unique(unlist(hdp2ref))
n=0
rownames_mat=c()
sigs_per_patient=list()

# set the minimum threshold for the component contribution
exp_thresh=0.05
min_sample_num=5
patient <- patients[1]
for(patient in patients){
  
  # make a patient list
  print(patient)
  sigs_per_patient[[patient]]=list()
  
  # get the samples of interest
  samples_patient=key_table[key_table$Patient==patient, "Sample"]
  select=grepl(paste(samples_patient,collapse="|"),rownames(hdp_exposures))
  
  if(sum(select)){
    n=n+1
    # check the contribution per sample
    hdp_sigs_present=colnames(hdp_exposures)[colSums(hdp_exposures[select,]>exp_thresh)>min(min_sample_num,sum(select)-1)]
    ref_sigs_present=unique(c(unlist(hdp2ref[hdp_sigs_present])))
    sigs_per_patient[[patient]]=ref_sigs_present
    
    # if present above threshold reset to 1
    sigs_per_patient_matrix[n,ref_sigs_present]=1
    rownames_mat=c(rownames_mat, patient)
  }
}

rownames(sigs_per_patient_matrix)=rownames_mat
pdf(paste0(output.dir, "/heatmap_sigs_",exp_thresh*100,"pc_",min_sample_num,".pdf"),width=12,height=10)
heatmap(sigs_per_patient_matrix,scale='none',cexRow = 0.6,margins = c(5,13))
dev.off()
write.table(sigs_per_patient_matrix,paste0(output.dir, "/sigs_",exp_thresh*100,"pc_",min_sample_num,".txt"))

# save the signatures per patient
saveRDS(sigs_per_patient,paste0(output.dir, "/sigs_per_patient.Rdata"))

#####################################################################################
# FIG THE EXTRACTED SIGNATURES BACK INTO SIGFIT
#####################################################################################

# sigs_per_patient <- readRDS(paste0(output.dir, "/sigs_per_patient.Rdata"))

# 
sample_counts=t(mutations)
sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
full_vec = paste0(rep(c("A","C","G","T"),each=4),"[",rep(sub_vec,each=16),"]",rep(c("A","C","G","T"),times=4))
rownames(sample_counts)= full_vec

# determine the signatures which we want to fit
sigs_fitting=unique(unlist(sigs_per_patient))
sample_signatures=sample_signatures_lower=sample_signatures_upper=matrix(0,ncol=length(sigs_fitting),nrow=ncol(sample_counts))
rownames(sample_signatures)=rownames(sample_signatures_lower)=rownames(sample_signatures_upper)=colnames(sample_counts)
colnames(sample_signatures)=colnames(sample_signatures_lower)=colnames(sample_signatures_upper)=sigs_fitting

dir.create(paste0(output.dir, "/pars_objects/"))

# patients <- patients[-c(1:4)]
# iterate over each patient and fit the found signatures to the observed data
for(patient in patients){
  
  # which patient?
  print(patient)
  ref_sigs_present=sigs_per_patient[[patient]]
  samples_patient=key_table[grep(patient, key_table$Sample), "Sample"]
  
  ## Check whether there are signatures abd the analysis has not been run before
  if(!is.null(ref_sigs_present) &! file.exists(paste0(output.dir, "/pars_objects/", patient, "_2023_08_20.Rdata"))){
    
    counts=t(sample_counts[,samples_patient])
    counts=t(sample_counts[,grep(patient, colnames(sample_counts))])
    colnames(counts)=colnames(cosmic_signatures_v3.2)
    
    #ref_sigs_present=names(sort(sig_order[ref_sigs_present]))
    ref_sigs_present=unique(c(ref_sigs_present))
    
    fit=fit_signatures(counts=counts, signatures = t(final_sigs[,ref_sigs_present]),
                       iter = 20000, warmup = 10000, model="poisson", chains = 2)
    
    # retrieve the fitted exposures
    pars <- retrieve_pars(fit,par = "exposures", hpd_prob = 0.95)
    saveRDS(pars, paste0(output.dir, "/pars_objects/", patient,"_2023_08_20.Rdata"))
    
    # Plot everything, super complete:
    plot_all(fit, out_path=paste(output.dir, "/", patient, "_all_spectra",sep=""), exp_cex_names=0.9, exp_margin_bottom=14)
    
    sample_signatures[rownames(pars$mean),colnames(pars$mean)]=as.matrix(pars$mean)
    sample_signatures_lower[rownames(pars$lower),colnames(pars$lower)]=as.matrix(pars$lower)
    sample_signatures_upper[rownames(pars$upper),colnames(pars$upper)]=as.matrix(pars$upper)
  }
}

write.table(sample_signatures,paste0(output.dir, "/sample_signatures_Lung_2023_08_20.txt"))
write.table(sample_signatures_lower,paste0(output.dir, "/sample_signatures_lower_Lung_2023_08_20.txt"))
write.table(sample_signatures_upper,paste0(output.dir, "/sample_signatures_upper_Lung_2023_08_20.txt"))

# bsub -G team154-grp -o out.out -e err.err -q gpu-huge -gpu "num=1:mode=exclusive_process:mps=yes" -n 4 -R 'select[mem>=32000] rusage[mem=32000]' -M32000 "python SigprofilerExtractor.py matrix ./ NTCU_mouse_trinuc_mut_mat.txt 2 10 100"
