rm(list=ls())
library(TwoSampleMR)
library(MRPRESSO)
library(tidyverse)
library(data.table)
library(ggplot2)
library(openxlsx)


setwd("your work place")
gut<-fread("Data_MBG.allHits.p1e4.txt",data.table = F, header = T)
gut_p1e5<-gut[which(gut$P.weightedSumZ < 1e-5),]
gut_all<-gut_p1e5
bac.list<-read.xlsx("GutList.xlsx",sheet = "AllBac",colNames = FALSE)

###############################MR for each exposure
for(i in 1:dim(bac.list)[1]){
####exposure for gut microbiota
bacname<-bac.list[i,]
gut_sin<-gut_all[which(gut_all$bac==bacname),]
exp_dat<-format_data(dat=gut_sin,type='exposure',phenotype_col = "bac",snp_col = "rsID",beta_col = "beta",se_col = "SE",effect_allele_col ="eff.allele",other_allele_col = "ref.allele",pval_col = "P.weightedSumZ",samplesize_col="N",chr_col="chr")
dim(exp_dat)

####outcome for BMI
out_dat<-extract_outcome_data(snps=exp_dat$SNP,outcomes='ieu-b-40',proxies=FALSE,maf_threshold=0.01)
dim(out_dat)

####Data harmonization
mydata <- harmonise_data(exposure_dat=exp_dat,outcome_dat=out_dat,action=2)
mydata <- mydata[which(mydata$mr_keep==TRUE),]
dim(mydata)

####clumping process
exp_dat_harm<-filter(exp_dat,exp_dat$SNP %in% mydata$SNP)
dim(exp_dat_harm)
exp_dat_clump <-clump_data(exp_dat_harm,clump_r2=0.01,clump_kb=500,pop = "EUR")
dim(exp_dat_clump)
mydata_clump<-filter(mydata,mydata$SNP %in% exp_dat_clump$SNP)

####removing outcome-related SNPs
mydata_filteroutcome<-mydata_clump[which(mydata_clump$pval.outcome>=0.05),]

####MR-PRESSO, removing the outliers
####The SNPs were sorted in ascending order according to the P-values of MR-PRESSO outlier test, and the remaining SNPs were eliminated one by one until there was no pleiotropy.
####When the number of SNPs is less than or equal to 3, it is indicated that there are not enough intrumental variables for MR-PRESSO.
mydata_outcomen<-dim(mydata_filteroutcome)[1]
if(mydata_outcomen<=3){
	presso_pval<-"NA"
	mydata_presso<-mydata_filteroutcome
	}else{
		presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = mydata_filteroutcome, NbDistribution = 1000,  SignifThreshold = 0.05)
		presso_pval<-presso$`MR-PRESSO results`$`Global Test`$Pvalue
		presso_snp<-presso$`MR-PRESSO results`$`Outlier Test`
		if(presso_pval>=0.05){
			mydata_presso<-mydata_filteroutcome
		}else{
			####removing the outliers
			out_order<-order(presso_snp$Pvalue)
			out_order
			for(ii in 1:length(out_order)){
				snp_ii<-out_order[1:ii]
				mydata_presso<-mydata_filteroutcome[-snp_ii,]
				mydata_presn<-dim(mydata_presso)[1]
				if(mydata_presn<=3){
					presso_pval<-"NA"
					}else{
					presso2 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = mydata_presso, NbDistribution = 1000,  SignifThreshold = 0.05)
					presso_pval<-presso2$`MR-PRESSO results`$`Global Test`$Pvalue
					mydata_presso<-mydata_presso
					}
				print(ii);
				if(presso_pval>=0.05){
					break;
				}
			}
		}
}
dim(mydata_presso)


####Steiger filtering
steiger_snpi<-steiger_filtering(mydata_presso)
steiger_filtering_snpi<-steiger_snpi[which(steiger_snpi$steiger_dir==TRUE),]
steiger_filtering_snpid<-steiger_filtering_snpi$SNP
mydata_presso<-mydata_presso[which(mydata_presso$SNP %in% steiger_filtering_snpid),]

####MR Steiger directionality test
out <- directionality_test(mydata_presso)

####F statistic
mydata_presn<-dim(mydata_presso)[1]
k<-mydata_presn
N<-18340
r2.exposure<-out$snp_r2.exposure
F.exposure<-(r2.exposure/(1-r2.exposure))*((N-k-1)/k)
out$F.exposure<-F.exposure

####MR analysis
res<-mr(mydata_presso, method_list=c("mr_egger_regression","mr_ivw","mr_ivw_mre","mr_ivw_fe","mr_wald_ratio","mr_weighted_median","mr_two_sample_ml","mr_simple_mode","mr_weighted_mode"))
res_or<-generate_odds_ratios(res)

####Sensitivity analysis
####heterogeneity
het<-mr_heterogeneity(mydata_presso, method_list=c("mr_egger_regression", "mr_ivw"))
####pleiotropy
plt<-mr_pleiotropy_test(mydata_presso)
####leave one out
res_loo<-mr_leaveoneout(mydata_presso)

####plots
####scatter plot
p1<-mr_scatter_plot(res,mydata_presso)
####forest plot
res_single<-mr_singlesnp(mydata_presso)
p2<-mr_forest_plot(res_single)
####funnel plot
p3<-mr_funnel_plot(res_single)
####leave one out plot
p4<-mr_leaveoneout_plot(res_loo)

####Save results
setwd("your work place")
dir.create(bacname)
filepath<-paste("your work place",bacname,sep="")
setwd(filepath)
ggsave(p1[[1]], file="mr_scatter_plot.pdf")
ggsave(p2[[1]], file="mr_forest_plot.pdf")
ggsave(p3[[1]], file="mr_funnel_plot.pdf")
ggsave(p4[[1]], file="mr_leaveoneout_plot.pdf")
write.table(mydata_presso,"mydata_presso.txt",sep="\t",quote = FALSE,row.names = FALSE)#IVs
write.table(res,"res.txt",sep="\t",quote = FALSE,row.names = FALSE)#res
write.table(res_or,"res_or.txt",sep="\t",quote = FALSE,row.names = FALSE)#res_or
write.table(het,"het.txt",sep="\t",quote = FALSE,row.names = FALSE)#heterogeneity
write.table(plt,"plt.txt",sep="\t",quote = FALSE,row.names = FALSE)#pleiotropy
write.table(steiger_snpi,"steiger_snpi.txt",sep="\t",quote = FALSE,row.names = FALSE)#steiger_filtering
write.table(out,"out.txt",sep="\t",quote = FALSE,row.names = FALSE)#MR Steiger
}


###############################Other outcomes
####outcome for Obesity Class I
out_dat<-extract_outcome_data(snps=exp_dat$SNP,outcomes='ieu-a-90',proxies=FALSE,maf_threshold=0.01)
dim(out_dat)

####outcome for Obesity Class II
out_dat<-extract_outcome_data(snps=exp_dat$SNP,outcomes='ieu-a-91',proxies=FALSE,maf_threshold=0.01)
dim(out_dat)

####outcome for Obesity Class III
out_dat<-extract_outcome_data(snps=exp_dat$SNP,outcomes='ieu-a-92',proxies=FALSE,maf_threshold=0.01)
dim(out_dat)

####outcome for WHR
setwd("your work place")
whr<-fread("whr.giant-ukbb.meta-analysis.combined.23May2018.txt",data.table = F,header=T)
whr$phenotype <- 'whr'
whr$SNP<-gsub("\\:.*", "", whr$SNP)
out_dat <- format_data(dat=whr,type = "outcome",snps = exp_dat$SNP,phenotype_col = "phenotype",snp_col = "SNP",beta_col = "BETA",se_col = "SE",effect_allele_col = "Tested_Allele",other_allele_col = "Other_Allele",pval_col = "P",chr_col = "CHR",pos_col = "POS",eaf_col="Freq_Tested_Allele",samplesize_col="N")
out_dat<-out_dat[which(out_dat$eaf.outcome>0.01),]
dim(out_dat)

####outcome for WHRadjBMI
setwd("your work place")
whradjbmi<-fread("whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt",data.table = F,header=T)
whradjbmi$phenotype <- 'whradjbmi'
whradjbmi$SNP<-gsub("\\:.*", "", whradjbmi$SNP)
out_dat <- format_data(dat=whradjbmi,type = "outcome",snps = exp_dat$SNP,phenotype_col = "phenotype",snp_col = "SNP",beta_col = "BETA",se_col = "SE",effect_allele_col = "Tested_Allele",other_allele_col = "Other_Allele",pval_col = "P",chr_col = "CHR",pos_col = "POS",eaf_col="Freq_Tested_Allele",samplesize_col="N")
out_dat<-out_dat[which(out_dat$eaf.outcome>0.01),]
dim(out_dat)
