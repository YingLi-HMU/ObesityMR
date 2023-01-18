rm(list=ls())
library(TwoSampleMR)
library(MRPRESSO)
library(tidyverse)
library(data.table)
library(ggplot2)
library(openxlsx)


####exposure for BMI
exp_dat <-extract_instruments(outcomes='ieu-b-40',p1=5e-8,clump=TRUE,p2 = 5e-08,r2 = 0.001,kb = 10000,access_token = NULL)
dim(exp_dat)
N<-681275
###############################Other Exposures
####exposure for Obesity Class I
exp_dat <-extract_instruments(outcomes='ieu-a-90',p1=5e-8,clump=TRUE,p2 = 5e-08,r2 = 0.001,kb = 10000,access_token = NULL)
dim(exp_dat)
N<-98697
####exposure for Obesity Class II
exp_dat <-extract_instruments(outcomes='ieu-a-91',p1=5e-8,clump=TRUE,p2 = 5e-08,r2 = 0.001,kb = 10000,access_token = NULL)
dim(exp_dat)
N<-72546
####exposure for Obesity Class III
exp_dat <-extract_instruments(outcomes='ieu-a-92',p1=5e-8,clump=TRUE,p2 = 5e-08,r2 = 0.001,kb = 10000,access_token = NULL)
dim(exp_dat)
N<-50364
####exposure for WHR
setwd("your work place")
whr<-fread("whr.giant-ukbb.meta-analysis.combined.23May2018.txt",data.table = F,header=T)
whr$phenotype <- 'whr'
whr$SNP<-gsub("\\:.*", "", whr$SNP)
#The number of SNPS in the GWAS summary datasets of WHR and WHRadjBMI was much larger than that of other obesity traits. Therefore, the threshold is set at 5e-9.
whr<-whr[which(whr$P < 5e-9),]
exp_dat<-format_data(dat=whr,type='exposure',phenotype_col = "phenotype",snp_col = "SNP",beta_col = "BETA",se_col = "SE",effect_allele_col = "Tested_Allele",other_allele_col = "Other_Allele",pval_col = "P",chr_col = "CHR",pos_col = "POS",eaf_col="Freq_Tested_Allele",samplesize_col="N")
exp_dat_clump<-clump_data(exp_dat,clump_r2=0.001,clump_kb=10000)
dim(exp_dat_clump)
N<-697734
####exposure for WHRadjBMI
setwd("your work place")
whradjbmi<-fread("whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt",data.table = F,header=T)
whradjbmi$phenotype <- 'whradjbmi'
whradjbmi$SNP<-gsub("\\:.*", "", whradjbmi$SNP)
#The number of SNPS in the GWAS summary datasets of WHR and WHRadjBMI was much larger than that of other obesity traits. Therefore, the threshold is set at 5e-9.
whradjbmi<-whradjbmi[which(whradjbmi$P < 5e-9),]
exp_dat<-format_data(dat=whr,type='exposure',phenotype_col = "phenotype",snp_col = "SNP",beta_col = "BETA",se_col = "SE",effect_allele_col = "Tested_Allele",other_allele_col = "Other_Allele",pval_col = "P",chr_col = "CHR",pos_col = "POS",eaf_col="Freq_Tested_Allele",samplesize_col="N")
exp_dat_clump<-clump_data(exp_dat,clump_r2=0.001,clump_kb=10000)
dim(exp_dat_clump)
N<-694649

####outcome for phylum.Actinobacteria.id.400
setwd("your work place")
bac<-fread("phylum.Actinobacteria.id.400.summary.txt",data.table = F,header=T)
bacname<-"phylum.Actinobacteria.id.400"
out_dat <- format_data(dat=bac,type = "outcome",snps = exp_dat$SNP,phenotype_col = "bac",snp_col = "rsID",beta_col = "beta",se_col = "SE",effect_allele_col = "eff.allele",other_allele_col = "ref.allele",pval_col = "P.weightedSumZ",chr_col = "chr",pos_col = "bp",samplesize_col="N")
dim(out_dat)
###############################Other Outcomes
####outcome for order.Bifidobacteriales.id.432
setwd("your work place")
bac<-fread("order.Bifidobacteriales.id.432.summary.txt",data.table = F,header=T)
bacname<-"order.Bifidobacteriales.id.432"
####outcome for family.Bifidobacteriaceae.id.433
setwd("your work place")
bac<-fread("family.Bifidobacteriaceae.id.433.summary.txt",data.table = F,header=T)
bacname<-"family.Bifidobacteriaceae.id.433"
####outcome for genus.LachnospiraceaeUCG008
setwd("your work place")
bac<-fread("genus.LachnospiraceaeUCG008.id.11328.summary.txt",data.table = F,header=T)
bacname<-"genus.LachnospiraceaeUCG008.id.11328"
####outcome for genus..Eubacteriumnodatumgroup.id.11297
setwd("your work place")
bac<-fread("genus..Eubacteriumnodatumgroup.id.11297.summary.txt",data.table = F,header=T)
bacname<-"genus..Eubacteriumnodatumgroup.id.11297"

####removing outcome-related SNPs
out_dat<-out_dat[-which(out_dat$pval.outcome<0.05),]

####Data harmonization
mydata <- harmonise_data(exposure_dat=exp_dat,outcome_dat=out_dat,action=2)
mydata_clump <- mydata[which(mydata$mr_keep==TRUE),]
mydata_filteroutcome<-mydata_clump

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
