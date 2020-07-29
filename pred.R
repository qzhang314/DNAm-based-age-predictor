##################################################################################################
###                                                                                              #      
###        This Script is used to do age prediction based on DNA methylation data(450K)          #
###        Coefficients of the predictor are based on 13,566 training samples                    #
###        Two age predictors based on different methods(Elastic Net and BLUP) can be used       #
###        Qian Zhang 27-03-2018, Email: q.zhang@uq.edu.au                                       #
##################################################################################################




############# for each probe, change to missing value to the mean value across all individuals #############
addna<-function(methy){
        methy[is.na(methy)]<-mean(methy,na.rm=T)
        return(methy)
}


############# 1. get the parameters ##################
args<-commandArgs(TRUE)
args<- as.character(args)

infile<-c()                     ########## input file name
outfile<-c()                    ########## output file name
agefile<-c()                    ########## file with individual ID and age

for(i in 1:(length(args)-1)){
        if(args[i]=="-i"){
                infile<- args[i+1]
        }else if(args[i]=="-o"){
                outfile<- args[i+1]
        }else if(args[i]=="-a"){
                agefile<- args[i+1]
        }
}

if(is.null(infile) | is.null(outfile) | is.null(agefile)){
        print("Error: please check if the filenames or flags are right")
        quit()
}

############# 2. data loading and QC ##################
print("1. Data loading and QC")

print("1.1 Reading the data")
readRDS(infile)-> data        ########## IND * Probe, each row represents one individual, it should be "RAW BETA" DNA methylation value

if(nrow(data) > ncol(data)){
        print("I guess you are using Probe in the row, data will be transformed!!!")
        data<-t(data)
}

print("1.2 Replacing missing values with mean value")
dataNona<-apply(data,2,function(x) addna(x))   ###############  replace the NA with mean value for each probe 

dataNona<- dataNona[,apply(dataNona,2,function(x) sum(is.na(x)))!=nrow(dataNona)] ############ remove the probe when it has NA across all individuals
print(paste0(ncol(data) - ncol(dataNona)," probe(s) is(are) removed since it has (they have) NA across all individuals"))


print("1.3 Standardizing")
dataNona.norm<- apply(dataNona,1,scale)        ############### standardize the DNA methylation within each individual, remove the mean and divided by the SD of each individual     Probe * IND
rownames(dataNona.norm)<-colnames(dataNona)


############# 3. get the coefficients of each probe from Elastic Net/BLUP method, !!!!WE HAVE TWO PREDICTORS!!!#############

print("2. Loading predictors")
read.table("en.coef",stringsAsFactor=F,header=T)->encoef
read.table("blup.coef",stringsAsFactor=F,header=T)->blupcoef

en_int<-encoef[1,2]
blup_int<-blupcoef[1,2]

encoef<-encoef[-1,]
blupcoef<-blupcoef[-1,]

rownames(encoef)<-encoef$probe
rownames(blupcoef)<-blupcoef$probe

############# 4. get common probes between predictors and data ##############
print("3. Checking misssing probes")

encomm<- intersect(rownames(encoef),rownames(dataNona.norm))
blupcomm<- intersect(rownames(blupcoef),rownames(dataNona.norm))

endiff<- nrow(encoef) - length(encomm)
blupdiff<- nrow(blupcoef) - length(blupcomm)

print(paste0(endiff," probe(s) in Elastic Net predictor is(are) not in the data"))
print(paste0(blupdiff," probe(s) in BLUP predictor is(are) not in the data"))
print("BLUP can perform better if the number of missing probes is too large!")

############# 5. extract the common probes and do age prediction ###############
print("4. Predicting")

encoef<-encoef[encomm,]
blupcoef<-blupcoef[blupcomm,]
encoef$coef%*%dataNona.norm[encomm,]+en_int->enpred
blupcoef$coef%*%dataNona.norm[blupcomm,]+blup_int->blupred


############# 6. Save the predicted result ###########
read.table(agefile,header=T,stringsAsFactor=F)->age.raw
enpred<-enpred[,age.raw$ID]
blupred<-blupred[,age.raw$ID]

age.raw$enpred<-as.double(enpred)
age.raw$blupred<-as.double(blupred)


write.table(age.raw,file=outfile,row.names=F,quote=F)

################ DONE ########################################

print("Completed!!!")
quit()

		
############ This is the script for plotting figure and calculating prediction accuracy ##########
read.table(file=outfile,stringsAsFactor=F,header=T) -> age
library(ggplot2)
library(reshape2)
colnames(age)<-c("ID","age.raw","Elastic Net","BLUP")
corv<-c(round(cor(age[,2],age[,3]),2),round(cor(age[,2],age[,4]),2))
rmse<-c(round(sqrt(mean((age[,2]-age[,3])^2)),2),round(sqrt(mean((age[,2]-age[,4])^2)),2))
age<-melt(age,measure.vars=c("Elastic Net","BLUP"))
ggplot(data=age,aes(x=age.raw,y=value))+facet_wrap(~variable)+geom_abline(intercept=0,slope=1)+geom_point()+xlab("Chronological Age")+ylab("Predicted Age")+theme(axis.title=element_text(size=12,face="bold"),legend.text=element_text(size=10,face="bold"),legend.title=element_text(face="bold"),axis.text=element_text(size=10,face="bold"),strip.text.x=element_text(face="bold"))+annotate("text", label = paste0("Corr = ",corv), size = 4, x = min(age$age.raw)+5, y = max(age$age.raw)-5,hjust=0)+annotate("text", label = paste0("RMSE = ",rmse), size = 4, x = min(age$age.raw)+5, y = max(age$age.raw)-6.5, hjust=0)
ggsave("prediction.png", width = 12, height = 6, dpi = 400)


