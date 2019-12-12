#MODIFIED DEAP LINEAR MIXED EFFECT MODEL SCRIPT
#Script originally by Wesley Thompson, PhD at University of California San Diego (2017)
#Script downloaded from ABCD Study Data Exploration and Analysis Portal (DEAP)
#Modifications made by Max M Owens, PhD at University of Vermont (Fall 2019)

#you will need to install these libraries before script will work
library(gamm4)
library(rjson)
library(stargazer)
library(knitr)
library(MuMIn)
library(R.matlab)
library(tableone)
library('parallel')
library('doParallel')
library('yaml')
library('lme4')
library('rjson')
library('parallel')
library('tidyr')
library('dplyr')
library('broom')

#read in data
data1=read.csv("/home/max/Documents/linear_mixed_model_abcd/splithalf1.csv")
data2=read.csv("/home/max/Documents/linear_mixed_model_abcd/splithalf2.csv")
data=read.csv("/home/max/Documents/linear_mixed_model_abcd/splithalfcombined.csv")

#set any factor variables as factors so that they are read into models as dummy coded rather than continuous
data[3]<-as.factor(data$sex)
data[4]<-as.factor(data$race.4level)
data[10]<-as.factor(data$married.bl)
data[11]<-as.factor(data$hisp)
data[12]<-as.factor(data$rel_family_id)
data[13]<-as.factor(data$household.income.bl)
data[20]<-as.factor(data$ksads_14_853_p)
data[21]<-as.factor(data$ksads_16_897_p)
data[22]<-as.factor(data$ksads_16_898_p)
data[23]<-as.factor(data$ksads_1_840_p)
data[24]<-as.factor(data$ksads_22_969_p)

data1[3]<-as.factor(data1$sex)
data1[4]<-as.factor(data1$race.4level)
data1[10]<-as.factor(data1$married.bl)
data1[11]<-as.factor(data1$hisp)
data1[12]<-as.factor(data1$rel_family_id)
data1[13]<-as.factor(data1$household.income.bl)
data1[20]<-as.factor(data1$ksads_14_853_p)
data1[21]<-as.factor(data1$ksads_16_897_p)
data1[22]<-as.factor(data1$ksads_16_898_p)
data1[23]<-as.factor(data1$ksads_1_840_p)
data1[24]<-as.factor(data1$ksads_22_969_p)

data2[3]<-as.factor(data2$sex)
data2[4]<-as.factor(data2$race.4level)
data2[10]<-as.factor(data2$married.bl)
data2[11]<-as.factor(data2$hisp)
data2[12]<-as.factor(data2$rel_family_id)
data2[13]<-as.factor(data2$household.income.bl)
data2[20]<-as.factor(data2$ksads_14_853_p)
data2[21]<-as.factor(data2$ksads_16_897_p)
data2[22]<-as.factor(data2$ksads_16_898_p)
data2[23]<-as.factor(data2$ksads_1_840_p)
data2[24]<-as.factor(data2$ksads_22_969_p)

#set a backup data for each iteration of the loop
user_data=data

dvs="rsfmri_cor_network.gordon_default_network.gordon_dorsalattn"
ivs=readLines("/home/max/Documents/linear_mixed_model_abcd/rsvarnames.txt")

for (d in dvs)
{
for (i in ivs)
{
#set if you want full sample or one of the split-halves below
  data=user_data
  #data=data1
  #data=data2

#Use this input list for including covariates
inputs = list(dep.var.=list(d),ind.var.=list(i),usercov.=list(),cov.fixed=list("race.4level","sex","high.educ.bl","age","household.income.bl"),smo.var=list(""),log.var=list(""),int.var=list(""),sq.var=list(""),gr.var=list(),ws.var=list(""),rand.var=list("rel_family_id","mri_info_device.serial.number"))

#Use this input list to not include covariates
#inputs = list(dep.var.=list(d),ind.var.=list(i),usercov.=list(),cov.fixed=list(""),smo.var=list(""),log.var=list(""),int.var=list(""),sq.var=list(""),gr.var=list(),ws.var=list(""),rand.var=list("rel_family_id","mri_info_device.serial.number"))


# DO NOT MAKE MODIFICATIONS AFTER THIS POINT (other than to change the output file at the end)
# all comments after this point (other than file output with sink at end) were in the script when I downloaded it from the NDA DEAP

##########################
##  user customization  ##
##########################
TEST = FALSE

extract.variables = function(a){
  rt = c()
  for (l in 1:length(a) ){ 
    if(length(a[[l]]) > 0){
      for(e in 1:length(a[[l]])){
        if(unlist(a[[l]][[e]]) != "")
          rt = c(rt, unlist(a[[l]][[e]]))
      } 
    }
  }
  rt_inster = c()
  for( item in 1:length(rt)){
    if(!is.character(rt[item])){
      next
    }
    if(length(unlist(strsplit(rt[item],"[*]"))) > 1){
      rt_inster = c(rt_inster, unlist(strsplit(rt[item],"[*]")))
    }
    else if(length(unlist(strsplit(rt[item],"[+]"))) > 1){
      rt_inster = c(rt_inster, unlist(strsplit(rt[item],"[+]")))
    }
    
    else if(length(unlist(strsplit(rt[item],"^2", fixed=TRUE))) > 1){
      rt_inster = c(rt_inster, unlist(strsplit(rt[item],"^2", fixed=TRUE)))
    }
    else if( length(regmatches(rt[item], gregexpr("(?<=\\().*?(?=\\))", rt[item], perl=T))[[1]] ) > 0 ){
      
      rt_inster = c(rt_inster, regmatches(rt[item], gregexpr("(?<=\\().*?(?=\\))", rt[item], perl=T))[[1]])
    }else{
      rt_inster = c(rt_inster, rt[item])
    }
  }
  rt = rt_inster
  return(rt);
}

if(length(unlist(inputs[['dep.var.']])) == 0 ){
  stop("Dependent variable is empty.")
}



varList.initial = extract.variables(inputs)
#if smooth by variables, need to split first
if(sum(grepl("by =",varList.initial))>=1){
  s.by.split = varList.initial[grepl("by =",varList.initial)]
  split.vars = unlist(strsplit(s.by.split, ", by = "))
  varList.initial = c(varList.initial,split.vars)
  varList.initial = varList.initial[!duplicated(varList.initial)]
}

vars.in.data = varList.initial[varList.initial %in% names(data)]
vars.keep = c("src_subject_id","rel_family_id","mri_info_deviceserialnumber",vars.in.data)
vars.keep = vars.keep[!duplicated(vars.keep)]

data = data[,vars.keep]

##################
##  functions   ##
##################

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

sep.vars = function(x){
  x = gsub(" ", "", x, fixed = TRUE)
  x = unlist(strsplit(x,"+",fixed=T))
  return(x)
}

censor =  function(x, fraction=.005){
  if(length(fraction) != 1 || fraction < 0 || fraction > 1){
    stop("bad value for 'fraction'")
  }
  lim <- quantile(x, probs=c(fraction/2, 1-fraction/2), na.rm = T)
  x[ x < lim[1] ] <- NA
  x[ x > lim[2] ] <- NA
  x
}

#########################
##  data  extraction   ##
#########################

trigger.warning = F

#exit script silently if @inputs is empyt
empty = T;
for (key in names(inputs)){

  if(length(inputs[[key]]) != 0 && inputs[[key]] != ""){
    empty = F;
  }
}

if(empty){
  options(warn=-1)
  opt <- options(show.error.messages=FALSE) 
  on.exit(options(opt)) 
  stop() 
}

#By using Rserve gamm4 is already loaded

### censor/windsorize first
wsVar = unlist(inputs[['ws.var']])
wsVar = sep.vars(wsVar)
if(length(wsVar)>0){
  for(ii in 1:length(wsVar)){
    if(class(data[,wsVar[ii]]) == "numeric"){
      data[,wsVar[ii]] = censor(data[,wsVar[ii]])
    }
  }
}

dependendVar   = unlist(inputs[['dep.var.']]);
dependendVar.name = NULL

if(length(dependendVar) == 0 ){
  stop("Dependent variable is empty.")
}


##if y is log-transformed...
if(substring(dependendVar,1,4) == "log("){
  dependendVar.name = substring(dependendVar,5,nchar(dependendVar)-1)
  if(sum(data[[dependendVar.name]] <= 0, na.rm=T) > 0){
    data[[dependendVar.name]][data[[dependendVar.name]] <= 0] = NA
    trigger.warning = T
    warning.logging.0 = paste0("1 or more log transformed variable contains values <=0. All <=0 values replaced with NA.")
  }
  data$Y.log = log(data[[dependendVar.name]])
  new.name = paste0("log.",dependendVar.name)
  names(data)[names(data) == "Y.log"] = new.name
  dependendVar = new.name
}

if( dependendVar %in% names(data)){
  if(is.factor(data[[dependendVar]])){
    if(nlevels(data[[dependendVar]]) > 2){
      stop("Categorical variables with more than 2 levels are not supported as dependent variables. \nConsider converting your categorical variable into a continuous variable.")
    }
  }
} else {
  stop(paste("Dependent variable <", dependendVar ,">does not exist in the database"));
}


if( dependendVar %in% names(data)){
  if(is.factor(data[[dependendVar]])){
    if(nlevels(data[[dependendVar]]) > 2){
      stop("Categorical variables with more than 2 levels are not supported as dependent variables. \nConsider converting your categorical variable into a continuous variable.")
    }
  }
} else {
  stop(paste("Dependent variable <", dependendVar ,">does not exist in the database"));
}


independendVar = unlist(inputs[['ind.var.']]);
usercovVar     = paste(unlist(inputs[['usercov.']]),  sep='+')

smoothVar.all = unlist(inputs[['smo.var']])
smoothVar.all = sep.vars(smoothVar.all)

logVar = unlist(inputs[['log.var']])
logVar = sep.vars(logVar)
#check if 0's in logged vars. if so, add 0.0001
strip.log = substring(logVar,5,nchar(logVar)-1)
if(length(strip.log)>0){
  if(sum(data[,strip.log]<=0 , na.rm=T) > 0){
    trigger.warning = T
    warning.logging.0 = paste0("1 or more log transformed variable contains values <=0. All <=0 values replaced with NA.")
    for(ii in 1:length(strip.log)){
      log.var_i = strip.log[ii]
      if(sum(data[,log.var_i]<=0 , na.rm=T) > 0){
        #data[[log.var_i]][data[[log.var_i]] == 0] = 0.0001
        data[[log.var_i]][data[[log.var_i]] <= 0] = NA
        
      }
    }
  }
}

interactionVar = unlist(inputs[['int.var']])

sqVar = unlist(inputs[['sq.var']])
sqVar = sep.vars(sqVar)
sqVar = substring(sqVar,1,nchar(sqVar)-2)

sqVar_SQUARED = NULL
if(length(sqVar)>0){
  for(ii in 1:length(sqVar)){
    sqVar_SQUARED[ii] = paste0(sqVar[ii],"_SQUARED")
    data[,sqVar_SQUARED[ii]] = data[,sqVar[ii]]^2
  }
}

groupVar = unlist(inputs[['gr.var']])
if(length(groupVar) > 0){
  if(is.character(groupVar) & nchar(groupVar) ==0){
    groupVar = NULL 
  } 
}

subsetVar   = unlist(inputs[['fl.var']]);
if(length(subsetVar) > 0){
  if(is.character(subsetVar) & nchar(subsetVar) ==0){
    subsetVar = NULL 
  } 
}
#Interacting grouping variable with independent variable
if(length(groupVar)>0){
  #may need to change it instead of searching for the string, it strips string first and looks for exact match
  is.smooth = grepl(independendVar, smoothVar.all)
  is.log =    independendVar %in% substring(logVar,5,nchar(logVar)-1)
  is.square = independendVar %in% sqVar
  ## If independent var is smooth
  if( sum(is.smooth) > 0 ){
    #replace s(independendVar) with s(independendVar, by = groupVar)
    smoothVar.all[is.smooth] = paste0("s(", independendVar,",by=",groupVar, ")")
  } else if(is.square){
    #else if squared, add var*groupvar and var_SQUARED*groupvar
    interactionVar = c(interactionVar, paste0(independendVar,"*",groupVar), paste0(independendVar,"_SQUARED*",groupVar) )
  } else if(is.log){    
    # else if log, interact with log(var)
    log.independent = logVar[independendVar == substring(logVar,5,nchar(logVar)-1)]
    interactionVar = c(interactionVar, paste0(log.independent,"*",groupVar) )
  } else {
    # else make normal interaction
    interactionVar = c(interactionVar, paste0(independendVar,"*",groupVar) )
  }
}  #may need to do another if with ^2 independent variables
smoothVarInt.ind = grepl(",by=", smoothVar.all)
smoothVarInt = smoothVar.all[smoothVarInt.ind]
smoothVar =    smoothVar.all[!smoothVarInt.ind]

if(0 %in% nchar(sqVar)) sqVar = character()

print(sqVar)
print(paste("length sqVar",length(sqVar)))

covfixedVar    = paste(unlist(inputs[['cov.fixed']]), sep='+')

smoothVarInt.stripped.term1 = ""
smoothVarInt.stripped.term2 = ""


if(length(smoothVarInt)>0){
  smoothVarInt.stripped.term1 = unlist(lapply( strsplit(smoothVarInt,","), function(x)x[[1]]))
  smoothVarInt.stripped.term1 = substring(smoothVarInt.stripped.term1,3,nchar(smoothVarInt.stripped.term1))
  smoothVarInt.stripped.term2 = unlist(lapply( strsplit(smoothVarInt,"by="), function(x)x[[2]]))
  smoothVarInt.stripped.term2 = substring(smoothVarInt.stripped.term2,1,nchar(smoothVarInt.stripped.term2)-1)
}

### if usercovVar's have been transformed, then they are stored in usercovVar as well as the 
### transformed vars (smoothVar, logVar, etc.), in which they will need to be removed from
### usercovVar before putting into formula
print(paste("before remove",usercovVar))
print(c( substring(smoothVar,3,nchar(smoothVar)-1),
         smoothVarInt.stripped.term1,
         smoothVarInt.stripped.term2,
         substring(logVar,5,nchar(logVar)-1),
         sqVar ))
cov.ind.remove = usercovVar %in% c( substring(smoothVar,3,nchar(smoothVar)-1),
                                    smoothVarInt.stripped.term1,
                                    smoothVarInt.stripped.term2,
                                    substring(logVar,5,nchar(logVar)-1),
                                    sqVar )
if(sum(cov.ind.remove)>0 ){
  usercovVar = usercovVar[!cov.ind.remove]
}
print(paste("after remove",usercovVar))

#remove covfixedVar when it is a smooth, smooth interaction (first variable), or log -- should usually only be age of covfixedVar
#ok to keep it in as interaction, squared, or smooth interaction term
covfixedVar.ind.remove = covfixedVar %in% c( substring(smoothVar,3,nchar(smoothVar)-1),
                                             smoothVarInt.stripped.term1,
                                             substring(logVar,5,nchar(logVar)-1))
if(sum(covfixedVar.ind.remove)>0 ){
  covfixedVar = covfixedVar[!covfixedVar.ind.remove]
}


form_arr = c(independendVar, usercovVar,covfixedVar, smoothVar.all, logVar, interactionVar, sqVar, sqVar_SQUARED);

### similarly for the independent variable...
#if independent variable is a smooth variable, log variable, or squared variable remove independendendVar from form_arr
if(independendVar %in% c(  substring(smoothVar,3,nchar(smoothVar)-1),
                           smoothVarInt.stripped.term1,
                           substring(logVar,5,nchar(logVar)-1),
                           sqVar ) ){
  form_arr = c(usercovVar,covfixedVar, smoothVar.all, logVar, interactionVar, sqVar, sqVar_SQUARED);
  #form_arr = form_arr[form_arr != independendVar]
}
form_arr = form_arr[form_arr!=""]
#########################
##  remove duplicates  ##
#########################
#take out duplicate variables
form_arr = form_arr[!duplicated(form_arr)]
formulastr = paste(dependendVar," ~ ",paste(form_arr,collapse='+'))
#get variables involve in the formula
varList = all.vars(as.formula(formulastr));
varList.independent = varList[-1]

#take out duplicated variables according to their values
if(length(varList.independent)>1){
  var.combin = combn(varList.independent,2)
  identical.values = c()
  for(ii in 1:ncol(var.combin)){
    vars_i = var.combin[,ii]
    #if two variables perfectly correlated, store both variables
    if(cor(as.numeric(data[,vars_i[1]]),as.numeric(data[,vars_i[2]]), use = "complete.obs") %in% c(1,-1)){
      # identical.values[ii] = list(vars_i)
      identical.values = c(identical.values,list(vars_i))
    }
  }
  vars.duplicate.remove = c()
  if(length(identical.values) > 0){
    for(ii in 1:length(identical.values)){
      vars_both_i = identical.values[[ii]]
      #if independendVar is in vars_both, remove other variable, otherwise doesn't matter which is removed
      if(independendVar %in% vars_both_i){
        vars.duplicate.remove[ii] = vars_both_i[vars_both_i != independendVar]
      } else{
        vars.duplicate.remove[ii] = vars_both_i[2]
      }
    }
    trigger.warning = T
    warning.duplicates = paste0("perfectly correlated variables detected, removed variable '", vars.duplicate.remove, "'")
    form_arr = form_arr[!form_arr %in% vars.duplicate.remove]
    formulastr = paste(dependendVar," ~ ",paste(form_arr,collapse='+'))
    varList = all.vars(as.formula(formulastr));
    varList.independent = varList[-1]
  }
}

print(varList)
print(formulastr)

#########################
##  data manipulation  ##
#########################
#data = data[data$eventname == "baseline_year_1_arm_1",]
print(summary(data[[independendVar]]))
#if independent variable has 5 or less unique values change it to character/factor variable
categorical.independent = FALSE

if(class(data[[independendVar]]) != "numeric"){
  categorical.independent = TRUE
}

print(summary(data[[independendVar]]))
#determine if logistic regression or not
categorical.dependent = FALSE
#data[[dependendVar]][data[[dependendVar]] == ""] = NA
if( length(table(data[[dependendVar]])) == 2 ){
  data[[dependendVar]] = as.factor(data[[dependendVar]])
  categorical.dependent = TRUE
} else{
  data[[dependendVar]] = as.numeric(as.character(data[[dependendVar]]))
}

if("household.income.bl" %in% names(data)){
  data$household.income.bl = as.character(data$household.income.bl)
  data$household.income.bl[data$household.income.bl == "[>=50K & <100K]"] = "[>=50K& <100K]"
  data$household.income.bl = as.factor(data$household.income.bl)
}

##################
##  subset data ##
##################
if(length(subsetVar)>0){
  json_data = rjson::fromJSON(file = subsetVar);
  subset = data.frame(src_subject_id=unlist(lapply(json_data[[1]]$set,function(d){ d[1] })), eventname=unlist(lapply(json_data[[1]]$set,function(d){d[2]})));
  data = merge(subset, data, all.x = T, all.y = F)
}

###################################################

data = data[complete.cases(data),]

print(dim(data))

###################################################
##  remove columns from model if 1 unique value  ##
###################################################

#less.than.2.levels = sapply(data , function(x) length(unique(x)) ) < 2
less.than.2.levels = sapply(data[,varList] , function(x) length(unique(x)) ) < 2

#trigger.warning = F
if(sum(less.than.2.levels) >0 ){
  trigger.warning = T
  #vars.to.remove = names(data)[less.than.2.levels]
  vars.to.remove = names(less.than.2.levels)[less.than.2.levels]
  
  
  #stop script if dependent variable has less than 2 unique values
  if(dependendVar %in% vars.to.remove){
    stop(paste0("'", dependendVar,"'" ," variable has <2 unique values"))
  }
  if(independendVar %in% vars.to.remove){
    stop(paste0("'", independendVar,"'" ," variable has <2 unique values"))
  }
  
  #remove vars.to.remove from form_arr, formulastr, varList, and varList.independent
  form_arr = form_arr[!form_arr %in% vars.to.remove]
  formulastr = paste(dependendVar," ~ ",paste(form_arr,collapse='+'))
  
  varList = all.vars(as.formula(formulastr));
  varList.independent = varList[-1]
  
  warning.1.unique.value = paste0("'", vars.to.remove,"'" ," variable has <2 unique values - removed from model")
}

########################################################################
##  remove independendVar (& groupVar) from formula to get delta R^2  ##
########################################################################

run.effect.size = T

if(run.effect.size){
  
  
  logVar.stripped    = substring(logVar,5,nchar(logVar)-1)
  smoothVar.stripped = substring(smoothVar,3,nchar(smoothVar)-1)
  
  if(length(interactionVar)>0){
    ##INTERACTION
    interaction.stripped.list = strsplit(interactionVar,"*",fixed=T)
  } else{
    interaction.stripped.list = character()
  }
  
  
  if(independendVar %in% smoothVar.stripped){
    ##SMOTH
    smoothVar.remove = smoothVar[smoothVar.stripped == independendVar]
    form_arr2 = form_arr[!(form_arr %in% smoothVar.remove)]
  } else if(independendVar %in% smoothVarInt.stripped.term1){
    ##SMOTH INTERACTION
    smoothVarInt.remove = smoothVarInt[smoothVarInt.stripped.term1 == independendVar]
    form_arr2 = form_arr[!(form_arr %in% smoothVarInt.remove)]
  } else if(independendVar %in% logVar.stripped){
    ##LOG
    logVar.remove = logVar[logVar.stripped == independendVar]
    form_arr2 = form_arr[!(form_arr %in% logVar.remove)]
    if(length(groupVar)>0){
      logVarinteraction.remove = c(paste0(logVar.remove,"*",groupVar) , groupVar)
      form_arr2 = form_arr2[!(form_arr2 %in% logVarinteraction.remove)]
    }
  } else if(independendVar %in% sqVar){
    ##SQUARED
    sqVar.remove = c(independendVar, paste0(independendVar,"_SQUARED"))
    form_arr2 = form_arr[!(form_arr %in% sqVar.remove)]
    if(length(groupVar)>0){
      sqVarinteraction.remove = c(paste0(sqVar.remove,"*",groupVar), groupVar)
      form_arr2 = form_arr2[!(form_arr2 %in% sqVarinteraction.remove)]
    }
  } else if(independendVar %in% unlist(interaction.stripped.list)){
    ##INTERACTION
    independ.location = c()
    for(ii in 1:length(interaction.stripped.list)){
      int.stripped.list_i = interaction.stripped.list[[ii]]
      independ.location[ii] = sum(independendVar == int.stripped.list_i)
    }
    interaction.remove = interactionVar[independ.location == 1]
    interactionVars.remove = unlist(strsplit(interaction.remove,"*",fixed=T))
    all.interactionVars.remove = c(interaction.remove,interactionVars.remove)
    form_arr2 = form_arr[!(form_arr %in% all.interactionVars.remove)]
  } else {
    #NO INDEPENDENT VAR TRANSFORMATION
    form_arr2 = form_arr[!(form_arr %in% independendVar)]
  }
  
  if(length(form_arr2) > 0){
    formulastr2 = paste(dependendVar," ~ ",paste(form_arr2,collapse='+'))
  } else{
    formulastr2 = paste(dependendVar," ~ 1")
  }
  
  
}

######################
##  random effects  ##
######################
include.random.scanner = "mri_info_deviceserialnumber" %in% unlist(inputs$rand.var)
include.random.site    = "abcd_site" %in% unlist(inputs$rand.var)
include.random.subject = "src_subject_id" %in% unlist(inputs$rand.var)

#IF less than 2 unique "eventname" & "include.random.subject" toggled on, 
#   need to toggle off random subject effect and produce warning (or error)
if(length(unique(data$eventname)) <2 & include.random.subject){
  include.random.subject = FALSE
  trigger.warning = T
  warning.1.unique.eventname = paste0("Your model data contains a single data point per subject (event name “baseline_year_1_arm_1”) only. To prevent problems with model convergence suggest to remove SUBJECT from the list of random effects. It is not required with this particular cross-sectional analysis.")
}
#IF 2+ unique "eventname" & "include.random.subject" toggled off, 
#   need to toggle on random subject effect and produce warning (or error)
if(length(unique(data$eventname)) > 1 & !include.random.subject){
  # include.random.subject = TRUE
  trigger.warning = T
  warning.2.unique.eventname = paste0("Your model data contains multiple data points per subject (event name “baseline_year_1_arm_1”). Suggest to include SUBJECT from the list of random effects.")
}

formula.random.str = c()
if(include.random.site & !include.random.subject){
  formula.random.str = "(1|abcd_site/rel_family_id)"
} else if(include.random.scanner & !include.random.subject){
  formula.random.str = "(1|mri_info_deviceserialnumber/rel_family_id)"
} else{
  if(include.random.site) formula.random.str[1] = "(1|abcd_site)"
  if(include.random.scanner) formula.random.str[1] = "(1|mri_info_deviceserialnumber)"
  formula.random.str[2] = "(1|rel_family_id)"
  if(include.random.subject) formula.random.str[3] = "(1|src_subject_id)"
}
formula.random.str = formula.random.str[!is.na(formula.random.str)]
formula.random.str = paste0("~",paste(formula.random.str, collapse = "+"))
formula.random = formula.random.str

#################
##  run model  ##
#################

if(categorical.dependent){ 
  formula_list = c(formulastr, formulastr2, paste0(dependendVar, " ~1"));
  cl <- makeCluster(3)
  registerDoParallel(3)
  model_list = foreach(ii= 1:length(formula_list)) %dopar% {
    gamm4(as.formula(formula_list[ii]), random = as.formula(formula.random), family = binomial, data = data)
  }
  print(model_list)
  model = model_list[[1]]
  model2 = model_list[[2]]
  model3 = model_list[[3]]
  stopCluster(cl)
  
  
} else{ #linear regression
  formula_list = c(formulastr, formulastr2, paste0(dependendVar, " ~1"));
  cl <- makeCluster(3)
  registerDoParallel(3)
  model_list = foreach(ii= 1:length(formula_list),.packages='gamm4') %dopar% {
    gamm4(as.formula(formula_list[ii]), random = as.formula(formula.random), data = data)
  }
  print(model_list)
  model = model_list[[1]]
  model2 = model_list[[2]]
  model3 = model_list[[3]]
  stopCluster(cl)
}

print(summary(model))
print(summary(model2$gam))

########################
##  other statistics  ##
########################
#N subjects, aic, bic, & r^2
n = summary(model$gam)$n
aic = round(AIC(model$mer),2)
bic = round(BIC(model$mer),2)
r2       = round(as.numeric(r.squaredLR(model$mer,model3$mer)),5)
r2_delta = round(as.numeric(r.squaredLR(model$mer,model2$mer)),5)
aic2 = round(AIC(model2$mer),2)
bic2 = round(BIC(model2$mer),2)


#Print results of model to text file
#CHANGE NAME OF OUTPUT FILE HERE

sink("/home/max/Documents/linear_mixed_model_abcd/rs_analysis_output.txt",append=TRUE)
print(summary(model$gam))
#Turn on line below to print model without IV as well
#print(summary(model2$gam))
sink()

  }
}
