#####build mixed model function####
mixed_model_nocov=function(x,y,data){
  stat_holder <- data.frame()
  stat_names <- c('B','SE','t','p','R2')
  for (k in stat_names) stat_holder[k] <- as.double()
  
  nullmod <- formula(paste(y, "~", 1 ))
  form_cov_only <- formula(paste(y, "~", x ))
  model <- gamm4(form_cov_only, data=data, 
                 random =~ (1|mri_info_device.serial.number/rel_family_id) )
  nlmodel <- gamm4(nullmod, data=data, 
                   random =~ (1|mri_info_device.serial.number/rel_family_id) )
  r2_delta <- round(as.numeric(r.squaredLR(model$mer, nlmodel$mer)), 5)
  sg <- summary(model$gam)
  
  for (statnum in 1:4) {stat_holder[1, statnum] <- sg$p.table[2, statnum]}
  stat_holder[1, 5] <- r2_delta
  return(stat_holder)
}  