#' @title Normalize MS TMT data.
#' @description Normalize MS TMT data which contain no zero rows.
#' @details Input dataframe from read_maxquant_prot() function, then return a data set after normalization with SL, TMM and IRS.
#' @param prot_dat A DataFrame from read_maxquant_prot() with only numeric value.
#' @param condition_num A integer of condition number.
#' @param repeat_num A integer of repeat number.
#' @return A DataFrame.
#' @return Plots.
#' @export
#' @import limma
#' @import edgeR
#' @import tidyverse
#' @import psych
#' @import RColorBrewer
#' @examples prot_norm_example <- data_norm(prot_dat_example[1:15], 5, 3)

data_norm <- function(prot_dat,condition_num,repeat_num,reference=FALSE,ref_col = c(1,7,13)){

  if (reference == TRUE){
    data_graphic_draw_norm(prot_dat[,-ref_col],condition_num,repeat_num,'Before Normalization')

    prot_sl <- data_norm_sl(prot_dat)
    data_graphic_draw_norm(prot_sl[,-ref_col],condition_num,repeat_num,'SL Normalization')

    prot_sl_tmm <- data_norm_tmm(prot_sl)
    data_graphic_draw_norm(prot_sl_tmm[,-ref_col],condition_num,repeat_num,'SL/TMM Normalization')

    prot_sl_tmm_ers <- data_norm_ers(prot_sl_tmm,condition_num,repeat_num,ref_col)
    data_graphic_draw_norm(prot_sl_tmm_ers[,-ref_col],condition_num,repeat_num,'SL/TMM/ERS Normalization')
    prot_norm <- prot_sl_tmm_ers

    data_pearson_draw(prot_dat, condition_num + 1, repeat_num,'Before Normalization')
    data_pearson_draw(prot_norm, condition_num + 1, repeat_num,'SL/TMM/ERS Normalization')

    prot_norm <- prot_norm[,-ref_col]
  }
  else{
    data_graphic_draw_norm(prot_dat,condition_num,repeat_num,'Before Normalization')

    prot_sl <- data_norm_sl(prot_dat)
    data_graphic_draw_norm(prot_sl,condition_num,repeat_num,'SL Normalization')

    prot_sl_tmm <- data_norm_tmm(prot_sl)
    data_graphic_draw_norm(prot_sl_tmm,condition_num,repeat_num,'SL/TMM Normalization')

    prot_sl_tmm_irs <- data_norm_irs(prot_sl_tmm,condition_num,repeat_num)
    data_graphic_draw_norm(prot_sl_tmm_irs,condition_num,repeat_num,'SL/TMM/IRS Normalization')
    prot_norm <- prot_sl_tmm_irs

    data_pearson_draw(prot_dat, condition_num, repeat_num,'Before Normalization')
    data_pearson_draw(prot_norm, condition_num, repeat_num,'SL/TMM/IRS Normalization')
  }

  return(prot_norm)
}

#' @title Normalize MS TMT data with Sample Loading.
#' @description Normalize MS TMT data with sample loading.
#' @details Input dataframe from read_maxquant_prot() function, then return a data set after normalization with SL.
#' @param prot_dat A DataFrame from read_maxquant_prot() with only numeric value.
#' @return A DataFrame.
#' @export
data_norm_sl <- function(prot_dat){
  target <- mean(colSums(prot_dat))
  norm_facs <- target / colSums(prot_dat)
  prot_dat_sl <- sweep(prot_dat,2,norm_facs,FUN = '*')
  return(prot_dat_sl)
}

#' @title Normalize MS TMT data with TMM.
#' @description Normalize MS TMT data with Trimmed Mean of -values.
#' @details Input dataframe from read_maxquant_prot() function, then return a data set after normalization with TMM.
#' @param prot_dat A DataFrame from read_maxquant_prot() with only numeric value.
#' @return A DataFrame.
#' @import edgeR
#' @export
data_norm_tmm <- function(prot_dat){
  target <- calcNormFactors(prot_dat)
  prot_dat_tmm <- sweep(prot_dat,2,target,FUN='/')
  return(prot_dat_tmm)
}

#' @title Normalize MS TMT data with IRS.
#' @description Normalize MS TMT data with Internal Reference Signal.
#' @details Input dataframe from read_maxquant_prot() function, then return a data set after normalization with IRS.
#' @param prot_dat A DataFrame from read_maxquant_prot() with only numeric value.
#' @return A DataFrame.
#' @export
data_norm_irs <- function(prot_dat,condition_num,repeat_num){
  i <- 1
  irs <- cbind(as.numeric(rownames(prot_dat)),rowSums(prot_dat[,1:condition_num]))
  while (i < repeat_num){
    start_site <- i * condition_num + 1
    stop_site <- i * condition_num + condition_num
    irs <- cbind(irs, rowSums(prot_dat[,start_site:stop_site]))
    i <- i + 1
  }
  irs <- irs[,-1]
  irs$average <- apply(irs,1,function(x){exp(mean(log((x))))})

  i <- 1
  prot_dat_irs <- cbind(rownames(prot_dat),(prot_dat[,1:condition_num] * irs$average / rowSums(prot_dat[,1:condition_num])))
  while (i < repeat_num){
    start_site <- i * condition_num + 1
    stop_site <- i * condition_num + condition_num
    prot_dat_irs <- cbind(prot_dat_irs,(prot_dat[,start_site:stop_site] * irs$average / rowSums(prot_dat[,start_site:stop_site])))
    i <- i + 1
  }

  rownames(prot_dat_irs) <- prot_dat_irs[,1]
  prot_dat_irs <- prot_dat_irs[,-1]

  return(prot_dat_irs)
}

#' @title Normalize MS TMT data with ERS.
#' @description Normalize MS TMT data with External Reference Signal.
#' @details Input dataframe from read_maxquant_prot() function, then return a data set after normalization with ERS.
#' @param prot_dat A DataFrame from read_maxquant_prot() with only numeric value.
#' @return A DataFrame.
#' @export
data_norm_ers <- function(prot_dat,condition_num,repeat_num,ref_col = c(1,7,13)){
  prot_ref <- prot_dat[,ref_col]
  temp <- prot_ref
  temp$rowsum <- rowSums(temp)
  prot_ref <- temp$rowsum / prot_ref

  i <- 1
  while(i <= ncol(prot_ref)){
    start_site <- (i - 1) * (condition_num + 1) + 1
    stop_site <- i * (condition_num + 1)
    prot_dat[,start_site : stop_site] <- prot_dat[,start_site : stop_site] * prot_ref[,i]
    i <- i + 1
  }

  return(prot_dat)
}

#' @title Draw Dataset Summary.
#' @description Draw dataset summary with density distribution, pca and cv.
#' @details Input dataframe, then return a dataset summary with density distribution, pca and cv.
#' @param prot_dat A DataFrame with only numeric value.
#' @param condition_num A integer of condition number.
#' @param repeat_num A integer of repeat number.
#' @param state A title
#' @export
#' @import limma
#' @import edgeR
#' @import tidyverse
#' @import psych
#' @import RColorBrewer
#' @examples prot_norm_example <- data_norm(prot_dat_example[1:15], 5, 3)
data_graphic_draw_norm <- function(prot_dat, condition_num, repeat_num, state){
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  par(mfrow = c(2, 2))
  boxplot(log2(prot_dat),col = rep(col_vector[9:(repeat_num+8)],each=condition_num), main = state)
  plotDensities(log2(prot_dat),col = rep(col_vector[9:(repeat_num+8)],condition_num), main = state)
  plotMDS(log2(prot_dat), col = col_vector[9:(condition_num+8)], main = state)
  data_calCV_draw(prot_dat,condition_num,repeat_num, state)
  par(mfrow = c(1, 1)) # reset to default
}

data_rearrange <- function(prot_dat,condition_num,repeat_num){
  i <- 1
  order_seq <- seq(1,(condition_num * repeat_num - condition_num + 1),condition_num)
  while(i < condition_num){
    i <- i + 1
    temp <- seq(i,(condition_num * repeat_num - condition_num + i),condition_num)
    order_seq <- c(order_seq,temp)
  }
  prot_dat_reorder <- prot_dat[,order_seq]
  return(prot_dat_reorder)
}

data_calCV_draw <- function(prot_dat,condition_num,repeat_num,state){
  prot_dat_reorder <- data_rearrange(prot_dat,condition_num,repeat_num)

  i <- 1
  mean_list <- apply(prot_dat_reorder[,1 : repeat_num], 1, mean)
  sd_list <- apply(prot_dat_reorder[,1 : repeat_num], 1, sd)
  cv_list <- sd_list / mean_list * 100
  prot_dat_cv <- cbind(as.numeric(rownames(prot_dat)),cv_list)
  while (i < condition_num) {
    mean_list <- apply(prot_dat_reorder[,(i * repeat_num + 1) : (i * repeat_num + repeat_num)], 1, mean)
    sd_list <- apply(prot_dat_reorder[,(i * repeat_num + 1) : (i * repeat_num + repeat_num)], 1, sd)
    cv_list <- sd_list / mean_list * 100
    prot_dat_cv <- cbind(prot_dat_cv,cv_list)
    i <- i + 1
  }
  prot_dat_cv <- prot_dat_cv[,-1]
  cond_seq <- seq(1,condition_num,1)
  cond_seq_title <- paste('cond',cond_seq,sep = '')
  colnames(prot_dat_cv) <- cond_seq_title

  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  boxplot(prot_dat_cv, col = col_vector[9:(condition_num+8)], main = state)
}

#' @title Draw Pearson Plot.
#' @description Draw Pearson Plot.
#' @details Input dataframe, then return a set of Pearson plots.
#' @param prot_dat A DataFrame from with only numeric value.
#' @param condition_num A integer of condition number.
#' @param repeat_num A integer of repeat number.
#' @param state A title.
#' @export
#' @import limma
#' @import edgeR
#' @import tidyverse
#' @import psych
data_pearson_draw <- function(prot_dat, condition_num, repeat_num,state=''){
  prot_dat_reorder <- data_rearrange(prot_dat,condition_num,repeat_num)
  i <- 1
  pairs.panels(log2(prot_dat_reorder[,1:repeat_num]), lm = TRUE,main=state)
  while (i < condition_num) {
    start_site <- i * repeat_num + 1
    stop_site <- i * repeat_num + repeat_num
    pairs.panels(log2(prot_dat_reorder[,start_site:stop_site]), lm = TRUE,main=state)
    i <- i + 1
  }
}

#' @title Zero Substitution in 3 Repeats.
#' @description Zero value substitution of 1 repeat in 3 repeats with average.
#' @details Input dataframe, and column start and end index of each repeat, then return a set of dataframe.
#' @param protdat A DataFrame from with only numeric value.
#' @param s1 A integer of start index of repeat 1.
#' @param e1 A integer of end index of repeat 1.
#' @param s2 A integer of start index of repeat 2.
#' @param e2 A integer of end index of repeat 2.
#' @param s3 A integer of start index of repeat 3 need substitution.
#' @param e3 A integer of end index of repeat 3 need substitution.
#' @param condition_num A integer of condition numbers.
#' @param least A integer for zero substitution for repeat 1 and 2.
#' @export
data_zero_substitution_3repeats <- function(protdat,s1,e1,s2,e2,s3,e3,condition_num,least=30){
  temp_1 <- protdat[,s1:e1]
  temp_1[temp_1 == 0] <- least
  temp_2 <- protdat[,s2:e2]
  temp_2[temp_2 == 0] <- least
  temp_3 <- protdat[,s3:e3]

  temp_12 <- cbind(temp_1,temp_2)
  temp_12_norm <- data_norm_irs(temp_12,condition_num,2)

  temp_1 <- temp_12_norm[,1:condition_num]
  temp_2 <- temp_12_norm[,(condition_num + 1):(condition_num * 2)]

  i <- 1
  while(i <= condition_num){
    temp_3[,i] <- (temp_1[,i] + temp_2[,i]) / 2
    i <- i + 1
  }

  protdat[,s1:e1] <- temp_1
  protdat[,s2:e2] <- temp_2
  protdat[,s3:e3] <- temp_3

  return(protdat)
}

#' @title Scale Dataframe.
#' @description Scale dataframe by row to mean of 0 and sd of 1.
#' @details Input dataframe, then return a dataframe after scaled.
#' @param prot A DataFrame with only numeric value.
#' @export
data_scale_mat <- function(prot){
  m <- apply(prot, 1, mean)
  s <- apply(prot, 1, sd)
  prot <- (prot - m) / s
  return(prot)
}

