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

data_norm <- function(prot_dat,condition_num,repeat_num,reference=TRUE,ref_col = c(1,7,13)){

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

data_norm_sl <- function(prot_dat){
  target <- mean(colSums(prot_dat))
  norm_facs <- target / colSums(prot_dat)
  prot_dat_sl <- sweep(prot_dat,2,norm_facs,FUN = '*')
  return(prot_dat_sl)
}

data_norm_tmm <- function(prot_dat){
  target <- calcNormFactors(prot_dat)
  prot_dat_tmm <- sweep(prot_dat,2,target,FUN='/')
  return(prot_dat_tmm)
}

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

