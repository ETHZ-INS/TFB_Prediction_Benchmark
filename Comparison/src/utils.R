getPRStats <- function(dt,
                       scores="score",
                       labels="cond",
                       models=NULL,
                       posClass="pos",
                       negClass="neg",
                       subSample=FALSE,
                       aggregate=FALSE,
                       seed=42){
  set.seed(seed)
  dt <- as.data.table(dt)
  
  if(is.null(models)){
    dt$model <- "all"
    models <- "model"
  }
  
  dt <- copy(dt)
  setnames(dt, scores, "scores")
  setnames(dt, labels, "labels")
  
  posFrac <- sum(dt$labels==posClass)/nrow(dt)
  
  setorder(dt, -scores)
  if(subSample)
  {
    dt <- dt[,.SD[sample(.N, min(5e5,.N))], by=c(models)]
  }
  
  dt[,tpr:=cumsum(labels==posClass)/sum(labels==posClass), by=c(models)]
  dt[,fpr:=cumsum(labels==negClass)/sum(labels==negClass), by=c(models)]
  dt[,fdr:=cumsum(labels==negClass)/seq_len(.N), by=c(models)]
  dt[,ppv:=cumsum(labels==posClass)/seq_len(.N), by=c(models)]
  dt[,p:=seq_len(.N), by=c(models)]
  dt[,idx:=1:.N, by=c(models)]
  
  # compute distance to top right
  dt[,dist:=sqrt((1-ppv)^2+(1-tpr)^2)]
  
  # find sparsification threshold
  dt[,min_dist:=min(dist), by=c(models)]
  dt[, thr:=data.table::first(scores[dist==min_dist]), by = c(models)]
  
  # get precision and recall at optimum
  dt[,tpr_thr:=sum(scores>=thr & labels==posClass)/sum(labels==posClass),
     by=c(models)]
  dt[,ppv_thr:=sum(scores>=thr & labels==posClass)/sum(scores>=thr),
     by=c(models)]
  
  dt[,sum_pos:=sum(labels==posClass), by=c(models)]
  dt[,sum_neg:=sum(labels==negClass), by=c(models)]
  dt <- subset(dt, sum_pos>0 & sum_neg>0)
  
  
  dt[, is_closest_0.01 := {
    sel <- which(tpr > 0.01)
    .I[ sel[ which.min(abs(tpr[sel] - 0.01)) ] ]
  }, by = models]
  dt[,ppv_0.01:=ppv[is_closest_0.01], by=c(models)]
  
  dt[, is_closest_0.05 := {
    sel <- which(tpr > 0.05)
    .I[ sel[ which.min(abs(tpr[sel] - 0.05)) ] ]
  }, by = models]
  dt[,ppv_0.05:=ppv[is_closest_0.05], by=c(models)]

  dt[, is_closest_0.1 := {
    sel <- which(tpr > 0.1)
    .I[ sel[ which.min(abs(tpr[sel] - 0.1)) ] ]
  }, by = models]
  dt[,ppv_0.1:=ppv[is_closest_0.1], by=c(models)]
  
  if(nrow(dt)>0){
    dt[,auc_pr_mod:=PRROC::pr.curve(scores.class0=scores,
                                    weights.class0=as.integer(labels),
                                    curve=FALSE)$auc.integral, # actually needs to be defined whats changing and whats not
       by=c(models)]
    
    setnames(dt, c("labels"), c(labels))
    
    if(aggregate){
      dt <- dt[,.(auc_pr_mod=data.table::first(auc_pr_mod),
                  pr=data.table::first(ppv_thr),
                  thr=data.table::first(thr),
                  recall=data.table::first(tpr_thr),
                  ppv_0.01=data.table::first(ppv_0.01),
                  ppv_0.05=data.table::first(ppv_0.05),
                  ppv_0.1=data.table::first(ppv_0.1)),
                  by=c(models)]
    }
    dt$pos_frac <- posFrac
    return(dt)
  }
  else
  {
    return(NULL)
  }
}
