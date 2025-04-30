cal_w2_pval <- function(count_per_gene, meta_individual, meta_phenotype,
                        perm_num=200, unif_round_unit=0.2,
                        weight=1, shrink=FALSE){
  cur_individual=unique(meta_individual)
  phenotype=meta_phenotype[match(cur_individual,meta_individual)]

  n=length(cur_individual)

  #count_per_gene=round(count_per_gene,unif_round_unit)
  cur_range=seq(min(count_per_gene),(max(count_per_gene)+unif_round_unit),by=unif_round_unit) #cur range +1
  unif_round_unit_digit=min(BSDE:::decimalplaces(unif_round_unit),5)
  cur_range=floor(cur_range*10^unif_round_unit_digit)/10^unif_round_unit_digit
  cur_range=c(cur_range, (max(cur_range)+cur_range[2]-cur_range[1])) #cur range expand +2
  len=length(cur_range)-1
  cur_burden=matrix((rep(cur_range[1:len],times=len)-rep(cur_range[1:len],each=len))^2,len,len)
  #cur_burden=cur_burden/(10^unif_round_unit)
  cur_distr=matrix(0,nrow=len,ncol=n) #row for ranges col for individuals
  #background_image=marix(0,nrow=len,ncol=(len-1))
  colnames(cur_distr)=cur_individual
  rownames(cur_distr)=cur_range[1:len]

  for(i_n in 1:n){
    cur_ind=cur_individual[i_n]
    cur_ind_count=count_per_gene[meta_individual==cur_ind]
    cur_freq=graphics::hist(cur_ind_count,breaks=cur_range,plot=FALSE)$counts
    cur_distr[,i_n]=cur_freq/sum(cur_freq)
  }

  if(shrink){ #however, typically, we don't shrink, since it may just by chance that there are no numbers between
    #remove zero counts, to fasten the calculation.
    shrink_index=which(apply(cur_distr,1,sum)>0)
    cur_burden=cur_burden[shrink_index,shrink_index]
    cur_distr=cur_distr[shrink_index,]
    cur_range=cur_range[shrink_index]
    len=length(cur_range)-1
  }

  # `%dorng%` <- doRNG::`%dorng%`
  # `%dopar%` <- foreach::`%dopar%`
  # mcoptions <- list(preschedule = FALSE)# load balance
  #w2_res=foreach::foreach(ip=0:perm_num) %dorng% {
  w2_res <- vector("list", perm_num + 1)
  for (ip in 0:perm_num) {
  #for(ip in 0:10) {
    #print(ip)
    if(ip>0){
      cur_phenotype=phenotype[sample.int(n,n)]
    }
    if(ip==0){
      cur_phenotype=phenotype
    }
    case_distr=cur_distr[,cur_phenotype==1]
    ctrl_distr=cur_distr[,cur_phenotype==0]

    w_case=1
    w_ctrl=1
    if(length(weight)!=1){
      stopifnot("weight length must match # of subjects"=length(weight)==length(cur_phenotype))
      w_case=weight[cur_phenotype==1]
      w_ctrl=weight[cur_phenotype==0]
      w_case =  w_case / sum(w_case)
      w_ctrl =  w_ctrl / sum(w_ctrl)
    }

    py <- reticulate::import_main()
    py$case_wass=NULL
    py$ctrl_wass=NULL
    py$case_distr=reticulate::r_to_py(case_distr)
    py$ctrl_distr=reticulate::r_to_py(ctrl_distr)
    py$w_case=reticulate::r_to_py(w_case)
    py$w_ctrl=reticulate::r_to_py(w_ctrl)

    reticulate::py_run_string("case_wass = cal_bary_wass(case_distr,w=w_case)")
    reticulate::py_run_string("ctrl_wass = cal_bary_wass(ctrl_distr,w=w_ctrl)")
    case_bc=py$case_wass
    ctrl_bc=py$ctrl_wass
    names(case_bc)=cur_range[1:len]
    names(ctrl_bc)=cur_range[1:len]
    #calculate w2

    w2=Barycenter::Greenkhorn(as.matrix(case_bc),(as.matrix(ctrl_bc)),costm=abs(cur_burden), lambda=0.01)$Distance
    #w2
    w2_res[[ip + 1]] <- list(w2,case_bc,ctrl_bc)
  }
  #w2_ob=w2_res[[1]]
  #w2_perm=sapply(w2_res,unlist)[-1]
  w2_ob=w2_res[[1]][[1]]
  w2_perm=sapply(w2_res[-1],function(x)unlist(x[[1]]))
  case_bc_ob=w2_res[[1]][[2]]
  ctrl_bc_ob=w2_res[[1]][[3]]

  # pval=mean(w2_perm>=w2_ob,na.rm=TRUE)
  pval=mean(c(w2_perm, w2_ob)>=w2_ob,na.rm=TRUE) # add one in both numerator and denominator
  return(list(pval,case_bc_ob,ctrl_bc_ob))
}
