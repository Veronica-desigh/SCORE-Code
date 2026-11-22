# =============================================================================
# Function_13_Report_tables_M01_v6.2.R | 2026-11-10
# 只从M01及S01保存的统计摘要生成27张补充表；保留用户v6.2 Word原生版式。
# 依赖：data.table, xml2, zip（Excel审计导出另需openxlsx）。
# 数据接口：S0/SH/ST/SHT为上游内部码；S0/SN/ST/SNT为Word显示标签。
# 修改方法学和数值前必须更新上游模型；报表层绝不重算CIF、插补或模型方差。
# =============================================================================
rp13_word_ns <- c(w='http://schemas.openxmlformats.org/wordprocessingml/2006/main')

rp_msg <- function(...) cat(sprintf('[%s] ',format(Sys.time(),'%H:%M:%S')), ..., '\n', sep='')

rp_dir <- function(p) {dir.create(p,showWarnings=FALSE,recursive=TRUE); invisible(p)}

rp_md5 <- function(p) if(file.exists(p)) unname(as.character(tools::md5sum(p))) else NA_character_

rp_require <- function(p, label='输入') if(!file.exists(p)) stop(label,'不存在：',p,call.=FALSE) else invisible(p)

rp_assert <- function(x, cols, label='数据') {miss<-setdiff(cols,names(x));if(length(miss))stop(label,'缺失字段：',paste(miss,collapse=', '),call.=FALSE); invisible(TRUE)}

rp_guard <- function(path){if(file.exists(path)&&!isTRUE(getOption('M01_REPORT_OVERWRITE',TRUE)))stop('输出文件已存在，overwrite=FALSE：',path);invisible(TRUE)}

rp_save_csv <- function(x, path) {rp_dir(dirname(path));rp_guard(path);data.table::fwrite(x,path,na='',bom=TRUE);invisible(path)}

rp_file <- function(analysis_dir, board, file, layer='pooled') file.path(analysis_dir,board,layer,paste0(file,'.csv'))

rp_read <- function(path, required=NULL, validate_status=FALSE) {
  rp_require(path);d<-data.table::fread(path,encoding='UTF-8')
  if(!is.null(required))rp_assert(d,required,path)
  if(validate_status && 'status'%in%names(d) && any(is.na(d$status)|d$status!='estimable'))
    stop('存在不可估计的 pooled 行：',path,call.=FALSE)
  d
}

rp_check_config <- function(cfg, style) {
  rp_assert(as.list(cfg),c('m','times','confidence_level'),'07 study_config')
  if(!all(style$horizons %in% cfg$times))stop('报告时间点不在 07 冻结 cfg$times 中')
  if(!(style$primary_horizon %in% cfg$times))stop('主时间点未由上游分析估计')
  if(!identical(sort(style$scenario_order),sort(c('S0','SH','ST','SHT'))))stop('RM 情景名称未匹配 07')
  invisible(TRUE)
}

rp_load_config <- function(p,style) {
  rp_require(p,'study_config.rds'); cfg<-readRDS(p);rp_check_config(cfg,style);cfg
}

rp_preflight <- function(paths) {
  data.table(source=names(paths),path=as.character(paths),exists=file.exists(unlist(paths)),
             bytes=as.numeric(file.info(unlist(paths))$size),
             md5=vapply(paths,rp_md5,character(1)))
}

rp_num <- function(x,digits=2L,group=FALSE) {
  v<-rep('—',length(x));ok<-is.finite(as.numeric(x));val<-as.numeric(x[ok]);
  val[abs(val)<0.5*10^(-digits)]<-0
  v[ok]<-formatC(val,format='f',digits=digits,big.mark=if(group)','else '',decimal.mark='.')
  v
}

rp_n <- function(x)rp_num(x,0L,TRUE)

rp_pct <- function(x,digits=1L,prob=TRUE)rp_num(if(prob)100*x else x,digits)

rp_ci <- function(estimate,lcl,ucl,digits=2L,scale=100) {
  if(length(estimate)!=length(lcl)||length(lcl)!=length(ucl))stop('CI 长度不一致')
  a<-rp_num(estimate*scale,digits);b<-rp_num(lcl*scale,digits);c<-rp_num(ucl*scale,digits)
  ifelse(a=='—'|b=='—'|c=='—','—',paste0(a,' (',b,'–',c,')'))
}

rp_count_pct <- function(n,proportion,digits=1L) {
  a<-rp_n(n);b<-rp_pct(proportion,digits);ifelse(a=='—'|b=='—','—',paste0(a,' (',b,')'))
}

rp_pvalue <- function(x,digits=3L) {
  ifelse(!is.finite(x),'—',ifelse(x<10^(-digits),paste0('<',formatC(10^(-digits),format='f',digits=digits)),rp_num(x,digits)))
}

rp_row <- function(d,filters,required=TRUE,unique=TRUE) {
  z<-as.data.table(d)
  for(nm in names(filters)){if(!nm%in%names(z))stop('来源缺少筛选字段：',nm);v<-filters[[nm]];z<-z[get(nm)%in%v]}
  if(required && nrow(z)==0L)stop('缺少目标结果：',paste(names(filters),filters,sep='=',collapse='; '),call.=FALSE)
  if(unique && nrow(z)>1L)stop('目标结果重复：',paste(names(filters),filters,sep='=',collapse='; '),call.=FALSE)
  z
}

rp_check_ci <- function(d,label) {
  if(is.null(d)||!all(c('estimate','lcl','ucl')%in%names(d)))return(invisible(NULL))
  z<-d[is.finite(estimate)&is.finite(lcl)&is.finite(ucl)]
  if(nrow(z) && any(z$lcl>z$estimate+1e-10 | z$estimate>z$ucl+1e-10))stop('置信区间不包含点估计：',label)
  invisible(TRUE)
}

rp_check_unique <- function(d,keys,label) {
  rp_assert(d,keys,label);if(anyDuplicated(d,by=keys))stop('结果键存在重复：',label,' / ',paste(keys,collapse=', '));TRUE
}

rp13_word_rows <- function(tbl) xml2::xml_find_all(tbl, './w:tr', ns = rp13_word_ns)

rp13_word_cells <- function(row) xml2::xml_find_all(row, './w:tc', ns = rp13_word_ns)

rp13_word_text <- function(cell) {
  t <- xml2::xml_find_all(cell, './/w:t', ns = rp13_word_ns)
  if(!length(t)) return('')
  paste(xml2::xml_text(t), collapse = '')
}

rp13_word_put <- function(cell, value) {
  # 只替换目标单元格文字节点；保留全部 tblPr/tcPr、单元格底纹、边框、对齐和字号。
  if(length(value) != 1L || is.na(value)) return(FALSE)
  t <- xml2::xml_find_all(cell, './/w:t', ns = rp13_word_ns)
  if(!length(t)){if(value=='—')return(FALSE);p<-xml2::xml_find_first(cell,'./w:p',ns=rp13_word_ns);if(inherits(p,'xml_missing'))p<-xml2::xml_add_child(cell,'w:p');r<-xml2::xml_add_child(p,'w:r');tt<-xml2::xml_add_child(r,'w:t');xml2::xml_text(tt)<-value;return(TRUE)}
  xml2::xml_text(t[[1L]]) <- as.character(value)
  if(length(t) > 1L) for(j in seq.int(2L, length(t))) xml2::xml_text(t[[j]]) <- ''
  TRUE
}

rp13_verify_stage10_manifest <- function(analysis_dir,manifest_path){
  q<-rp_read(manifest_path,c('source','relative_path','md5','bytes','contract','m'))
  required<-c('ht_counts','ht_risk','rm_domains','rm_overlap_compat',
    'risk_requirements','intervention_requirements','hvh_risk','hvh_paired',
    'reclassification')
  if(!identical(sort(q$source),sort(required))||
     anyDuplicated(q$source)||any(q$contract!='classification_unified_v6')||
     uniqueN(q$m)!=1L||any(is.na(q$md5)))
    stop('13: invalid unified step-10 export manifest; rerun 10 fully')
  if(any(grepl('\\',q$relative_path,fixed=TRUE)) ||
     any(startsWith(q$relative_path,'/')) || any(grepl('..',q$relative_path,fixed=TRUE)))
    stop('13: manifest relative_path contains an unsafe path segment')
  paths<-file.path(analysis_dir,q$relative_path)
  if(any(!file.exists(paths)))stop('13: unified step-10 output missing: ',
    paste(q$relative_path[!file.exists(paths)],collapse=', '))
  actual<-vapply(paths,rp_md5,character(1))
  sizes<-as.numeric(file.info(paths)$size)
  if(any(is.na(sizes)|sizes!=q$bytes|actual!=q$md5))
    stop('13: step-10 CSV mixed/stale versions or incomplete export; rerun 10 fully. Files: ',
      paste(q$relative_path[is.na(sizes)|sizes!=q$bytes|actual!=q$md5],collapse=', '))
  invisible(TRUE)
}

rp61_build <- function(src,blueprint,style) {
  out<-list();s<-style
  put<-function(t,r,c,value,key){
    if(length(value)!=1L||is.na(value)||grepl('TBC|x[.]xx',value,ignore.case=TRUE))stop('Unfilled cell ',t,'/',r,'/',c)
    out[[length(out)+1L]]<<-data.table(table=t,row=r,col=c,value=as.character(value),source_key=key)
  }
  one<-function(name,filters){
    z<-rp_row(src[[name]],filters)
    if('status'%in%names(z)&&!z$status%in%c('estimable','not_applicable'))stop('Invalid estimate: ',name,' ',paste(filters,collapse='/'))
    z
  }
  count<-function(z,n='n',p='proportion')rp_count_pct(z[[n]],z[[p]],2)
  ci<-function(name,f){z<-one(name,f);rp_ci(z$estimate,z$lcl,z$ucl,2)}
  metrics<-c('risk_natural','risk_2019','risk_2025','rd_2019_vs_natural','rd_2025_vs_natural','rd_2025_vs_2019')
  scope_keys<-setNames(names(s$labels$scope),unname(s$labels$scope))
  tr_keys<-setNames(names(s$labels$transition),unname(s$labels$transition))
  scenario_keys<-c(S0='S0',SN='SH',ST='ST',SNT='SHT')
  # 6: actual method identifiers, keeping source definitions and post-MI text.
  method_vars<-c(`1`='ethnicity',`2`='education',`3`='household_income',`4`='employment_status',
    `5`='townsend_deprivation_index',`6`='rm_source_inactive',`7`='alcohol_score',`8`='diet_score',
    `9`='bmi',`10`='waist_circumference',`11`='triglycerides',`12`='rm_source_family_proxy',
    `13`='chronic_liver_disease',`14`='chronic_respiratory_disease',`15`='rm_source_inflammatory',
    `16`='rm_source_psychiatric',`17`='rm_source_af',`18`='rm_source_stressor_proxy',
    `19`='rm_source_hiv',`20`='rm_source_sleep_apnoea',`21`='rm_source_fatty_liver',
    `22`='rm_source_premature_menopause',`23`='rm_source_hdp',`24`='rm_lpa_status',`25`='rm_source_crp')
  for(r in names(method_vars)){
    v<-method_vars[[r]];z<-rp_row(src$methods,list(variable=v),required=FALSE)
    if(r=='22'){
      # 完全忠实于v6.2 Word：未插补的绝经变量作为固定观察派生量。
      if(nrow(z) && !(as.character(z$method) %in% c('not_imputed','fixed_observed','Complete','')))
        stop('MISSING_SOURCE: premature menopause 插补记录与Word定义冲突：',z$method)
      put(4,22L,2,'Not imputed; fixed observed derivation',
          'methods/rm_source_premature_menopause/fixed_observed');next
    }
    if(!nrow(z))stop('MISSING_SOURCE: actual MI method: ',v)
    method<-z$method
    if(is.na(method)||method=='')method<-'Complete'
    put(4,as.integer(r),2,rp52_method_label(method),paste('methods',v,sep='/'))
  }
  # 7: explicit source variable and category keys; no dependence on old table numbers.
  parent<-c(`1`='education',`5`='household_income',`11`='employment_status',`15`='townsend_deprivation_index',
    `16`='rm_deprivation',`19`='ethnicity',`22`='rm_inactive',`25`='alcohol_score',`28`='diet_score',
    `31`='bmi',`32`='waist_circumference',`33`='rm_obesity',`36`='triglycerides',
    `37`='rm_family_proxy',`40`='rm_inflammatory',`43`='rm_psychiatric',`46`='rm_af',`49`='rm_stressor_proxy',
    `52`='rm_hiv',`55`='rm_sleep_apnoea',`58`='rm_fatty_liver',`61`='rm_hdp',`64`='rm_lpa',
    `67`='rm_source_crp',`68`='rm_crp')
  cont<-c('townsend_deprivation_index','bmi','waist_circumference','triglycerides','rm_source_crp')
  levels<-list(education=c('low','medium','high'),household_income=c('less_than_18000','18000_to_30999','31000_to_51999','52000_to_100000','greater_than_100000'),employment_status=c('employed_or_self_employed','retired','other_not_employed'))
  raw<-src$baseline
  raw<-raw[table==6 | (table==3 & !variable%in%raw[table==6,unique(variable)])]
  d<-unique(raw[,.(variable,kind,level,imputation_id,source,n_total,n_observed,n_missing,count,proportion,q1,median,q3)])
  desc<-function(v,kind,level=NULL,source){
    kk<-kind;ss<-source;lev_arg<-level
    q<-d[variable==v & kind==kk & source==ss]
    if(!nrow(q) && kk=='missing')q<-unique(d[variable==v & source==ss,.(imputation_id,n_total,n_observed,n_missing,count=NA_real_,proportion=NA_real_,q1=NA_real_,median=NA_real_,q3=NA_real_)])
    if(!is.null(level)){
      all_levels<-q;q<-q[as.character(q$level)==lev_arg]
      if(!nrow(q) && lev_arg=='0' && nrow(all_levels) && all(as.character(all_levels$level)=='1')){
        q<-copy(all_levels);q[,`:=`(count=n_observed-count,proportion=1-proportion)]
      }
    }
    if(!nrow(q)||anyDuplicated(q$imputation_id))stop('Missing/duplicate baseline summary ',v,' ',kind,' ',source)
    q[,lapply(.SD,mean),.SDcols=c('n_total','n_observed','n_missing','count','proportion','q1','median','q3')]
  }
  iq<-function(q)paste0(rp_num(q$median,2),' (',rp_num(q$q1,2),'–',rp_num(q$q3,2),')')
  v<-NULL;idx<-0L
  for(r in 1:70){
    label<-blueprint[table==9 & row==r & col==0,text]
    if(startsWith(label,'Panel'))next
    if(as.character(r)%in%names(parent)){
      v<-parent[[as.character(r)]];idx<-0L;k<-if(v%in%cont)'continuous'else 'missing'
      before<-desc(v,k,source='observed')
      put(9,r,2,rp_count_pct(before$n_missing,before$n_missing/before$n_total,2),paste('baseline',v,'missing',sep='/'))
      if(k=='continuous'){
        put(9,r,1,iq(before),paste('baseline',v,'observed',sep='/'))
        put(9,r,3,iq(desc(v,k,source='completed')),paste('baseline',v,'completed',sep='/'))
      }
    }else{
      idx<-idx+1L;lev<-if(v%in%names(levels))levels[[v]][idx]else as.character(idx-1L)
      for(c in c(1L,3L)){
        source<-if(c==1)'observed'else 'completed';q<-desc(v,'category',lev,source)
        put(9,r,c,count(q,'count'),paste('baseline',v,lev,source,sep='/'))
      }
    }
  }
  # 8: percentages, mean (range) across MI, in template category order.
  lv<-list(education=c('low','medium','high'),household_income=levels$household_income,
    employment_status=c('other_not_employed','retired','employed_or_self_employed'))
  for(spec in list(list(v='education',r=2:4),list(v='household_income',r=6:10),list(v='employment_status',r=12:14))){
    for(j in seq_along(spec$r))for(c in 1:2){
      g<-if(c==1)'low'else 'high';q<-src$lca[variable==spec$v & value==lv[[spec$v]][j] & ses_class==g]
      ids<-sort(unique(src$lca$imputation_id))
      if(anyDuplicated(q$imputation_id)||length(ids)!=unique(src$main$m))stop('LCA source incomplete')
      q<-merge(data.table(imputation_id=ids),q,by='imputation_id',all.x=TRUE);q[is.na(percent),percent:=0]
      put(10,spec$r[j],c,paste0(rp_num(mean(q$percent),2),' (',rp_num(min(q$percent),2),'–',rp_num(max(q$percent),2),')'),paste('lca',spec$v,lv[[spec$v]][j],g,sep='/'))
    }
  }
  # 9--19 contexts are read from the template labels, never its example values.
  for(t in 11:21){
    sc<-'overall';year<-'2019';h<-15;comparison<-NULL
    for(r in sort(unique(blueprint[table==t & row>0,row]))){
      label<-blueprint[table==t & row==r & col==0,text]
      next_label<-blueprint[table==t & row==r+1L & col==0,text]
      is_group<-isTRUE(blueprint[table==t & row==r & col==0,merged]) || (label=='Hypertension' && identical(next_label,'No hypertension'))
      is_scope<-label%in%names(scope_keys) && !is_group
      if(label%in%names(scope_keys) && !(label=='Hypertension' && identical(next_label,'No hypertension')))sc<-scope_keys[[label]]
      if(grepl('^2019',label))year<-'2019'
      if(grepl('^2025',label))year<-if(grepl('2025−2019',label,fixed=TRUE))'2025-2019'else '2025'
      if(label=='10-year')h<-10
      if(label=='15-year')h<-15
      if(t==21 && grepl(' vs ',label,fixed=TRUE)){
        comparison<-gsub(' vs ','_vs_',label,fixed=TRUE)
        comparison<-gsub('SNT','SHT',comparison,fixed=TRUE);comparison<-gsub('SN','SH',comparison,fixed=TRUE)
      }
      f<-list(scope=sc,scenario='S0');key<-paste(t,label,sc,year,h,comparison,sep='/')
      emit<-function(c,v)put(t,r,c,v,key)
      if(t==11 && is_scope){
        ct<-one('transitions',c(f,list(transition=s$transition_order[1])))
        emit(1,rp_n(ct$n_scope))
        for(j in 1:2){q<-one('requirements',c(f,list(guideline=c('2019','2025')[j],metric='n_intervention')));emit(j+1,count(q,'MI_mean','MI_mean_proportion'))}
        for(j in 1:6)emit(j+3,count(one('transitions',c(f,list(transition=s$transition_order[j])))))
      }
      if(t%in%c(12,14,18,19) && label%in%names(tr_keys)){
        tr<-tr_keys[[label]];ft<-c(f,list(transition=tr));q<-one('transitions',ft)
        if(t==12){
          emit(1,count(q));for(j in 1:5)emit(j+1,rp_num(q[[c('mean_baseline_ldl','T_2019','T_2025','mean_D_2019','mean_D_2025')[j]]],2))
        }
        if(t==14){
          emit(1,count(q))
          # v6.2 Word在≥60岁新增干预(n=2)行，仅保留人数，风险与RD均为—。
          # 属于样表明确的稀疏单元格展示规则；绝不写回不稳定的示例估计。
          sparse_suppressed <- sc=='age_ge_60' && tr=='newly_intervention'
          if(!sparse_suppressed){
            emit(2,ci('transition_observed',c(ft,list(time=15,estimand='observed_transition_risk'))))
            mm<-c('risk_natural','rd_2019_vs_natural','rd_2025_vs_natural','rd_2025_vs_2019')
            for(j in 1:4)emit(j+2,ci('transition_effects',c(ft,list(time=15,estimand=mm[j]))))
            emit(7,ci('contribution',c(ft,list(time=15,estimand='absolute_RD_contribution'))))
          }
        }
        if(t==18)for(j in 1:4){ft$scenario<-s$scenario_order[j];emit(j,count(one('transitions',ft)))}
        if(t==19)for(j in 1:4){ft$scenario<-s$scenario_order[j];ft$scope<-'overall';ft$guideline<-year;emit(j,rp_num(one('mean_reduction',ft)$mean_D,2))}
      }
      if(t==14 && label=='Total'){
        q<-src$transitions[scope==sc & scenario=='S0'];emit(1,rp_count_pct(sum(q$n),1,2))
        emit(2,ci('observed',c(f,list(time=15,estimand='observed_scope_risk'))))
        for(j in 1:4)emit(j+2,ci('s0',c(f,list(time=15,estimand=c('risk_natural','rd_2019_vs_natural','rd_2025_vs_natural','rd_2025_vs_2019')[j]))))
        emit(7,ci('s0',c(f,list(time=15,estimand='rd_2025_vs_2019'))))
      }
      if(t==13 && is_scope)for(j in 1:6)emit(j,ci('s0',c(f,list(time=h,estimand=metrics[j]))))
      if(t==13 && label%in%c('Age','Sex','SES','Hypertension') && (label!='Hypertension'||is_group)){
        v<-c(Age='age',Sex='sex',SES='ses',Hypertension='hypertension')[[label]]
        ff<-list(scope=paste0('heterogeneity_',v),time=h)
        if(t==13){ff$scenario<-'S0';ff$estimand<-'delta_rd_2025_vs_2019';nm<-'s0_het'}else{ff$rm_comparison<-comparison;nm<-'rm_het'}
        emit(if(t==13)7 else 4,rp_pvalue(one(nm,ff)$p_value,2))
      }
      if(t==15){
        domains<-c(`3`='rm_deprivation',`4`='rm_stress',`5`='rm_psychiatric',`6`='rm_south_asian',
          `8`='rm_obesity',`9`='rm_inactive',`10`='rm_sleep_apnoea',`11`='rm_fatty_liver',
          `13`='rm_family',`14`='rm_lpa',`16`='rm_inflammatory',`17`='rm_hiv',`18`='rm_crp',
          `20`='rm_premature_menopause',`21`='rm_hdp',`23`='rm_af',`26`='any_H',`27`='any_T',`28`='any_HT')
        for(j in 1:2){
          y<-c('2019','2025')[j]
          if(r==30L){
            z<-one('s01_domains',list(guideline=y,scenario_id='SNT_DOMAIN_GE2',metric='at_least_2'))
            put(t,r,j+1L,count(z),paste('s01_domains',y,'SNT_DOMAIN_GE2','at_least_2',sep='/'))
          }else if(as.character(r)%in%names(domains)){
            v<-domains[[as.character(r)]]
            if(blueprint[table==t & row==r & col==j+1L,text]!='—'){
              z<-one('domains',list(scope='overall',guideline=y,domain=v,level=1))
              put(t,r,j+1L,count(z),paste('domains',y,v,sep='/'))
            }
          }
        }
      }
      if(t==21 && is_scope){
        dimensions<-c(age_lt_60='age',female='sex',high_ses='ses',no_hypertension='hypertension')
        if(sc%in%names(dimensions)){
          v<-dimensions[[sc]]
          z<-one('rm_het',list(scope=paste0('heterogeneity_',v),time=15,rm_comparison=comparison))
          put(t,r,4L,rp_pvalue(z$p_value,2),paste('rm_het',comparison,v,15,sep='/'))
        }
      }
      if(t==16 && (startsWith(label,'N+')||startsWith(label,'N−')||startsWith(label,'Any '))){
        cell<-if(startsWith(label,'Any '))'all'else c('N+T−'='H1T0','N−T+'='H0T1','N+T+'='H1T1')[[label]]
        q<-one('new_counts',list(guideline=year,scenario='SHT',cell=cell))
        emit(1,rp_count_pct(q$n_upgrade,q$n_upgrade/q$n_initial_M,2))
        if(cell=='all')obs<-ci('upgraded',list(scope=paste(year,'overall',sep='__'),scenario='SHT',time=15,estimand='observed_upgraded_risk'))
        else obs<-ci('ht_risk',list(scope='overall',scenario='S0',guideline=year,time=15,h=as.integer(substr(cell,2,2)),t=as.integer(substr(cell,4,4)),estimand='HT_stratum_risk'))
        emit(2,obs);emit(3,rp_count_pct(q$n_new,q$proportion_new,2))
        emit(4,ci('new_risk',list(scope='overall',scenario='SHT',guideline=year,cell=cell,time=15)))
      }
      if(t==17 && label%in%names(scenario_keys)){
        k<-scenario_keys[[label]];q<-one('requirements',list(scope='overall',scenario=k,guideline=year,metric='n_intervention'))
        emit(1,count(q,'MI_mean','MI_mean_proportion'))
        if(k=='S0'){emit(2,'—');emit(3,'—')}else{
          q<-one('new_counts',list(guideline=year,scenario=k,cell='all'))
          emit(2,rp_count_pct(q$n_new,q$proportion_new,2))
          emit(3,ci('new_risk',list(scope='overall',scenario=k,guideline=year,cell='all',time=15)))
        }
      }
      if(t==20 && label%in%names(scenario_keys)){
        emit(1,ci('observed',c(f,list(time=15,estimand='observed_scope_risk'))))
        f$scenario<-scenario_keys[[label]]
        for(j in 1:6)emit(j+1,ci('all_effects',c(f,list(time=15,estimand=metrics[j]))))
      }
      if(t==21 && is_scope)for(j in 1:3){
        emit(j,ci('rm_effects',list(scope=sc,rm_comparison=comparison,time=15,estimand=c('additional_reduction_2019','additional_reduction_2025','additional_reduction_difference_2025_vs_2019')[j])))
      }
    }
  }
  rbindlist(out)
}

rp51_main_cells <- function(d){
  out<-list();put<-function(r,c,v,key)out[[length(out)+1L]]<<-data.table(table=1,row=r,col=c,value=v,source_key=key)
  groups<-c('overall','no_observed_ascvd','observed_ascvd')
  for(c in 1:3){
    g<-groups[c];z<-d[scope==g];n<-unique(z$n_total)
    if(nrow(z)!=16L||length(n)!=1L||anyDuplicated(z$row))stop('Main Table 1 schema failed')
    put(0,c,paste0(c('Overall (N=','No observed ASCVD (n=','Observed ASCVD (n=')[c],rp_n(n),')'),paste('main',g,'N',sep='/'))
    for(r in 1:16){q<-z[row==r]
      v<-if(q$kind=='continuous')paste0(rp_num(q$median,2),' (',rp_num(q$q1,2),'–',rp_num(q$q3,2),')')else rp_count_pct(q$count,q$proportion,2)
      put(r,c,v,paste('main',g,q$variable,sep='/'))
    }
  }
  rbindlist(out)
}

rp52_method_label <- function(method){
 labels<-c(logreg='Binary logistic regression',polr='Ordinal logistic regression (proportional odds)',
  polyreg='Multinomial logistic regression',pmm='Predictive mean matching',
  rm_tg='Assay-bound-constrained log-normal regression',rm_crp='Assay-bound-constrained log-normal regression',
  rm_female_hdp='Sex-restricted logistic regression (women only)',
  not_imputed='Not imputed; fixed observed derivation',
  rm_south_asian='Ethnicity-restricted binary logistic regression',Complete='Complete')
 if(!method%in%names(labels))stop('Unmapped MI method: ',method)
 unname(labels[[method]])
}

rp52_method_audit <- function(methods){
 q<-methods[variable%in%c('rm_source_premature_menopause','rm_source_hdp'),.(variable,method)]
 if(!'rm_source_premature_menopause'%in%q$variable)
   q<-rbind(q,data.table(variable='rm_source_premature_menopause',method='not_imputed'))
 q[,template_method:=ifelse(variable=='rm_source_premature_menopause',
   'Not imputed; fixed observed derivation','Sex-restricted logistic regression (women only)')]
 q[,actual_method:=vapply(method,function(v)
   if(v%in%c('not_imputed','fixed_observed',''))'Not imputed; fixed observed derivation'
   else rp52_method_label(v),character(1))]
 q[,matches_template:=actual_method==template_method]
 q
}

rp61_sensitivity_cells <- function(src, blueprint, style) {
  checked_unique <- function(data, keys) {
    data <- unique(data)
    rp_check_unique(data, keys, 'S01 full precision source')
    data
  }
  counts <- checked_unique(src$s01_counts, c('scenario_id','guideline','metric'))
  transitions <- checked_unique(src$s01_transitions[grouping==6], c('scenario_id','group'))
  effects <- checked_unique(src$s01_effects[scope=='overall' & time==15], c('scenario','estimand'))
  dictionary <- src$s01_dictionary
  rp_check_unique(dictionary,'rm_scenario','S01 scenario dictionary')
  years <- c('2019','2025')
  cells <- list()
  put <- function(table,row,col,value,key,applicable=TRUE) {
    if(length(value)!=1L || is.na(value) || (applicable && value=='—'))
      stop('MISSING_SOURCE: ',table,'/',row,'/',col,' ',key)
    cells[[length(cells)+1L]] <<- data.table(table=table,row=row,col=col,value=value,
      source_key=key,status=if(applicable)'READY'else 'NOT_APPLICABLE')
  }
  for(t in 22:26) {
    year <- NULL
    for(r in sort(unique(blueprint[table==t & row>0,row]))) {
      label <- blueprint[table==t & row==r & col==0,text]
      if(t==22 && grepl('^(2019|2025)',label)) year <- substr(label,1L,4L)
      scenario_label <- if(t==22) blueprint[table==t & row==r & col==1,text] else label
      if(!length(scenario_label) || !nzchar(scenario_label) || grepl('^Analysis [1-7]$',scenario_label)) next
      id <- rp_row(dictionary,list(rm_scenario=scenario_label))$scenario_id
      applicable_years <- setdiff(years,as.character(counts[scenario_id==id & metric=='not_applicable',guideline]))
      emit_count <- function(y,metric,c) {
        ok <- y%in%applicable_years
        key <- paste('s01_counts',id,y,metric,sep='/')
        value <- if(ok) {
          q <- rp_row(counts,list(scenario_id=id,guideline=y,metric=metric))
          rp_count_pct(q$n,q$value,2)
        } else '—'
        put(t,r,c,value,key,ok)
      }
      if(t==22) {
        metrics <- c('score_based_high_very_high','rm_upgraded','final_high_very_high')
        for(j in seq_along(metrics)) emit_count(year,metrics[j],j+1L)
      }
      if(t==23) for(j in seq_along(years)) {
        emit_count(years[j],'upgrades_lost_vs_snt',2L*j-1L)
        emit_count(years[j],'upgrades_lost_vs_snt_initial_moderate',2L*j)
      }
      if(t==24) for(j in seq_along(years)) emit_count(years[j],'intervention_eligible',j)
      if(t==25) for(j in seq_along(style$transition_order)) {
        group <- style$transition_order[j]
        ok <- length(applicable_years)==2L
        value <- if(ok) {
          q <- rp_row(transitions,list(scenario_id=id,group=group))
          rp_count_pct(q$n,q$proportion,2)
        } else '—'
        put(t,r,j,value,paste('s01_transitions',id,'6',group,sep='/'),ok)
      }
      if(t==26) {
        metrics <- c('risk_2019','risk_2025','rd_2025_vs_2019','delta_rd_vs_s0','delta_rd_vs_snt')
        for(j in seq_along(metrics)) {
          ok <- if(j<=2L) years[j]%in%applicable_years else length(applicable_years)==2L
          value <- if(ok) {
            q <- rp_row(effects,list(scenario=id,estimand=metrics[j]))
            if(q$status!='estimable' || q$m_valid!=q$m) stop('Unusable S01 estimate: ',id,'/',metrics[j])
            rp_ci(q$estimate,q$lcl,q$ucl,2)
          } else '—'
          put(t,r,j,value,paste('s01_effects',id,'overall','15',metrics[j],sep='/'),ok)
        }
      }
    }
  }
  rbindlist(cells)
}

rp62_evalues_cells <- function(e, style) {
  rp_assert(e,c('scope','scenario','time','estimand','estimate','lcl','ucl',
    'evalue_point','ci_nearest_null','evalue_ci'),'09 E-values pooled')
  names <- c('rr_2019_vs_natural','rr_2025_vs_natural','rr_2025_vs_2019')
  o <- list()
  for (i in seq_along(names)) {
    z <- rp_row(e,list(scope='overall',scenario='S0',
      time=style$primary_horizon,estimand=names[i]))
    if('status' %in% names(z) && !all(z$status=='estimable'))
      stop('MISSING_SOURCE: Table 27 estimate not estimable: ',names[i])
    vals <- c(rp_ci(z$estimate,z$lcl,z$ucl,style$digits$rr,scale=1),
      rp_num(z$evalue_point,style$digits$rr),
      rp_num(z$ci_nearest_null,style$digits$rr),
      rp_num(z$evalue_ci,style$digits$rr))
    if(any(vals=='—'))stop('MISSING_SOURCE: Table 27 incomplete E-values: ',names[i])
    for(j in seq_along(vals))o[[length(o)+1L]] <- data.table(
      table=27L,row=as.integer(i),col=as.integer(j),value=vals[j],
      source_key=paste('evalues',names[i],c('RR_CI','point','nearest_null','CI')[j],sep='/'),
      status='READY')
  }
  rbindlist(o,fill=TRUE)
}

rp62_source_candidates <- function(t) {
  m <- list(
    `4`='methods',`9`='baseline',`10`='lca',
    `11`=c('requirements','transitions'),`12`='transitions',
    `13`=c('s0','s0_het'),
    `14`=c('transitions','transition_observed','transition_effects','contribution','observed','s0'),
    `15`=c('domains','s01_domains'),
    `16`=c('new_counts','ht_risk','new_risk','upgraded'),
    `17`=c('requirements','new_counts','new_risk'),
    `18`='transitions',`19`='mean_reduction',
    `20`=c('observed','all_effects'),`21`=c('rm_effects','rm_het'),
    `22`='s01_counts',`23`='s01_counts',`24`='s01_counts',
    `25`='s01_transitions',`26`='s01_effects',`27`='evalues')
  q <- m[[as.character(t)]]
  if(is.null(q))character(0) else q
}

rp62_blueprint_read <- function(path) {
  rp_require(path,'v6.2蓝图CSV')
  b <- rp_read(path,required=c('table','row','col','text','kind','merged','demo_numeric'))
  b[,`:=`(table=as.integer(table),row=as.integer(row),col=as.integer(col),
         text=as.character(text),merged=as.logical(merged),
         demo_numeric=as.logical(demo_numeric))]
  b[is.na(text),text:='']
  if(nrow(b)==0L||anyNA(b[,.(table,row,col,merged,demo_numeric)])||
     anyDuplicated(b,by=c('table','row','col')))
    stop('v6.2蓝图存在空键、缺失或重复格')
  actual<-sort(unique(b$table))
  if(!identical(actual,seq_len(27L)))
    stop('v6.2蓝图必须包含S01–S27，实际：',paste(actual,collapse=','))
  expected_rows <- c(8,11,34,26,17,18,23,42,71,15,14,64,29,73,31,
                     11,11,64,22,46,41,143,74,78,78,74,4)
  expected_cols <- c(4,3,4,4,4,4,4,6,4,3,10,7,8,8,4,5,4,5,5,8,5,5,5,3,7,6,5)
  for(k in seq_len(27L)) {
    z<-b[table==k]
    if(max(z$row)!=expected_rows[k]-1L || max(z$col)!=expected_cols[k]-1L ||
      uniqueN(z$row)!=expected_rows[k])
      stop('v6.2蓝图表形状不匹配：S',sprintf('%02d',k))
  }
  b
}

rp62_word_cell <- function(tbl, row0, col0) {
  rows<-xml2::xml_find_all(tbl,'./w:tr',ns=rp13_word_ns)
  if(row0+1L>length(rows)||row0<0L)stop('Word row越界：',row0)
  tcs<-rp13_word_cells(rows[[row0+1L]])
  pos<-0L
  for(tc in tcs) {
    span_node<-xml2::xml_find_first(tc,'./w:tcPr/w:gridSpan',ns=rp13_word_ns)
    span<-if(inherits(span_node,'xml_missing'))1L else {
      attrs<-xml2::xml_attrs(span_node)
      v<-suppressWarnings(as.integer(unname(attrs)[1L]))
      if(is.na(v))stop('Word表格有非法gridSpan')
      v
    }
    if(col0==pos) return(list(cell=tc,span=span))
    if(col0>pos && col0<pos+span)
      stop('尝试填入合并单元格的内部网格列：',row0,'/',col0)
    pos<-pos+span
  }
  stop('Word column越界：',row0,'/',col0)
}

rp62_template_check <- function(doc, blueprint) {
  tabs<-xml2::xml_find_all(doc,'./w:body/w:tbl',ns=rp13_word_ns)
  if(length(tabs)!=27L)stop('v6.2模板必须有27个表；当前为',length(tabs))
  for(k in seq_len(27L)) {
    tr<-rp13_word_rows(tabs[[k]])
    b<-blueprint[table==k]
    if(length(tr)!=max(b$row)+1L)
      stop('表格S',k,'行数与蓝图不符：',length(tr))
    for(i in seq_len(nrow(b))) {
      c<-rp62_word_cell(tabs[[k]],b$row[i],b$col[i])
      x<-rp13_word_text(c$cell)
      # 数值样本只检查位置和数值型，而不把样本当作最终数据。
      if(isTRUE(b$demo_numeric[i])) {
        if(!grepl('^[[:space:]]*[-−+]?[0-9]|^[[:space:]]*<0[.]',x))
          stop('S',k,' R',b$row[i],' C',b$col[i],': 示例数值格发生结构性变化')
      } else if(!identical(trimws(gsub('[[:space:]]+',' ',x)),
                                trimws(gsub('[[:space:]]+',' ',b$text[i])))) {
        stop('S',k,' R',b$row[i],' C',b$col[i],': 静态文字与v6.2蓝图不同，请更新模板/蓝图')
      }
    }
  }
  tabs
}

rp62_validate_cell_coverage <- function(cells,blueprint) {
  rp_assert(cells,c('table','row','col','value','source_key'),'R生成的Word动态单元格')
  x<-copy(as.data.table(cells))
  x[,`:=`(table=as.integer(table),row=as.integer(row),col=as.integer(col))]
  if(anyNA(x[,.(table,row,col,value)]) || anyDuplicated(x,by=c('table','row','col')))
    stop('Word填值存在空字段或重复定位（禁止以最后一次覆盖处理）')
  b<-blueprint[,.(table,row,col,demo_numeric,text)]
  q<-merge(x,b,by=c('table','row','col'),all.x=TRUE)
  if(anyNA(q$demo_numeric))stop('R有目标单元格不在v6.2 Word蓝图')
  bad4<-q[table==4L & !(row%in%1:25 & col==2L)]
  if(nrow(bad4))stop('Table 4只允许修改第3列的25项实际插补方法')
  method_conflict<-q[table==4L & value!=text]
  if(nrow(method_conflict))stop('MISSING_SOURCE: Table 4实际插补方法与v6.2 Word不一致；首个为行',
    method_conflict$row[1],'：Word=',method_conflict$text[1],
    '；上游=',method_conflict$value[1], '。请核实插补方法后更新Word及蓝图。')
  bad<-q[!demo_numeric & !(text %in% c('—','')) & table!=4L]
  if(nrow(bad))stop('不允许覆盖非数值静态文字单元格，首个：',
    paste(bad[1,.(table,row,col)],collapse='/'))
  missing<-b[demo_numeric==TRUE][!x,on=.(table,row,col)]
  if(nrow(missing))stop('MISSING_SOURCE: v6.2 Word仍有样表数值未被R覆盖，首个：S',
      missing$table[1],' / R',missing$row[1],' / C',missing$col[1],
      '；总数=',nrow(missing))
  # 非适用单元格可写—；已展示数值的单元格不得以—替代来源。
  invalid<-q[demo_numeric & value=='—']
  if(nrow(invalid))stop('MISSING_SOURCE: 数值格被填为—：',
    paste(invalid[1,.(table,row,col)],collapse='/'))
  x
}

rp62_word_write <- function(template, output, cells, blueprint, overwrite=FALSE) {
  if(!requireNamespace('xml2',quietly=TRUE)||!requireNamespace('zip',quietly=TRUE))
    stop('生成v6.2 DOCX需要R包 xml2、zip')
  rp_require(template,'v6.2 Word模板')
  rp_guard(output)
  tmp<-tempfile('m01_v62_docx_');dir.create(tmp)
  on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
  utils::unzip(template,exdir=tmp)
  docpath<-file.path(tmp,'word','document.xml')
  if(!file.exists(docpath))stop('DOCX包缺少word/document.xml')
  doc<-xml2::read_xml(docpath)
  tabs<-rp62_template_check(doc,blueprint)
  cells<-rp62_validate_cell_coverage(cells,blueprint)
  for(i in seq_len(nrow(cells))) {
    z<-cells[i]
    tc<-rp62_word_cell(tabs[[z$table]],z$row,z$col)
    # 只改动态值的<w:t>文本；不动w:tcPr、w:tblPr、页码或OMML。
    rp13_word_put(tc$cell,z$value)
  }
  xml2::write_xml(doc,docpath,encoding='UTF-8')
  rp_dir(dirname(output))
  ziptmp<-tempfile('m01_v62_packed_',tmpdir=dirname(output),fileext='.docx')
  on.exit(unlink(ziptmp),add=TRUE)
  files<-list.files(tmp,recursive=TRUE,all.files=TRUE,include.dirs=FALSE,no..=TRUE)
  zip::zipr(ziptmp,files=files,root=tmp,compression_level=6L,mode='mirror')
  entries<-utils::unzip(ziptmp,list=TRUE)$Name
  needed<-c('[Content_Types].xml','_rels/.rels','word/document.xml',
            'word/_rels/document.xml.rels')
  if(!all(needed%in%entries))stop('输出DOCX丢失OOXML文件')
  if(!file.copy(ziptmp,output,overwrite=overwrite))stop('无法写出文件：',output)
  rp_msg('v6.2 Word生成：',output)
  invisible(output)
}

rp62_provenance <- function(cells, paths,style) {
  d<-copy(as.data.table(cells))
  d[,`:=`(table=as.integer(table),row=as.integer(row),col=as.integer(col))]
  if(!'status'%in%names(d))d[,status:='READY']
  src_names<-vapply(d$table,function(t)
    paste(intersect(rp62_source_candidates(t),names(paths)),collapse=';'),character(1))
  d[,candidate_sources:=src_names]
  d[,candidate_paths:=vapply(strsplit(src_names,';',fixed=TRUE),function(keys)
      paste(unlist(paths[keys[nzchar(keys)]]),collapse=';'),character(1))]
  hashes<-vapply(paths,rp_md5,character(1))
  d[,candidate_md5:=vapply(strsplit(src_names,';',fixed=TRUE),function(keys)
      paste(hashes[keys[nzchar(keys)]],collapse=';'),character(1))]
  d[,source_code_version:='M01_13_v6.2']
  d[,display_strategy_order:=paste(c('S0','SN','ST','SNT'),collapse='>')]
  setcolorder(d,c('table','row','col','value','source_key','status',
    'candidate_sources','candidate_paths','candidate_md5','source_code_version',
    'display_strategy_order'))
  d
}

# Invoke the OOXML writer through the same stage-13 entry point.
rp63_generate <- function(project_root, code_dir, analysis_dir, overwrite) {
  template_dir <- Sys.getenv("M01_REPORT_TEMPLATE_DIR", file.path(
    dirname(project_root), "03_研究方案", "02 研究方案与预期结果", "00_方案更新",
    "Risk modifiers修改方案_v2.0_20261102", "示例表", "附件表格"))
  template <- Sys.getenv("M01_SUPPLEMENTARY_TEMPLATE", file.path(
    template_dir, "Supplementary_Tables_v6.3_20261116.docx"))
  lock <- Sys.getenv("M01_SUPPLEMENTARY_STYLE_LOCK", file.path(
    template_dir, "Supplementary_Tables_样式锁定清单_v6.3_20261116.json"))
  legacy <- Sys.getenv("M01_REPORT_LEGACY_OUTPUT", file.path(
    analysis_dir, "publication_outputs_M01", "v6.2"))
  report <- Sys.getenv("M01_REPORT_OUTPUT", file.path(
    analysis_dir, "publication_outputs_M01", "v6.3"))
  output <- file.path(report, "tables", "M01_Supplementary_Tables_v6.3.docx")
  options(M01_REPORT_OVERWRITE = overwrite)
  rp_guard(output)
  rp_guard(sub("\\.docx$", ".cells.csv", output))
  python <- rp63_python()
  args <- c(file.path(code_dir, "report_v63_docx.py"), template, lock, analysis_dir,
            file.path(project_root, "LJW_outcome", "data_intervention"), legacy,
            file.path(code_dir, "config", "report_v62", "M01_supplementary_v6.2_blueprint.csv"),
            output)
  table2_template <- Sys.getenv("M01_MAIN_TABLE2_TEMPLATE",
    file.path(dirname(template_dir), "Table_2_External_cohort.docx"))
  rp_require(table2_template, "Main Table 2 template")
  rp_guard(file.path(dirname(output),"M01_Main_Table_2_v6.3.docx"))
  code <- system2(python, shQuote(c(args, "--table2-template", table2_template)))
  if (code != 0L) stop("Supplementary table export failed (exit ", code, ")")
  invisible(output)
}




rp63_python <- function() {
  configured <- Sys.getenv("M01_REPORT_PYTHON")
  candidates <- if (nzchar(configured)) configured else unique(c(
    file.path(Sys.getenv("USERPROFILE"), "anaconda3", "python.exe"),
    unname(Sys.which(c("python", "python3")))))
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  for (python in candidates) {
    if (grepl("[/\\\\]WindowsApps[/\\\\]python[0-9]*\\.exe$", python, ignore.case = TRUE)) next
    result <- tryCatch(suppressWarnings(system2(
      python, c("-c", shQuote("import lxml.etree; print('M01_PYTHON_READY')")),
      stdout = TRUE, stderr = TRUE)), error = function(e) character())
    if ("M01_PYTHON_READY" %in% result && is.null(attr(result, "status"))) return(python)
  }
  stop("Python with lxml was not found. Set M01_REPORT_PYTHON to the full path of python.exe.",
       if (nzchar(configured)) paste0(" Configured path: ", configured), call. = FALSE)
}


