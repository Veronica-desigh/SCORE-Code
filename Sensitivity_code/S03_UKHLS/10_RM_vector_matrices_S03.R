# 10 | Publication-style vector PDF matrices (no risks, no confidence intervals).
# Each PDF has one page per RM set x strategy; 2019-ST pairing is explicitly 2019 S0.
# =========================== USER CONFIG ===========================
DATASET_ID <- "S03_UKHLS"   # Fixed UKHLS only
THREADS <- 5L               # reserved for consistency; drawing is sequential and light
TEST_N <- 0L
SMOKING_MODE <- "same_wave"
PROJECT_ROOT <- ""
CELL_DIGITS <- 1L
# ===================================================================
DATASET_ID <- "S03_UKHLS"
THREADS<-as.integer(Sys.getenv("S03_RM_THREADS",as.character(THREADS)))
TEST_N<-as.integer(Sys.getenv("S03_RM_TEST_N",Sys.getenv("S03_ANALYSIS_TEST_N",as.character(TEST_N))))
SMOKING_MODE<-Sys.getenv("S03_SMOKING_MODE",SMOKING_MODE)
PROJECT_ROOT<-Sys.getenv("S03_RM_PROJECT_ROOT",Sys.getenv("S03_PROJECT_ROOT",PROJECT_ROOT))
if(!is.finite(THREADS)||THREADS<1)stop("Invalid settings")
a<-grep("^--file=",commandArgs(),value=TRUE);fr<-Filter(Negate(is.null),lapply(sys.frames(),function(z)z$ofile))
start<-if(length(fr))dirname(tail(fr,1L)[[1L]])else if(length(a))dirname(sub("^--file=","",a[1L])) else getwd()
root<-if(nzchar(PROJECT_ROOT))normalizePath(PROJECT_ROOT,winslash="/",mustWork=TRUE)else normalizePath(start,winslash="/",mustWork=TRUE)
if(!nzchar(PROJECT_ROOT))while(!dir.exists(file.path(root,"LJW_outcome"))){up<-dirname(root);if(up==root)stop("Project root missing");root<-up}
out<-file.path(root,"LJW_outcome","data_intervention",DATASET_ID)
if(TEST_N>0L)out<-file.path(out,paste0("test_",TEST_N))
if(SMOKING_MODE=="wave2_prior")out<-file.path(out,"wave2_smoking_sensitivity")
out<-file.path(out,"external_RM_analyses")
figdir<-file.path(out,"04_vector_figures");dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
logdir<-file.path(out,"logs");dir.create(logdir,recursive=TRUE,showWarnings=FALSE)
t0<-Sys.time();logfile<-file.path(logdir,"10_RM_figures.log")
log<-function(...){msg<-sprintf("[%s +%.1fs] %s",format(Sys.time(),"%H:%M:%S"),as.numeric(difftime(Sys.time(),t0,units="secs")),paste0(...,collapse=""));message(msg);cat(msg,"\n",file=logfile,append=TRUE)}
read_pooled<-function(name){path<-file.path(out,"03_pooled",name);if(!file.exists(path))stop("Run 09 first: ",path);read.csv(path,fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE,check.names=FALSE)}
risk_path<-file.path(out,"03_pooled","04_risk_matrix_MI_mean.csv")
latest_08<-file.path(out,"quality_control","08_completed_chains.csv")
if(file.exists(latest_08)&&file.exists(risk_path)&&
   file.info(latest_08)$mtime>file.info(risk_path)$mtime)
  stop("Step 08 output is newer than 09 pooled matrices; rerun step 09 before drawing")
risk<-read_pooled("04_risk_matrix_MI_mean.csv")
groups<-read_pooled("05_six_group_matrix_MI_mean.csv")
risk_levels<-c("L_score","M_score","H_rm","H_score","VH_score")
transition_levels<-c("unchanged_untreated","newly_intervention","intervention_removed","treatment_intensified","unchanged_intensity","treatment_deintensified")
risk_labels19<-c("L (SCORE)","M (SCORE)","H (RM upgrade)","H (SCORE)","VH (SCORE)")
risk_labels25<-c("L (SCORE2)","M (SCORE2)","H (RM upgrade)","H (SCORE2)","VH (SCORE2)")
group_labels<-c("UU","NI","IR","TI","UI","TD")
group_descriptions<-c("UU: unchanged untreated","NI: newly intervention","IR: intervention removed",
                       "TI: treatment intensified","UI: unchanged intensity","TD: treatment deintensified")
format_n<-function(n){if(abs(n-round(n))<1e-8)format(round(n),big.mark=",",scientific=FALSE,trim=TRUE)else formatC(n,format="f",digits=CELL_DIGITS,big.mark=",")}
format_pct<-function(p){if(is.na(p))return("NA");if(p>0&&p<0.05)return("<0.1%");sprintf("%.1f%%",p)}
blues<-grDevices::colorRampPalette(c("#f8fafc","#dce9f6","#7badd4","#356f9f"))(101)
plot_heatmap<-function(mat,nr,nc,row_names,col_names,title_txt,subtitle,footnote,axisnames){
  if(nrow(mat)!=nr*nc)stop("Matrix is not exhaustive: ",title_txt)
  graphics::plot.new()
  graphics::plot.window(xlim=c(-2.2,nc+0.30),ylim=c(-1.50,nr+1.05),xaxs="i",yaxs="i",asp=1)
  maxp<-max(mat$pct_total,na.rm=TRUE)
  for(i in seq_len(nr))for(j in seq_len(nc)){
    cell<-mat[mat$row_id==i & mat$col_id==j,,drop=FALSE]
    if(nrow(cell)!=1L||!is.finite(cell$n)||cell$n<0||!is.finite(cell$pct_total))stop("Invalid cell")
    col<-if(cell$n==0)"#F5F5F7"else blues[pmax(1L,pmin(101L,1L+round(100*cell$pct_total/max(1e-8,maxp))))]
    graphics::rect(j-1,nr-i,j,nr-i+1,col=col,border="#CDD4DE",lwd=0.8)
    if(cell$n>0){
      graphics::text(j-0.5,nr-i+0.57,labels=format_n(cell$n),cex=if(nc==6)0.76 else 0.88,font=2,col="#18283a")
      graphics::text(j-0.5,nr-i+0.30,labels=paste0("(",format_pct(cell$pct_total),")"),cex=if(nc==6)0.67 else 0.74,col="#27394a")
    }else graphics::text(j-0.5,nr-i+0.48,labels="-",cex=1,col="#9ca3af")
  }
  for(i in seq_len(nr))graphics::text(-0.18,nr-i+0.5,row_names[i],adj=1,cex=if(nr==6)0.76 else 0.85,col="#29384b")
  for(j in seq_len(nc))graphics::text(j-0.5,nr+0.17,col_names[j],cex=if(nc==6)0.80 else 0.83,col="#29384b",font=2)
  graphics::text(-2.1,nr/2,axisnames[1],srt=90,cex=0.95,font=2,col="#2b4055")
  graphics::text(nc/2,nr+0.75,axisnames[2],cex=0.95,font=2,col="#2b4055")
  # Emphasize boundary between below H and H/VH; matches reference figure quadrant logic.
  if(nr==5&&nc==5){graphics::segments(2,0,2,nr,lwd=2,col="#536477");graphics::segments(0,nr-2,nc,nr-2,lwd=2,col="#536477")}
  graphics::text(nc/2,-0.33,subtitle,cex=0.79,col="#344054")
  graphics::text(nc/2,-0.80,footnote,cex=0.64,col="#667085")
  graphics::title(main=title_txt,cex.main=1.15,font.main=2,line=0.5)
  graphics::box()
}

make_matrix<-function(data,rows,cols,values){
  n<-length(rows);m<-length(cols)
  data$row_id<-match(data[[values[1]]],rows)
  data$col_id<-match(data[[values[2]]],cols)
  if(anyNA(data$row_id)||anyNA(data$col_id)||anyDuplicated(data[c("row_id","col_id")]))stop("Matrix category mismatch")
  data
}

pairs<-unique(risk[c("RM_set","strategy")]);pairs<-pairs[order(match(pairs$RM_set,c("strict","expanded")),match(pairs$strategy,c("S0","SN","ST","SNT"))),,drop=FALSE]
log("Risk matrix figure: ",nrow(pairs)," pages")
grDevices::pdf(file.path(figdir,"Figure_RM_risk_reclassification.pdf"),width=12.8,height=9.3,onefile=TRUE,useDingbats=FALSE,family="sans")
tryCatch({
  for(i in seq_len(nrow(pairs))){
    set<-pairs$RM_set[i];sc<-pairs$strategy[i]
    d<-risk[risk$RM_set==set & risk$strategy==sc,,drop=FALSE]
    d<-make_matrix(d,risk_levels,risk_levels,c("origin_2019","origin_2025"))
    sub<-if(sc=="ST")"2019 S0 versus 2025 ST (2019 ST not applicable)"else paste0("2019 ",if(sc=="SNT")"SN (same as SNT)"else sc," versus 2025 ",sc)
    plot_heatmap(d,5,5,risk_labels19,risk_labels25,
      paste0(DATASET_ID," | ",toupper(set)," RM | ",sc," | Risk reclassification"),sub,
      "Cell: MI-mean estimated number (percentage of all participants); no ASCVD outcome estimates",
      c("2019 risk categories","2025 risk categories"))
    log("Risk matrix page ",i,"/",nrow(pairs),": ",set,"/",sc)
  }
},finally=grDevices::dev.off())

pairs<-unique(groups[c("RM_set","strategy")]);pairs<-pairs[pairs$strategy!="S0",,drop=FALSE]
pairs<-pairs[order(match(pairs$RM_set,c("strict","expanded")),match(pairs$strategy,c("SN","ST","SNT"))),,drop=FALSE]
log("Six-group matrix figure: ",nrow(pairs)," pages")
grDevices::pdf(file.path(figdir,"Figure_RM_six_group_transitions.pdf"),width=12.8,height=9.6,onefile=TRUE,useDingbats=FALSE,family="sans")
tryCatch({
  for(i in seq_len(nrow(pairs))){
    set<-pairs$RM_set[i];sc<-pairs$strategy[i]
    d<-groups[groups$RM_set==set & groups$strategy==sc,,drop=FALSE]
    d<-make_matrix(d,transition_levels,transition_levels,c("score_based_group","RM_based_group"))
    sub<-if(sc=="ST")"Baseline S0 six groups -> 2019 S0 vs 2025 ST six groups"else paste0("Baseline S0 six groups -> ",sc," six groups")
    plot_heatmap(d,6,6,group_labels,group_labels,
      paste0(DATASET_ID," | ",toupper(set)," RM | ",sc," | Treatment-decision transitions"),sub,
      paste(group_descriptions,collapse="  |  "),c("Score-based (S0) group","RM-based group"))
    log("Treatment matrix page ",i,"/",nrow(pairs),": ",set,"/",sc)
  }
},finally=grDevices::dev.off())
log("10 completed; vector PDF results in ",figdir)
