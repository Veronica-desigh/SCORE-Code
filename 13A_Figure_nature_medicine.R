# M01 main Figures 3 and 4: unified journal-style publication artwork.
# Based on 13A_Figure_nature_medicine(4).R. All original Figure 3 inputs,
# fixed-S0 definitions, and risk/LDL-C estimands are preserved.
# New artwork: five integrated panels, no large section-title banners;
# honest bar lengths with data-adjacent labels; consistent CI symbols.
# Figure 3: preserve the targeted A/B/C axis corrections from the latest supplied script.
# Figure 4: preserve the previously agreed v6.11 table–forest design.
# Export both figures with stable filenames (no version suffix).
# Run after M01 14 has exported the v6.3 figure-data CSVs.
required <- c("ggplot2", "patchwork", "svglite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Install R packages: ", paste(missing, collapse=", "))
library(ggplot2)
library(patchwork)

find_project_root <- function() {
  configured <- Sys.getenv("M01_PROJECT_ROOT", "")
  if (nzchar(configured)) return(normalizePath(configured, winslash="/", mustWork=TRUE))
  # Support source() in RStudio as well as Rscript from any working directory.
  script_paths <- unlist(lapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) character() else as.character(frame$ofile)
  }), use.names=FALSE)
  cli_paths <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))
  starts <- unique(c(getwd(), dirname(c(script_paths, cli_paths))))
  for (start in starts) {
    candidate <- normalizePath(start, winslash="/", mustWork=FALSE)
    repeat {
      if (dir.exists(file.path(candidate, "LJW_code")) &&
          dir.exists(file.path(candidate, "LJW_outcome"))) return(candidate)
      parent <- dirname(candidate)
      if (identical(parent, candidate)) break
      candidate <- parent
    }
  }
  stop("Cannot locate project root. Set M01_PROJECT_ROOT to the 04_Code directory.")
}
root <- find_project_root()
analysis <- Sys.getenv("M01_ANALYSIS_OUTPUT", file.path(root, "LJW_outcome",
                                                        "data_intervention", "M01_ASCVD", "analysis_results"))
report <- file.path(analysis,"publication_outputs_M01","v6.3")
input <- file.path(report,"figure_data")
output <- file.path(report,"figures","nature_medicine")
message("M01 figure-data directory: ", input)
font <- Sys.getenv("M01_NATURE_FONT","Arial")
palette <- c(S0="#767676",SN="#2166AC",ST="#E66101",SNT="#1B7837")
shapes <- c(S0=1,SN=1,ST=1,SNT=1) # Uniform open circles
stroke <- 0.15                # axes/zero reference (mm in ggplot linewidth)
ci_stroke <- 0.30             # stronger confidence-interval line (ggplot mm)
ci_point_size <- 1.60         # hollow point diameter (mm)
ci_point_stroke <- 0.20       # hollow point border (ggplot stroke)
bar_outline <- 0.106         # black bar outline: about 0.30 pt in ggplot mm
ci_cap_fig4 <- 0.11          # short caps in Figure 4 table–forest
ci_cap_fig3 <- 0.085         # short caps, data units in horizontal dodged plots
# One visual system across Figures 3 and 4; all lengths in ggplot units unless stated.
bar_dodge_width <- 0.78
bar_width <- 0.68
fig3_row_ratio <- c(0.49, 0.51) # Taller A/B row; B decision groups are less compressed.
fig3_size_mm <- c(width=240, height=208) # 195-mm artwork + 13-mm legend
fig4_size_mm <- c(width=190, height=210)
groups <- setNames(c("Unchanged untreated","New intervention","Intervention removed",
                     "Treatment intensified","Unchanged intensity","Treatment de-intensified"),
                   c("unchanged_untreated","newly_intervention","intervention_removed",
                     "treatment_intensified","unchanged_intensity","treatment_deintensified"))
# Style lock v1.1: sizes are points at the final 180 mm artwork width.
type <- list(base=8.0, label=7.7, axis=8.1, title=9.2, tag=9.5, body=7.3, small=7.0)
# Figure 4 uses its own text-size scale, independent of Figure 3.
type4 <- list(body=7.6, subgroup=8.0, dimension=8.6, header=8.4,
              axis=8.0, tick=7.7, direction=7.0)
theme_nm <- theme_classic(base_size=type$base,base_family=font) +
  theme(axis.line=element_line(linewidth=stroke),axis.ticks=element_line(linewidth=stroke),
        axis.ticks.length=grid::unit(1,"mm"),axis.text=element_text(size=type$label,colour="black"),
        axis.title=element_text(size=type$axis),plot.title=element_text(size=type$title,face="bold",hjust=0),
        axis.title.y = element_text(
          size = type$axis,
          margin = margin(r = 2)
        ),
        plot.tag=element_text(size=type$tag,face="bold"),
        plot.tag.position=c(0,1),plot.margin=margin(3,3,3,1),
        strip.background=element_blank(),strip.text=element_text(size=type$label),
        legend.text=element_text(size=type$label),legend.key.size=grid::unit(3,"mm"),
        legend.spacing.x=grid::unit(2,"mm"),panel.spacing=grid::unit(2.0,"mm"))
read_panel <- function(figure, name, n) {
  path <- file.path(input,sprintf("Figure_%02d_%s.csv",figure,name))
  if (!file.exists(path)) stop("Missing M01 14 output: ",path)
  x <- read.csv(path,fileEncoding="UTF-8-BOM")
  stopifnot(nrow(x)==n,all(x$scenario %in% names(palette)))
  x$scenario <- factor(x$scenario,levels=names(palette))
  x
}
group_factor <- function(x) factor(groups[x],levels=rev(unname(groups)))
a <- read_panel(3,"classification",6L)
b <- read_panel(3,"reclassification",18L)
c <- read_panel(3,"ldl_change",30L)
# Fixed S0 groups compare the two frameworks. Plot the RM-induced change
# in that between-framework contrast, not two separate within-year panels.
# Linearity permits subtraction of the paired MI means (same people and MIs).
c19 <- c[c$guideline==2019,]
c25 <- c[c$guideline==2025,]
c_key <- function(x) paste(x$fixed_group, x$scenario, sep="|")
stopifnot(!anyDuplicated(c_key(c19)), !anyDuplicated(c_key(c25)),
          setequal(c_key(c19), c_key(c25)))
c19 <- c19[match(c_key(c25), c_key(c19)),]
stopifnot(all(c19$n==c25$n), all(c19$m==c25$m))
c <- c25
c$mean_delta_mmol_l <- c25$mean_delta_mmol_l-c19$mean_delta_mmol_l
c$guideline <- NULL
c$contrast <- "[(D2025-D2019)RM-(D2025-D2019)S0]"
d <- read_panel(3,"overall_rd",4L)
e <- read_panel(3,"fixed_group_rd",16L)
b$group <- group_factor(b$row_group)
c$group <- group_factor(c$fixed_group)
e$group <- group_factor(e$group)
d$group <- factor(d$scenario,rev(names(palette)))
# Keep the six/five/four S0 decision-group definitions unchanged.
# For horizontal grouped bars, a reversed factor is necessary to display
# SN, ST, SNT from TOP to BOTTOM when grouped by ggplot's dodge layout.
pd <- position_dodge(width=bar_dodge_width)
fill <- scale_fill_manual(values=palette,name=NULL,guide="none")
colour <- scale_colour_manual(values=palette,name=NULL,breaks=names(palette),drop=FALSE)
shape <- scale_shape_manual(values=shapes,name=NULL,breaks=names(palette),drop=FALSE)
zero <- geom_vline(xintercept=0,linewidth=stroke,colour="grey62",linetype=2)

# Consistent label construction and true-data positioning: never add an arbitrary
# strategy-specific x-offset (e.g., +7 pp / +0.16 mmol/L).
b$dodge_group <- factor(as.character(b$scenario),levels=c("SNT","ST","SN"))
c$dodge_group <- factor(as.character(c$scenario),levels=c("SNT","ST","SN"))
# Format honest bar-end labels. Near-zero bars get a readable label immediately
# to the right/left of x=0; this changes ONLY the label location, not bar length.
# No strategy-specific x offset and no nonlinear axis transformation.
# (diagonal N - original row N) / original row N * 100.
# Original membership stays fixed; incoming participants are excluded.
stopifnot(all(is.finite(b$pct_moved)),all(b$pct_moved>=0 & b$pct_moved<=100))
b$retention_change_pp <- -b$pct_moved
b$bar_label <- ifelse(b$pct_moved>0,
                      ifelse(b$pct_moved<1,sprintf("%.2f",b$retention_change_pp),
                             sprintf("%.1f",b$retention_change_pp)),"0")
b$label_x <- pmin(b$retention_change_pp,-1.25)-1.05
c$bar_label <- ifelse(is.finite(c$mean_delta_mmol_l),
                      ifelse(abs(c$mean_delta_mmol_l) < 0.02,
                             sprintf("%.3f", c$mean_delta_mmol_l),
                             sprintf("%.2f", c$mean_delta_mmol_l)), NA_character_)
c$label_x <- ifelse(c$mean_delta_mmol_l < 0,
                    pmin(c$mean_delta_mmol_l,-0.025)-0.012,
                    pmax(c$mean_delta_mmol_l,0.025)+0.014)
c$label_justify <- ifelse(c$mean_delta_mmol_l < 0, 1, 0)
# Bar labels are positioned using the SAME dodge object as the bars, so their
# y coordinates stay aligned with each relevant treatment strategy.

panel_title <- function() theme(
  plot.title=element_text(size=type$title,face="bold",hjust=0,
                          margin=margin(b=3)),
  plot.title.position="plot",
  plot.tag=element_blank())

# PANEL A. Preserve real original/upgraded percentages and both guidelines.
# The subtitle uses a quiet light/dark key; bar outlines remain fine black.
# A single plotting panel gives both guideline blocks one continuous baseline.
a$axis_x <- match(as.character(a$scenario),c("SN","ST","SNT")) +
  3.4*(match(as.character(a$guideline),c("2019","2025"))-1)
panel_a <- ggplot(a, aes(axis_x,fill=scenario)) +
  geom_col(aes(y=original_pct),width=.66,colour="black",linewidth=bar_outline) +
  geom_tile(aes(y=upgraded_center,height=upgraded_pct),width=.66,
            alpha=.38,colour=NA) +
  geom_tile(aes(y=upgraded_center,height=upgraded_pct),width=.66,
            fill=NA,colour="black",linewidth=bar_outline) +
  geom_text(aes(y=final_pct,label=sprintf("%.1f%%",final_pct)),
            vjust=-.45,size=type$label/ggplot2::.pt,family=font) +
  scale_x_continuous(breaks=c(1:3,1:3+3.4),
                     labels=rep(c("SN","ST","SNT"),2),
                     expand=expansion(add=.6)) +
  annotate("text",x=c(2,5.4),y=-9,label=c("2019","2025"),
           family=font,size=type$label/ggplot2::.pt) + fill +
  scale_y_continuous(breaks=c(0,20,40,60,80),expand=expansion(mult=c(0,.01))) +
  coord_cartesian(ylim=c(0,88),expand=c(bottom=FALSE),clip="off") +
  geom_hline(yintercept=0,colour="black",linewidth=.25,inherit.aes=FALSE) +
  labs(title="A   High/very-high risk classification",
       x=NULL,y="Participants (%)") + theme_nm + panel_title() +
  theme(legend.position="none",plot.margin=margin(3,3,17,1),
        plot.subtitle=element_blank(),
        strip.text=element_text(size=type$label,face="plain"),
        axis.title.y = element_text(
          size = type$axis,
          margin = margin(r = 2)
        ),
        strip.placement="outside",strip.switch.pad.grid=grid::unit(1.5,"mm"))

# PANEL B. Six genuine fixed-S0 groups with three coloured grouped bars each.
# Both bars and labels use identical dodge. Positive bar starts exactly at x=0.
panel_b <- ggplot(b,aes(retention_change_pp,group,fill=scenario,group=dodge_group)) +
  geom_col(position=pd,width=bar_width,colour="black",
           linewidth=bar_outline,show.legend=FALSE) +
  geom_text(aes(x=label_x,label=bar_label),position=pd,
            hjust=1,size=type$small/ggplot2::.pt,family=font,
            colour="grey25",show.legend=FALSE) +
  scale_x_continuous(breaks=c(-80,-60,-40,-20,0),expand=expansion(mult=c(0,0))) +
  scale_y_discrete(expand=expansion(add=c(.5,.5))) +
  # Keep the discrete y-axis half-category expansion: the x-axis baseline
  # and x=0 vertical axis then meet at the lower-left panel corner.
  # x remains zero-expanded by scale_x_continuous() above.
  coord_cartesian(xlim=c(min(-86,min(b$label_x)-6),0),
                  clip="off") +
  fill +
  zero +
  labs(title="B   Treatment-decision reclassification under RM strategies",
       x="Change in group membership (pp)",y=NULL) +
  theme_nm+panel_title() + theme(legend.position="none",
                                 # Use native x/y axes, rather than a geom_vline with a different endpoint.
                                 axis.line.y=element_line(linewidth=stroke,colour="black"),
                                 axis.line.x=element_line(linewidth=stroke,colour="black"),
                                 axis.text.y=element_text(size=type$label),
                                 plot.margin=margin(3,5,3,2))

# PANEL C. Keep the *existing* paired two-guideline contrast of RM-minus-S0
# change in required minimum LDL-C reduction; never paste mock-up numbers.
panel_c <- ggplot(c,aes(mean_delta_mmol_l,group,fill=scenario,group=dodge_group)) +
  geom_col(position=pd,width=bar_width,colour="black",
           linewidth=bar_outline,show.legend=FALSE) +
  geom_text(aes(x=label_x,label=bar_label,hjust=label_justify),
            position=pd,size=type$small/ggplot2::.pt,family=font,
            colour="grey25",show.legend=FALSE) +
  # Bar origin is the true zero; no additional x=0 line is needed for the
  # current nonnegative input. Keep the reference for any future negatives.
  {if (any(c$mean_delta_mmol_l < 0, na.rm=TRUE)) zero} +
  scale_x_continuous(expand=expansion(mult=c(0,0))) +
  # Start the visible x range at exactly zero when all changes are
  # nonnegative. Preserve the negative range if the real data contain it.
  coord_cartesian(xlim=c(min(0,c$mean_delta_mmol_l,c$label_x,na.rm=TRUE),
                         max(0,c$mean_delta_mmol_l,c$label_x,na.rm=TRUE)+.11),
                  clip="off") +
  fill +
  labs(title="C   Changes in mean required minimum LDL-C reduction",
       x="Change in 2025 − 2019 difference (mmol/L)",y=NULL) +
  theme_nm + panel_title() + theme(
    legend.position = "none",
    axis.text.y = element_text(size = type$label),
    axis.title.x = element_text(
      size = type$axis,
      margin = margin(t = 8)
    ),
    plot.margin = margin(3, 4, 4, 2)
  )

forest <- function(x,tag,title) {
  stopifnot(all(is.finite(x$effect_pp)),
            all(x$low_pp<=x$effect_pp & x$effect_pp<=x$high_pp))
  breaks <- if (tag=="D") c(-.05,0,.05,.10) else waiver()
  x_limits <- if (tag=="D") c(min(-.07,min(x$low_pp)-.012),max(.125,max(x$high_pp)+.012)) else
    range(c(x$low_pp,x$high_pp))+c(-.08,.08)
  ggplot(x,aes(effect_pp,group,colour=scenario,shape=scenario,group=scenario)) + zero +
    geom_errorbar(aes(xmin=low_pp,xmax=high_pp),orientation="y",
                  width=ci_cap_fig3,linewidth=ci_stroke,position=pd) +
    geom_point(size=ci_point_size,stroke=ci_point_stroke,position=pd) +
    colour+shape+scale_x_continuous(breaks=breaks) +
    coord_cartesian(xlim=x_limits,clip="off")+
    labs(title=paste0(tag,"   ",title),x="2025 − 2019 RD (pp; 95% CI)",y=NULL) +
    theme_nm+panel_title()+theme(
      axis.text.y=element_text(size=type$label),
      axis.title.x=element_text(size=type$axis, margin=margin(t=8))
    )
}

# Five panels arranged as approved: A/B on the taller first row, C/D/E below.
# NO separator line or section banners. Uniform title size in every panel.
# Figure 3 legends share one custom, single-line bottom strip. There is no
# guide in Panel A or an automatic guide collected from the forest panels.
# Important: the left-hand axis title of A was being pushed far away from
# its ticks by patchwork aligning A to the long category labels in C/E.
# Free ONLY A's left edge from that cross-row alignment; keep A/B dimensions.
if (!"free" %in% getNamespaceExports("patchwork")) {
  stop("Figure 3 v6.15 requires patchwork::free() (update the patchwork package).")
}
panel_a_layout <- patchwork::free(panel_a, side="l", type="panel")
top <- wrap_plots(panel_a_layout,panel_b,widths=c(.40,.60))
bottom <- wrap_plots(
  panel_c,
  forest(d,"D","Overall RD")+theme(legend.position="none"),
  forest(e,"E","RD between 2025 and 2019 under RM strategies") +
    theme(legend.position="none"),
  widths=c(.45,.19,.36))
fig3_body <- wrap_plots(top,bottom,ncol=1,heights=fig3_row_ratio)

# ONE bottom legend row with TWO clearly separated titled groups.
# Panel A classification components uses TWO groups of THREE compact swatches:
#   dark SN/ST/SNT = Score-based; light SN/ST/SNT = RM-upgraded.
# The second group shows the four colored open circles for RM strategies.
# Positions are specified along a dedicated 240-mm legend strip. Each item has
# a reserved, nonoverlapping horizontal region; NOTHING wraps to a second row.
# No title/text keys appear within Panel A itself.
light_component <- function(hex, alpha=.38) {
  channels <- as.numeric(grDevices::col2rgb(hex)[,1L])
  mixed <- round(alpha * channels + (1-alpha) * 255)
  grDevices::rgb(mixed[1L],mixed[2L],mixed[3L],maxColorValue=255)
}
# Shift the entire one-row Figure 3 legend towards the left.
# Units correspond to the legend strip's x scale (0–240).
# Keep the first title inside the drawing area to prevent clipping.
legend_shift_x <- -9
legend_left_safe_x <- 1
legend_x <- function(x) max(legend_left_safe_x, x + legend_shift_x)
legend_label <- function(x, label, bold=FALSE, size=type$small) {
  x <- legend_x(x)
  annotate("text",x=x,y=.5,label=label,hjust=0,vjust=.5,
           family=font,fontface=if (bold) "bold" else "plain",
           size=size/ggplot2::.pt,colour="black")
}
# Each three-swatch key fully represents SN (blue), ST (orange), SNT (green).
# Swatches are adjacent but separated by a small white gap; every swatch has
# the same fine black border as the corresponding Panel A bars.
legend_triplet <- function(x, light=FALSE) {
  x <- legend_x(x)
  strategies <- c("SN","ST","SNT")
  lapply(seq_along(strategies), function(i) {
    strategy <- strategies[i]
    color <- unname(palette[strategy])
    if (light) color <- light_component(color)
    left <- x + (i-1)*3.1
    annotate("rect",xmin=left,xmax=left+2.65,
             ymin=.385,ymax=.615,fill=color,
             colour="black",linewidth=bar_outline)
  })
}
legend_circle <- function(x, strategy) {
  x <- legend_x(x)
  annotate("point",x=x,y=.5,shape=21,size=1.65,
           stroke=ci_point_stroke,fill="white",
           colour=unname(palette[strategy]))
}
fig3_legend <- ggplot() +
  # Legend group 1 (left): classification source, each with the FULL palette.
  legend_label(3,"Panel A classification components",bold=TRUE) +
  legend_triplet(64,light=FALSE) +
  legend_label(75,"Score-based") +
  legend_triplet(97,light=TRUE) +
  legend_label(108,"RM-upgraded") +
  # Keep >= 12 mm of whitespace before the second titled legend group.
  legend_label(150,"RM strategy",bold=TRUE) +
  legend_circle(181,"S0") + legend_label(185,"S0") +
  legend_circle(196,"SN") + legend_label(200,"SN") +
  legend_circle(211,"ST") + legend_label(215,"ST") +
  legend_circle(226,"SNT") + legend_label(230,"SNT") +
  scale_x_continuous(limits=c(0,240),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,1),expand=c(0,0)) +
  coord_cartesian(clip="off") +
  theme_void(base_family=font) +
  theme(plot.margin=margin(0,0,0,0))
# Dedicated legend height does not squeeze the A-E plot body.
fig3 <- wrap_plots(fig3_body,fig3_legend,ncol=1,heights=c(195,13))


# FIGURE 4. Five-column subgroup table–forest, with pale heading bands.
f <- read_panel(4,"subgroup_risk_differences",32L)
stopifnot(!anyDuplicated(f[c("scope","scenario")]),length(unique(f$scope))==8L,
          all(f$rd_lo<=f$rd_est & f$rd_est<=f$rd_hi))
f <- f[order(-f$y),]
scope_order <- unique(f$scope)
dimension_order <- c("Age","Sex","SES","Hypertension")
# Use ONE shared physical row map for the five aligned columns. A three-tier
# hierarchy needs distinct vertical gaps: dimension -> subgroup -> strategy.
# Do not independently space the text and forest panels.
dimension_names <- c(Age="Age",Sex="Sex",SES="SES",Hypertension="Hypertension")
strategy_order <- c("S0", "SN", "ST", "SNT")
layout_rows <- vector("list", length(scope_order) * length(strategy_order))
subgroup_rows <- vector("list", length(scope_order))
dimension_rows <- vector("list", length(dimension_order))
section_rules <- numeric()
y_cursor <- 100
row_i <- 0L
sub_i <- 0L
for (dim_i in seq_along(dimension_order)) {
  dim_key <- dimension_order[dim_i]
  dimension_rows[[dim_i]] <- data.frame(label=unname(dimension_names[dim_key]), y=y_cursor)
  y_cursor <- y_cursor - 2.10  # space below the dimension heading
  dim_scopes <- scope_order[vapply(scope_order, function(sc)
    identical(as.character(f$dimension[match(sc,f$scope)]), dim_key), logical(1))]
  stopifnot(length(dim_scopes)==2L)
  for (j in seq_along(dim_scopes)) {
    sc <- dim_scopes[j]
    sub_i <- sub_i + 1L
    subgroup_rows[[sub_i]] <- data.frame(scope=sc, y=y_cursor)
    y_cursor <- y_cursor - 1.50  # room between the subgroup heading and S0
    for (scenario in strategy_order) {
      row_i <- row_i + 1L
      layout_rows[[row_i]] <- data.frame(scope=sc,scenario=scenario,y=y_cursor)
      y_cursor <- y_cursor - 1.08  # regular spacing within S0/SN/ST/SNT
    }
    if (j < length(dim_scopes)) {
      y_cursor <- y_cursor - 1.40  # additional gap between subgroup levels
    } else if (dim_i < length(dimension_order)) {
      section_rules <- c(section_rules, y_cursor - 1.02)
      y_cursor <- y_cursor - 2.35  # clear gap between top-level dimensions
    }
  }
}
row_map <- do.call(rbind,layout_rows)
stopifnot(nrow(row_map)==nrow(f), !anyDuplicated(row_map[c("scope","scenario")]))
f$y <- row_map$y[match(paste(f$scope,as.character(f$scenario)),
                       paste(row_map$scope,row_map$scenario))]
stopifnot(all(is.finite(f$y)))
headers <- merge(do.call(rbind,subgroup_rows),
                 f[f$scenario=="S0",c("scope","label","n")],by="scope",sort=FALSE)
headers <- headers[match(scope_order,headers$scope),]
headers$n_label <- formatC(headers$n,format="f",digits=0,big.mark=",")
headers$label <- sub(" · S0$","",headers$label)
dimensions <- do.call(rbind,dimension_rows)
# Publication-style dimension bands, not alternating shading for data rows.
# Use the SAME four y locations in all five aligned subplots. Bands are
# limited to headings and never cover subgroup or RM strategy values.
heading_fill <- "#EDF2F6"
heading_ink <- "#263A50"
heading_halfheight <- 0.72  # vertical row-map units; 2.10 to next subgroup
heading_bands <- transform(dimensions,
                           ymin=y-heading_halfheight,
                           ymax=y+heading_halfheight)
# Background must be drawn BEFORE all separators, text, zero line and CIs.
heading_background <- geom_rect(data=heading_bands,
                                aes(ymin=ymin,ymax=ymax),xmin=-Inf,xmax=Inf,
                                inherit.aes=FALSE,fill=heading_fill,colour=NA)
f$label <- as.character(f$scenario) # indentation is controlled by draw_x
f$n_label <- ""
# Conditional on each completed dataset's subgroup weights. Between-MI
# variance includes any changes in those weights across imputations.
mi <- read.csv(file.path(analysis,"strategy_risks_and_effects","per_imputation",
                         "all_scenario_risks_and_effects.csv"),fileEncoding="UTF-8-BOM")
mi <- subset(mi,time==15 & estimand=="rd_2025_vs_2019")
mi$scenario <- sub("^SHT$","SNT",sub("^SH$","SN",mi$scenario))
totals <- subset(mi,scope=="overall")
stopifnot(!anyDuplicated(totals[c("scenario","imputation_id")]))
denominator <- setNames(totals$n_standardized,paste(totals$scenario,totals$imputation_id))
contribution_ci <- function(scope_key,strategy) {
  d <- mi[mi$scope==scope_key & mi$scenario==strategy,]
  stopifnot(nrow(d)>1L,!anyDuplicated(d$imputation_id),
            setequal(d$imputation_id,totals$imputation_id[totals$scenario==strategy]))
  w <- 100*d$n_standardized/unname(denominator[paste(strategy,d$imputation_id)])
  q <- w*d$Q
  u <- w^2*d$U
  stopifnot(all(is.finite(q)),all(is.finite(u)&u>=0),all(w>0&w<=100))
  between <- (1+1/length(q))*var(q)
  variance <- mean(u)+between
  lambda <- if (variance>0) between/variance else 0
  complete_df <- min(d$dfcom)
  observed_df <- (complete_df+1)/(complete_df+3)*complete_df*(1-lambda)
  inverse_df <- lambda^2/(length(q)-1)
  if (is.finite(observed_df) && observed_df>0) inverse_df <- inverse_df+1/observed_df
  half <- qt(.975,df=1/inverse_df)*sqrt(variance)
  c(estimate=mean(q),lo=mean(q)-half,hi=mean(q)+half)
}
interval <- t(vapply(seq_len(nrow(f)),function(i)
  contribution_ci(f$scope[i],as.character(f$scenario[i])),numeric(3)))
stopifnot(max(abs(interval[,"estimate"]-f$contribution_pp))<1e-8)
format_interval <- function(est,lo,hi,digits) {
  num <- function(x,d) {
    value <- round(x,d)
    if(value==0) value <- 0
    sprintf(paste0("%.",d,"f"),value)
  }
  vapply(seq_along(est),function(i) {
    paste0(num(est[i],digits)," (",num(lo[i],digits),", ",num(hi[i],digits),")")
  },character(1))
}
for (prefix in c("rd"))
  f[[prefix]] <- format_interval(f[[paste0(prefix,"_est")]],f[[paste0(prefix,"_lo")]],
                                 f[[paste0(prefix,"_hi")]],2)
f$contribution <- format_interval(interval[,"estimate"],interval[,"lo"],interval[,"hi"],3)
limits <- c(min(f$y)-2.0,max(dimensions$y)+1.1)
rules <- section_rules
separator <- geom_hline(yintercept=rules,colour="grey85",linewidth=.15)
y_axis <- scale_y_continuous(limits=limits,expand=c(0,0),breaks=NULL)
# A table–forest hybrid: open plotting area, regular aligned rows, no shaded
# backgrounds or unnecessary y-axis frame. The same y coordinates are used in
# all five columns. Use wide space for the forest rather than its neighbors.
text_theme <- theme_nm + theme(
  axis.line=element_blank(),axis.ticks=element_blank(),
  axis.text.x=element_blank(),axis.title.x=element_blank(),
  panel.grid=element_blank(),
  plot.margin=margin(1,1.1,1,1.1),
  plot.title=element_text(size=type4$header,face="bold",hjust=0,margin=margin(b=2)))

column <- function(field,title) {
  is_left <- field %in% c("label","n_label")
  # Explicit x coordinate instead of depending on an in-aes conditional.
  plot_data <- f
  plot_data$draw_x <- if (field=="label") 0.12 else if (is_left) 0 else 1
  g <- ggplot(plot_data,aes(draw_x,y,label=.data[[field]])) +
    heading_background + separator +
    geom_text(hjust=if (is_left) 0 else 1,size=type4$body/ggplot2::.pt,
              family=font) + y_axis +
    scale_x_continuous(limits=c(0,1),expand=c(0,0)) +
    coord_cartesian(clip="off") + labs(title=title,x=NULL,y=NULL) + text_theme
  if (is_left) {
    g <- g + geom_text(data=transform(headers,draw_x=if (field=="label") 0.045 else 0),aes(label=.data[[field]]),
                       hjust=0,size=type4$subgroup/ggplot2::.pt,family=font,fontface="bold")
  }
  # Header alignment mirrors the numeric column: the two effect-value
  # headings are flush right; subgroup and N headings remain flush left.
  if (field %in% c("rd", "contribution")) {
    g <- g + theme(plot.title=element_text(size=type4$header,
                                           face="bold",hjust=1,margin=margin(b=2)))
  }
  if (field=="label") {
    g <- g + geom_text(data=transform(dimensions,draw_x=0),
                       aes(x=draw_x,y=y,label=label),inherit.aes=FALSE,hjust=0,
                       size=type4$dimension/ggplot2::.pt,fontface="bold",
                       colour=heading_ink,family=font)
  }
  g
}

# Four restrained full-width heading bands provide hierarchy; white subgroup
# and strategy rows remain untouched. Draw the band BEFORE the zero line/CIs.
forest4 <- ggplot(f,aes(rd_est,y,colour=scenario,shape=scenario)) +
  heading_background + separator + zero +
  geom_errorbar(aes(xmin=rd_lo,xmax=rd_hi),orientation="y",
                width=ci_cap_fig4,linewidth=ci_stroke) +
  geom_point(size=ci_point_size,stroke=ci_point_stroke) + y_axis + colour + shape +
  scale_x_continuous(breaks=c(-.1,0,.1,.2,.3),expand=expansion(mult=c(.02,.03))) +
  coord_cartesian(xlim=c(-.16,.38),clip="off")+
  labs(title="",x="2025 − 2019 RD (pp)",y=NULL)+
  theme_nm+theme(legend.position="none",
                 axis.text.x=element_text(size=type4$tick),
                 axis.title.x=element_text(size=type4$axis),
                 axis.line.y=element_blank(),
                 axis.ticks.y=element_blank(),axis.text.y=element_blank(),
                 plot.margin=margin(1,2,7,2),
                 plot.title=element_text(size=type4$header,face="bold",hjust=0,margin=margin(b=2)))
# Put the directional wording at x locations, not inside a padded x-axis title.
forest4 <- forest4 + annotate("text",x=-.075,y=min(limits)+.55,
                              label="Favors 2025",size=type4$direction/ggplot2::.pt,
                              family=font,hjust=.5,colour="grey35") +
  annotate("text",x=.22,y=min(limits)+.55,
           label="Favors 2019",size=type4$direction/ggplot2::.pt,
           family=font,hjust=.5,colour="grey35")
# Figure 4: five visual columns. Remove the 2025−NC and 2019−NC RD
# columns entirely, while retaining the 2025−2019 effect text and plot,
# both anchored on the same 32 row coordinates as the contribution column.
# Give the forest panel more usable horizontal space (+10 relative width units)
# without changing the total 190 mm canvas. Keep single-line headers.
fig4 <- wrap_plots(
  column("label","Subgroup"),
  column("n_label","N"),
  column("rd","2025 − 2019 RD (pp; 95% CI)"),
  forest4,
  column("contribution","Contribution (pp; 95% CI)"),
  nrow=1,widths=c(28,12,44,65,41))


# Stable, version-free exports for both figures.
dir.create(output,recursive=TRUE,showWarnings=FALSE)
figure_specs <- list(
  list(number=3L, plot=fig3, width=unname(fig3_size_mm["width"]), height=unname(fig3_size_mm["height"])),
  list(number=4L, plot=fig4, width=unname(fig4_size_mm["width"]), height=unname(fig4_size_mm["height"]))
)
qa_dir <- Sys.getenv("M01_NATURE_QA_DIR", "")
if (nzchar(qa_dir)) dir.create(qa_dir,recursive=TRUE,showWarnings=FALSE)
for (spec in figure_specs) {
  prefix <- file.path(output,sprintf("Figure_%02d_Nature_Medicine",spec$number))
  ggsave(paste0(prefix,".pdf"),spec$plot,device=cairo_pdf,
         width=spec$width,height=spec$height,units="mm",bg="white")
  ggsave(paste0(prefix,".svg"),spec$plot,device=svglite::svglite,
         width=spec$width,height=spec$height,units="mm",bg="white")
  if (nzchar(qa_dir)) {
    ggsave(file.path(qa_dir,sprintf("Figure_%02d.png",spec$number)),
           spec$plot,width=spec$width,height=spec$height,units="mm",
           dpi=180,bg="white")
  }
  message("Figure ",spec$number," written to ",prefix,".[pdf/svg]")
}
