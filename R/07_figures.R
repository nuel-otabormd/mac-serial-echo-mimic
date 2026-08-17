# 07_figures.R: Figures 1-4 from the results objects (PNG, TIFF and PDF at 600 dpi) into outputs/figures/.
if (!grepl("UTF-8", Sys.getlocale("LC_CTYPE"), ignore.case=TRUE)) for (loc in c("en_US.UTF-8","C.UTF-8")) if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break   # write en dashes and other symbols as UTF-8 even when the calling shell uses the C locale
suppressPackageStartupMessages({library(ggplot2); library(gridExtra); library(grid); library(cmprsk)}); source("R/00_common.R")
r1 <- readRDS("outputs/results/results_part1.rds"); r2 <- readRDS("outputs/results/results_part2.rds"); d <- readRDS("data/frame.rds")
OUT <- "outputs/figures/"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE); pdf(NULL)   # null device: gtable building never opens Rplots.pdf
# PNG and TIFF (600 dpi, LZW) are drawn by ragg. The vector PDF needs a device that embeds Unicode glyphs (en dash, >=): quartz on macOS, cairo_pdf elsewhere;
# the base pdf() device would silently substitute a hyphen for the en dash. If no such device is available the PDF is skipped with a message and the run continues.
pdf_dev <- function(filename, width, height, bg="white", ...) {
  if (Sys.info()[["sysname"]] == "Darwin") grDevices::quartz(file=filename, type="pdf", width=width, height=height, bg=bg)
  else if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf(filename=filename, width=width, height=height, bg=bg)
  else stop("no Unicode-capable PDF device (cairo) on this system") }
save_fig <- function(g, name, w, h){
  ragg::agg_png(paste0(OUT,name,".png"), width=w, height=h, units="in", res=600, background="white"); grid.draw(g); invisible(dev.off())
  ragg::agg_tiff(paste0(OUT,name,".tiff"), width=w, height=h, units="in", res=600, background="white", compression="lzw"); grid.draw(g); invisible(dev.off())
  tryCatch({ pdf_dev(paste0(OUT,name,".pdf"), width=w, height=h); grid.draw(g); invisible(dev.off()) }, error=function(e) message("PDF for ", name, " skipped: ", conditionMessage(e))) }

# ---------- Figure 1: cohort assembly ----------
fl <- r2$flow; m3 <- unique(r1$main[r1$main$domain=="D3", c("stratum","def","n_pt","events")]); ev <- function(st, de) format(m3$events[m3$stratum==st & m3$def==de], big.mark=",")
cb <- unique(r1$confirmed_baseline[, c("def","n_pt","events")])
box <- function(x, y, w, h, lab, fill="#EEF4FB", fs=9.5) { grid.roundrect(x=unit(x,"npc"), y=unit(y,"npc"), width=unit(w,"npc"), height=unit(h,"npc"), gp=gpar(fill=fill, col="#1F4E79", lwd=1), r=unit(2,"mm")); grid.text(lab, x=unit(x,"npc"), y=unit(y,"npc"), gp=gpar(fontfamily="sans", fontsize=fs, lineheight=0.95)) }
seg <- function(x0,y0,x1,y1,arr=TRUE) grid.lines(x=unit(c(x0,x1),"npc"), y=unit(c(y0,y1),"npc"), arrow=if (arr) arrow(length=unit(2,"mm"), type="closed") else NULL, gp=gpar(fill="black", col="black", lwd=1))
draw_flow <- function(){ grid.newpage(); pushViewport(viewport(width=0.94, height=0.97))
  bx <- 0.34
  box(bx,0.955,0.50,0.075, sprintf("Adults with a transthoracic echocardiography study\ncarrying the annular calcification field, n = %s", format(fl["source_adults"],big.mark=",")), fs=9)
  seg(bx,0.9175,bx,0.865); seg(bx,0.895,0.60,0.895)
  box(0.80,0.895,0.40,0.05, sprintf("Only one echocardiography episode, n = %s", format(fl["single_episode"],big.mark=",")), fill="#F2F2F2", fs=9)
  box(bx,0.825,0.50,0.075, sprintf("Adults with two or more episodes\nn = %s", format(fl["eligible"],big.mark=",")), fs=9)
  seg(bx,0.7875,bx,0.705)                                  # down to analysis cohort
  seg(bx,0.755,0.60,0.755)                                 # branch to exclusion box
  box(0.80,0.755,0.40,0.11, sprintf("Excluded, n = %d\nMitral prosthesis or annuloplasty ring\nby the end of the index episode, n = %d\nRheumatic mitral valve reported by the\nend of the index episode, n = %d", fl["excl_pros"]+fl["excl_rheum"], fl["excl_pros"], fl["excl_rheum"]), fill="#F2F2F2", fs=8.5)
  box(bx,0.665,0.40,0.075, sprintf("Analysis cohort\nn = %s", format(fl["analysed"],big.mark=",")))
  seg(bx,0.6275,bx,0.585,arr=FALSE); seg(0.14,0.585,0.78,0.585,arr=FALSE)
  for (x in c(0.14,0.46,0.78)) seg(x,0.585,x,0.535)
  box(0.14,0.44,0.27,0.18, sprintf("No reported calcification at index\n(onset cohort)\nn = %s\n\nFirst-observed events %s\nConfirmed events %s", format(fl["blank"],big.mark=","), ev("onset","first"), ev("onset","confirmed")), fs=9)
  box(0.46,0.44,0.27,0.18, sprintf("Mild calcification at index\n(progression cohort)\nn = %s\n\nFirst-observed events %s\nConfirmed events %s", format(fl["mild"],big.mark=","), ev("progression","first"), ev("progression","confirmed")), fs=9)
  box(0.78,0.44,0.34,0.18, sprintf("Moderate or severe calcification at index\n(advanced disease)\nn = %s\n\nDescribed together with patients reaching\nmoderate or severe during follow-up", format(fl["modsev"],big.mark=",")), fs=9)
  seg(0.14,0.35,0.14,0.30); box(0.14,0.21,0.27,0.17, sprintf("Stricter onset cohort:\nfirst two episodes both blank,\nfollowed from the second\nn = %s\nevents: %s first-observed,\n%s confirmed", format(cb$n_pt[cb$def=="first"],big.mark=","), format(cb$events[cb$def=="first"],big.mark=","), format(cb$events[cb$def=="confirmed"],big.mark=",")), fill="#F2F2F2", fs=9)
  popViewport() }
ragg::agg_png(paste0(OUT,"Figure1_flow.png"), width=7.5, height=7.5, units="in", res=600, background="white"); draw_flow(); invisible(dev.off())
ragg::agg_tiff(paste0(OUT,"Figure1_flow.tiff"), width=7.5, height=7.5, units="in", res=600, background="white", compression="lzw"); draw_flow(); invisible(dev.off())
tryCatch({ pdf_dev(paste0(OUT,"Figure1_flow.pdf"), width=7.5, height=7.5); draw_flow(); invisible(dev.off()) }, error=function(e) message("PDF for Figure1_flow skipped: ", conditionMessage(e))); cat("Figure 1 written\n")

# ---------- Figure 2: reliability of the reported grade ----------
# Type sizes (points): tick 7, annotation 7.5, axis title 8.5, legend 6.8, panel letter 10 bold. Okabe-Ito bars; ColorBrewer Blues heatmap.
FS <- c(tick=7, annot=7.5, axlab=8.5, legend=6.8, letter=10); pt <- function(x) x/ggplot2::.pt
th_base <- theme_classic(base_size=FS["annot"], base_family="sans") +
  theme(axis.text=element_text(colour="black", size=FS["tick"]), axis.title=element_text(size=FS["axlab"]), strip.background=element_blank(), strip.text=element_text(size=FS["annot"], colour="black"),
        plot.title=element_text(face="bold", size=FS["letter"], hjust=0), plot.title.position="plot", legend.position="bottom", legend.text=element_text(size=FS["legend"]),
        legend.key.size=unit(0.28,"cm"), legend.margin=margin(0,0,0,0), axis.line=element_line(linewidth=0.3), axis.ticks=element_line(linewidth=0.3), plot.margin=margin(4,8,4,4))
E <- r2$hmm$emission; L <- r2$hmm$emission_L; U <- r2$hmm$emission_U; lab <- c("None","Mild","Moderate","Severe")
em <- expand.grid(true=lab, reported=lab); em$p <- as.vector(E[1:4,1:4]); em$lo <- as.vector(L[1:4,1:4]); em$hi <- as.vector(U[1:4,1:4]); em$true <- factor(em$true, levels=rev(lab)); em$reported <- factor(em$reported, levels=lab)
em$lab1 <- sprintf("%.0f%%", 100*em$p); em$lab2 <- sprintf("(%.0f\u2013%.0f)", 100*em$lo, 100*em$hi); em$diag <- as.character(em$true)==as.character(em$reported)
em$col <- ifelse(em$p >= 0.62, "white", "grey10")   # white text only on the darkest tiles keeps contrast above 4.5:1
blues <- c("#F7FBFF","#DEEBF7","#C6DBEF","#9ECAE1","#6BAED6","#4292C6","#2171B5","#08519C","#08306B")
gA <- ggplot(em, aes(reported, true, fill=p)) + geom_tile(colour=NA) +
  geom_text(aes(label=lab1, fontface=ifelse(diag,"bold","plain")), colour=em$col, size=pt(FS["axlab"]), family="sans", nudge_y=0.13) +
  geom_text(aes(label=lab2), colour=em$col, size=pt(FS["annot"]), family="sans", nudge_y=-0.17) +
  scale_fill_gradientn(colours=blues, limits=c(0,1), guide="none") + scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  labs(x="Grade reported on the echocardiogram", y="Model-estimated underlying grade", title="A") + th_base +
  theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_text(size=FS["tick"], colour="black"), axis.title.x=element_text(margin=margin(t=6)), axis.title.y=element_text(margin=margin(r=8)))
fa <- r2$first_read_fate; fa <- fa[fa$group=="age_cat",]; fa$level <- factor(fa$level, levels=c("<65","65-74","75-84","85+")); fa$xlab <- paste0(as.character(fa$level), "\nn = ", fa$n)
fa$stratum <- factor(fa$stratum, levels=c("onset","progression"), labels=c("No reported calcification at index","Mild calcification at index"))
fl2 <- rbind(data.frame(fa[,c("stratum","level","xlab")], fate="Confirmed", v=fa$confirmed), data.frame(fa[,c("stratum","level","xlab")], fate="Refuted", v=fa$refuted), data.frame(fa[,c("stratum","level","xlab")], fate="Never re-examined", v=fa$last_echo))
fl2$fate <- factor(fl2$fate, levels=c("Never re-examined","Refuted","Confirmed")); fl2$xlab <- factor(fl2$xlab, levels=unique(fl2$xlab[order(fl2$stratum, fl2$level)]))
gB <- ggplot(fl2, aes(xlab, v, fill=fate)) + geom_col(width=0.68) + facet_wrap(~stratum, scales="free_x") +
  scale_fill_manual(values=c(`Never re-examined`="#D9D9D9", Refuted="#E69F00", Confirmed="#0072B2"), breaks=c("Confirmed","Refuted","Never re-examined"), name=NULL) +
  scale_y_continuous(labels=function(x) paste0(100*x,"%"), expand=c(0,0), limits=c(0,1), breaks=c(0,0.25,0.5,0.75,1)) +
  labs(x="Age at the qualifying echocardiogram, years", y="Share of first moderate or\nsevere reads", title="B") + th_base +
  theme(panel.spacing.x=unit(0.35,"cm"), axis.title.x=element_text(margin=margin(t=4)), axis.title.y=element_text(margin=margin(r=6)), legend.spacing.x=unit(0.25,"cm"), legend.box.margin=margin(2,0,0,0))
g2 <- arrangeGrob(gA, gB, ncol=1, heights=c(1,1))
save_fig(g2, "Figure2_reliability", 7.1, 7.9)
cat("Figure 2 written\n")
# ---------- Figure 3: hazard ratios (A) and ratio of hazard ratios (B) ----------
labS <- c(age10="Age, per 10 years", male="Male sex", av_catmild="AVC, mild vs none", av_catmod="AVC, moderate vs none", av_catsev="AVC, severe vs none",
          ef_catlt40="LVEF <40% vs \u226550%", ef_cat40to49="LVEF 40\u201349% vs \u226550%", zE="E/e', per SD", zla="LA dimension, per SD", zivs="Septal thickness, per SD",
          zbmi="BMI, per SD", af="AF or flutter", eg30to60="eGFR 30\u201359 vs \u226560", eglt30="eGFR <30 vs \u226560", esrd="Dialysis dependence")
th_forest <- th_base + theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_text(size=FS["annot"], colour="black"), axis.text.x=element_text(size=FS["tick"]),
                             panel.spacing.x=unit(0.3,"cm"), legend.key.width=unit(0.5,"cm"), legend.spacing.x=unit(0.3,"cm"))
m <- r1$main; d3 <- m[m$domain=="D3" & m$term %in% names(labS),]; d3$y <- as.numeric(factor(d3$term, levels=rev(names(labS))))
d3$stratum <- factor(d3$stratum, levels=c("onset","progression"), labels=c("Onset: no reported calcification at index","Progression: mild calcification at index"))
d3$def <- factor(d3$def, levels=c("first","confirmed"), labels=c("First-observed","Confirmed")); d3$yy <- d3$y + ifelse(d3$def=="First-observed", 0.17, -0.17)   # first-observed drawn above confirmed
gA <- ggplot(d3, aes(x=hr, y=yy, colour=def, shape=def, fill=def)) + geom_vline(xintercept=1, linetype="42", colour="grey55", linewidth=0.35) +
  geom_errorbar(aes(xmin=lo, xmax=hi), width=0, linewidth=0.45, orientation="y", lineend="round") + geom_point(size=1.9, stroke=0.5) +
  facet_wrap(~stratum) + scale_x_log10(breaks=c(0.25,0.5,1,2,4,8,16), labels=c("0.25","0.5","1","2","4","8","16"), limits=c(0.18,30)) +
  scale_y_continuous(breaks=seq_along(labS), labels=rev(unname(labS)), expand=expansion(add=c(0.6,0.6))) +
  scale_shape_manual(values=c(`First-observed`=16, Confirmed=21), name=NULL) + scale_colour_manual(values=c(`First-observed`="#0072B2", Confirmed="#D55E00"), name=NULL) + scale_fill_manual(values=c(`First-observed`="#0072B2", Confirmed="white"), name=NULL) +
  labs(x="Hazard ratio (95% CI)", y=NULL, title="A") + th_forest
it <- r1$interaction; it$term2 <- sub("^prog:","", sub(":prog$","", it$term))
labB <- c(age10="Age, per 10 years", male="Male sex", av_lin="AVC, per grade", ef_catlt40="LVEF <40% vs \u226550%", ef_cat40to49="LVEF 40\u201349% vs \u226550%", zE="E/e', per SD",
          zla="LA dimension, per SD", zivs="Septal thickness, per SD", zbmi="BMI, per SD", af="AF or flutter", eg30to60="eGFR 30\u201359 vs \u226560", eglt30="eGFR <30 vs \u226560", esrd="Dialysis dependence")
b <- it[it$def=="first" & it$term2 %in% names(labB),]; b$y <- as.numeric(factor(b$term2, levels=rev(names(labB)))); b$sig <- ifelse(b$lo > 1 | b$hi < 1, "excludes 1", "includes 1")
gB <- ggplot(b, aes(x=hr, y=y, colour=sig)) + geom_vline(xintercept=1, linetype="42", colour="grey55", linewidth=0.35) +
  geom_errorbar(aes(xmin=lo, xmax=hi), width=0, linewidth=0.5, orientation="y", lineend="round") + geom_point(aes(fill=sig), shape=21, size=2.1, colour="white", stroke=0.4) +
  annotate("text", x=0.985, y=length(labB)+0.85, label="Weaker for progression", hjust=1, size=pt(FS["annot"]), family="sans", fontface="italic", colour="grey38") +
  annotate("text", x=1.015, y=length(labB)+0.85, label="Stronger for progression", hjust=0, size=pt(FS["annot"]), family="sans", fontface="italic", colour="grey38") +
  scale_x_log10(breaks=c(0.4,0.5,0.7,1,1.5,2), labels=c("0.4","0.5","0.7","1","1.5","2"), limits=c(0.3,2.1), oob=scales::oob_keep) +
  scale_y_continuous(breaks=seq_along(labB), labels=rev(unname(labB))) +
  scale_colour_manual(values=c(`excludes 1`="#B2182B", `includes 1`="grey42"), guide="none") + scale_fill_manual(values=c(`excludes 1`="#B2182B", `includes 1`="grey42"), guide="none") +
  coord_cartesian(clip="off", ylim=c(0.5, length(labB)+1.1)) + labs(x="Ratio of hazard ratios, progression relative to onset (95% CI)", y=NULL, title="B") + th_forest + theme(axis.title.x=element_text(margin=margin(t=6)))
g <- arrangeGrob(gA, gB, ncol=1, heights=c(1.12,1))
save_fig(g, "Figure3_forest", 7.1, 8.4)
cat("Figure 3 written\n")
# ---------- Figure 4: cumulative incidence of calcific mitral stenosis, with 95% bands and numbers at risk ----------
th <- theme_classic(base_size=11, base_family="sans") + theme(axis.text=element_text(colour="black"), plot.title=element_text(face="bold", size=11), legend.position="top", legend.text=element_text(size=9), legend.title=element_text(size=9))
x <- d[d$ms0==0 & !is.na(d$ms0),]; grp <- ifelse(x$sev1==0,"No reported calcification", ifelse(x$sev1==1,"Mild","Moderate or severe"))
z4 <- mk(x$t_ms, x$t_death, x$t_last_study); tt <- z4$tt; ev <- z4$ev
ci <- cuminc(tt, ev, grp); df <- do.call(rbind, lapply(names(ci)[grepl(" 1$", names(ci))], function(n) data.frame(group=sub(" 1$","",n), time=ci[[n]]$time, est=ci[[n]]$est, var=ci[[n]]$var)))
df$group <- factor(df$group, levels=c("No reported calcification","Mild","Moderate or severe")); df <- df[df$time<=8,]; df$lo <- pmax(df$est-1.96*sqrt(df$var),0); df$hi <- pmin(df$est+1.96*sqrt(df$var),1)
cols <- c("#56B4E9","#E69F00","#D55E00")
g4 <- ggplot(df, aes(time, est, colour=group, fill=group, linetype=group)) + geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.15, colour=NA, stat="identity") + geom_step(linewidth=0.9) +
  scale_colour_manual(values=cols, name="Annular calcification at index") + scale_fill_manual(values=cols, name="Annular calcification at index") + scale_linetype_manual(values=c("dotted","dashed","solid"), name="Annular calcification at index") +
  scale_y_continuous(labels=function(v) paste0(100*v,"%"), limits=c(0,NA)) + scale_x_continuous(breaks=0:8, limits=c(0,8)) + labs(x="Years from index echocardiogram", y="Cumulative incidence of calcific mitral stenosis\n(mean gradient 5 mmHg or more; death as competing risk)") + th
nr <- sapply(levels(df$group), function(g) sapply(0:8, function(t) sum(tt[grp==g] >= t))); rownames(nr) <- 0:8
rt <- data.frame(group=factor(rep(colnames(nr), each=nrow(nr)), levels=rev(levels(df$group))), time=rep(0:8, times=ncol(nr)), n=as.vector(nr))
gR <- ggplot(rt, aes(time, group, label=format(n, big.mark=","))) + geom_text(size=2.7, family="sans") + scale_x_continuous(breaks=0:8, limits=c(0,8)) + labs(x=NULL, y=NULL, title="Number at risk") +
  theme_classic(base_size=9, base_family="sans") + theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(colour="black", size=8), plot.title=element_text(size=9, face="plain"), plot.margin=margin(0,8,4,4))
g4r <- arrangeGrob(g4, gR, ncol=1, heights=c(4.2,1))
save_fig(g4r, "Figure4_calcific_stenosis_CIF", 7, 6.2); write.csv(nr, "outputs/tables/Figure4_numbers_at_risk.csv")
cat("Figure 4 regenerated at 600 dpi\n")
invisible(dev.off())
