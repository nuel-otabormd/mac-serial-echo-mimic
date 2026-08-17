# 07_figures.R: Figures 1-4 from the results objects (PNG, TIFF and PDF at 600 dpi) into outputs/figures/.
suppressPackageStartupMessages({library(ggplot2); library(gridExtra); library(grid); library(cmprsk)})
r1 <- readRDS("outputs/results/results_part1.rds"); r2 <- readRDS("outputs/results/results_part2.rds"); d <- readRDS("data/frame.rds")
OUT <- "outputs/figures/"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE); pdf(NULL)   # null device: gtable building never opens Rplots.pdf
save_fig <- function(g, name, w, h){ ggsave(paste0(OUT,name,".png"), g, width=w, height=h, dpi=600, bg="white"); ggsave(paste0(OUT,name,".tiff"), g, width=w, height=h, dpi=600, bg="white", compression="lzw"); ggsave(paste0(OUT,name,".pdf"), g, width=w, height=h, bg="white") }

# ---------- Figure 1: cohort assembly ----------
fl <- r2$flow; m3 <- unique(r1$main[r1$main$domain=="D3", c("stratum","def","n_pt","events")]); ev <- function(st, de) format(m3$events[m3$stratum==st & m3$def==de], big.mark=",")
cb <- unique(r1$confirmed_baseline[, c("def","n_pt","events")])
box <- function(x, y, w, h, lab, fill="#EEF4FB", fs=9.5) { grid.roundrect(x=unit(x,"npc"), y=unit(y,"npc"), width=unit(w,"npc"), height=unit(h,"npc"), gp=gpar(fill=fill, col="#1F4E79", lwd=1), r=unit(2,"mm")); grid.text(lab, x=unit(x,"npc"), y=unit(y,"npc"), gp=gpar(fontfamily="sans", fontsize=fs, lineheight=0.95)) }
seg <- function(x0,y0,x1,y1,arr=TRUE) grid.lines(x=unit(c(x0,x1),"npc"), y=unit(c(y0,y1),"npc"), arrow=if (arr) arrow(length=unit(2,"mm"), type="closed") else NULL, gp=gpar(fill="black", col="black", lwd=1))
draw_flow <- function(){ grid.newpage(); pushViewport(viewport(width=0.94, height=0.96))
  # top box
  bx <- 0.40; box(bx,0.925,0.58,0.10, sprintf("Adults with two or more transthoracic echocardiography\nepisodes carrying the annular calcification field\nn = %s", format(fl["eligible"],big.mark=",")))
  seg(bx,0.875,bx,0.745)                                   # down to analysis cohort (arrowhead at box top)
  seg(bx,0.81,0.585,0.81)                                  # branch to exclusion box (arrowhead at its left edge)
  box(0.79,0.81,0.41,0.11, sprintf("Excluded, n = %d\nMitral prosthesis or annuloplasty ring\nat or before index, n = %d\nRheumatic mitral valve on any report, n = %d", fl["excl_pros"]+fl["excl_rheum"], fl["excl_pros"], fl["excl_rheum"]), fill="#F2F2F2", fs=9)
  box(bx,0.695,0.40,0.09, sprintf("Analysis cohort\nn = %s", format(fl["analysed"],big.mark=",")))
  seg(bx,0.65,bx,0.585,arr=FALSE); seg(0.14,0.585,0.78,0.585,arr=FALSE)
  for (x in c(0.14,0.46,0.78)) seg(x,0.585,x,0.535)
  box(0.14,0.44,0.27,0.18, sprintf("No reported calcification at index\n(onset cohort)\nn = %s\n\nFirst-observed events %s\nConfirmed events %s", format(fl["blank"],big.mark=","), ev("onset","first"), ev("onset","confirmed")), fs=9)
  box(0.46,0.44,0.27,0.18, sprintf("Mild calcification at index\n(progression cohort)\nn = %s\n\nFirst-observed events %s\nConfirmed events %s", format(fl["mild"],big.mark=","), ev("progression","first"), ev("progression","confirmed")), fs=9)
  box(0.78,0.44,0.34,0.18, sprintf("Moderate or severe calcification at index\n(advanced disease)\nn = %s\n\nDescribed together with patients reaching\nmoderate or greater during follow-up", format(fl["modsev"],big.mark=",")), fs=9)
  seg(0.14,0.35,0.14,0.30); box(0.14,0.21,0.27,0.17, sprintf("Stricter onset cohort:\nfirst two episodes both blank,\nfollowed from the second\nn = %s\nevents: %s first-observed,\n%s confirmed", format(cb$n_pt[cb$def=="first"],big.mark=","), format(cb$events[cb$def=="first"],big.mark=","), format(cb$events[cb$def=="confirmed"],big.mark=",")), fill="#F2F2F2", fs=9)
  popViewport() }
png(paste0(OUT,"Figure1_flow.png"), width=7.5, height=7.5, units="in", res=600, bg="white", family="sans"); draw_flow(); dev.off()
tiff(paste0(OUT,"Figure1_flow.tiff"), width=7.5, height=7.5, units="in", res=600, compression="lzw", bg="white", family="sans"); draw_flow(); dev.off()
pdf(paste0(OUT,"Figure1_flow.pdf"), width=7.5, height=7.5, family="sans"); draw_flow(); dev.off(); cat("Figure 1 written\n")

# ---------- Figure 2: reliability of the reported grade ----------
th2 <- theme_classic(base_size=13, base_family="sans") + theme(axis.text=element_text(colour="black"), legend.position="bottom", strip.background=element_blank(),
        strip.text=element_text(size=12), plot.title=element_text(face="bold", size=17, hjust=0), plot.title.position="plot", legend.text=element_text(size=12))
E <- r2$hmm$emission; L <- r2$hmm$emission_L; U <- r2$hmm$emission_U; lab <- c("None","Mild","Moderate","Severe")
em <- expand.grid(true=lab, reported=lab); em$p <- as.vector(E[1:4,1:4]); em$lo <- as.vector(L[1:4,1:4]); em$hi <- as.vector(U[1:4,1:4]); em$true <- factor(em$true, levels=rev(lab)); em$reported <- factor(em$reported, levels=lab)
em$lab1 <- sprintf("%.0f%%", 100*em$p); em$lab2 <- sprintf("(%.0f\u2013%.0f)", 100*em$lo, 100*em$hi); em$diag <- as.character(em$true)==as.character(em$reported); em$col <- ifelse(em$p>0.5,"white","black")
gA <- ggplot(em, aes(reported, true, fill=p)) + geom_tile(colour="white", linewidth=0.8) +
  geom_text(aes(label=lab1, fontface=ifelse(diag,"bold","plain")), colour=em$col, size=5.2, family="sans", nudge_y=0.15) + geom_text(aes(label=lab2), colour=em$col, size=4.2, family="sans", nudge_y=-0.15) +
  scale_fill_gradient(low="#F7FBFF", high="#08519C", limits=c(0,1), guide="none") + scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  labs(x="Grade reported on the echocardiogram", y="Model-estimated underlying grade", title="A") + th2 + theme(axis.line=element_blank(), axis.ticks=element_blank()) + coord_fixed(ratio=0.55)
fa <- r2$first_read_fate; fa <- fa[fa$group=="age_cat",]; fa$level <- factor(fa$level, levels=c("<65","65-74","75-84","85+")); fa$xlab <- paste0(as.character(fa$level), "\nn = ", fa$n)
fa$stratum <- factor(fa$stratum, levels=c("onset","progression"), labels=c("No reported calcification at index","Mild calcification at index"))
fl2 <- rbind(data.frame(fa[,c("stratum","level","xlab")], fate="Confirmed", v=fa$confirmed), data.frame(fa[,c("stratum","level","xlab")], fate="Refuted", v=fa$refuted), data.frame(fa[,c("stratum","level","xlab")], fate="Never re-examined", v=fa$last_echo))
fl2$fate <- factor(fl2$fate, levels=c("Never re-examined","Refuted","Confirmed")); fl2$xlab <- factor(fl2$xlab, levels=unique(fl2$xlab[order(fl2$stratum, fl2$level)]))
gB <- ggplot(fl2, aes(xlab, v, fill=fate)) + geom_col(width=0.62) + facet_wrap(~stratum, scales="free_x") +
  scale_fill_manual(values=c(`Never re-examined`="#D9D9D9", Refuted="#E69F00", Confirmed="#0072B2"), breaks=c("Confirmed","Refuted","Never re-examined"), name=NULL) +
  scale_y_continuous(labels=function(x) paste0(100*x,"%"), expand=c(0,0), limits=c(0,1)) +
  labs(x="Age at index, years", y="Share of first moderate or\ngreater reads", title="B") + th2 + theme(legend.key.size=unit(0.45,"cm"))
g2 <- arrangeGrob(gA, gB, ncol=1, heights=c(1.05,1))
save_fig(g2, "Figure2_reliability", 8, 9.5)
cat("Figure 2 written\n")
# ---------- Figure 3: hazard ratios (A) and ratio of hazard ratios (B) ----------
labS <- c(age10="Age, per 10 years", male="Male sex", av_catmild="AVC, mild vs none", av_catmod="AVC, moderate vs none", av_catsev="AVC, severe vs none",
          ef_catlt40="LVEF <40% vs \u226550%", ef_cat40to49="LVEF 40\u201349% vs \u226550%", zE="E/e', per SD", zla="LA dimension, per SD", zivs="Septal thickness, per SD",
          zbmi="BMI, per SD", af="AF or flutter", eg30to60="eGFR 30\u201359 vs \u226560", eglt30="eGFR <30 vs \u226560", esrd="Dialysis dependence")
th3 <- theme_classic(base_size=13, base_family="sans") + theme(axis.text=element_text(colour="black"), legend.position="bottom", strip.background=element_blank(),
        strip.text=element_text(size=12), plot.title=element_text(face="bold", size=17, hjust=0), plot.title.position="plot", legend.text=element_text(size=12))
m <- r1$main; d3 <- m[m$domain=="D3" & m$term %in% names(labS),]; d3$term <- factor(d3$term, levels=rev(names(labS)))
d3$stratum <- factor(d3$stratum, levels=c("onset","progression"), labels=c("Onset: no reported calcification at index","Progression: mild calcification at index"))
d3$def <- factor(d3$def, levels=c("first","confirmed"), labels=c("First-observed","Confirmed"))
gA <- ggplot(d3, aes(x=hr, y=term, shape=def, colour=def)) + geom_vline(xintercept=1, linetype="dashed", colour="grey50") +
  geom_errorbar(aes(xmin=lo, xmax=hi), width=0, position=position_dodge(width=0.6), linewidth=0.6, orientation="y") + geom_point(position=position_dodge(width=0.6), size=2.6, fill="white") +
  facet_wrap(~stratum) + scale_x_log10(breaks=c(0.25,0.5,1,2,4,8,16), labels=c("0.25","0.5","1","2","4","8","16")) + scale_y_discrete(labels=function(x) labS[x]) +
  scale_shape_manual(values=c(16,21), name=NULL) + scale_colour_manual(values=c("#0072B2","#D55E00"), name=NULL) + labs(x="Hazard ratio (95% CI)", y=NULL, title="A") + th3
it <- r1$interaction; it$term2 <- sub("^prog:","", sub(":prog$","", it$term))
labB <- c(age10="Age, per 10 years", male="Male sex", av_lin="AVC, per grade", ef_catlt40="LVEF <40% vs \u226550%", ef_cat40to49="LVEF 40\u201349% vs \u226550%", zE="E/e', per SD",
          zla="LA dimension, per SD", zivs="Septal thickness, per SD", zbmi="BMI, per SD", af="AF or flutter", eg30to60="eGFR 30\u201359 vs \u226560", eglt30="eGFR <30 vs \u226560", esrd="Dialysis dependence")
b <- it[it$def=="first" & it$term2 %in% names(labB),]; b$term2 <- factor(b$term2, levels=rev(names(labB))); b$sig <- ifelse(b$lo > 1 | b$hi < 1, "excludes 1", "includes 1")
gB <- ggplot(b, aes(x=hr, y=term2, colour=sig)) + geom_vline(xintercept=1, linetype="dashed", colour="grey50") +
  geom_errorbar(aes(xmin=lo, xmax=hi), width=0, linewidth=0.7, orientation="y") + geom_point(size=2.6) +
  annotate("text", x=0.98, y=length(labB)+0.85, label="Weaker for progression", hjust=1, size=4.2, family="sans", fontface="italic", colour="grey40") +
  annotate("text", x=1.02, y=length(labB)+0.85, label="Stronger for progression", hjust=0, size=4.2, family="sans", fontface="italic", colour="grey40") +
  scale_x_log10(breaks=c(0.4,0.5,0.7,1,1.5,2), labels=c("0.4","0.5","0.7","1","1.5","2"), limits=c(0.3,2.1), oob=scales::oob_keep) + scale_y_discrete(labels=function(x) labB[x]) +
  scale_colour_manual(values=c(`excludes 1`="#B2182B", `includes 1`="grey45"), guide="none") + coord_cartesian(clip="off", ylim=c(0.6, length(labB)+1.1)) +
  labs(x="Ratio of hazard ratios, progression relative to onset (95% CI)", y=NULL, title="B") + th3
g <- arrangeGrob(gA, gB, ncol=1, heights=c(1.15,1))
save_fig(g, "Figure3_forest", 8, 10.5)
cat("Figure 3 written\n")
# ---------- Figure 4: cumulative incidence of calcific mitral stenosis ----------
th <- theme_classic(base_size=12, base_family="sans") + theme(axis.text=element_text(colour="black"), plot.title=element_text(face="bold", size=11), legend.position="bottom")
x <- d[d$ms0==0 & !is.na(d$ms0),]; grp <- ifelse(x$sev1==0,"No reported calcification", ifelse(x$sev1==1,"Mild","Moderate or severe"))
ev <- ifelse(!is.na(x$t_ms),1, ifelse(!is.na(x$t_death) & (is.na(x$t_last) | x$t_death<=x$t_last+30),2,0)); tt <- pmax(ifelse(!is.na(x$t_ms), x$t_ms, ifelse(ev==2, x$t_death, x$t_last)),1)/365.25
ci <- cuminc(tt, ev, grp); df <- do.call(rbind, lapply(names(ci)[grepl(" 1$", names(ci))], function(n) data.frame(group=sub(" 1$","",n), time=ci[[n]]$time, est=ci[[n]]$est, var=ci[[n]]$var)))
df$group <- factor(df$group, levels=c("No reported calcification","Mild","Moderate or severe")); df <- df[df$time<=8,]
g4 <- ggplot(df, aes(time, est, colour=group, linetype=group)) + geom_step(linewidth=0.9) + scale_colour_manual(values=c("#56B4E9","#E69F00","#D55E00"), name="Annular calcification at index") + scale_linetype_manual(values=c("dotted","dashed","solid"), name="Annular calcification at index") +
  scale_y_continuous(labels=function(v) paste0(100*v,"%")) + scale_x_continuous(breaks=0:8) + labs(x="Years from index echocardiogram", y="Cumulative incidence of calcific mitral stenosis\n(mean gradient 5 mmHg or more, death as competing risk)") + th
nr <- sapply(levels(df$group), function(g) sapply(0:8, function(t) sum(tt[grp==g] >= t))); rownames(nr) <- 0:8
save_fig(g4, "Figure4_calcific_stenosis_CIF", 7, 5.2); write.csv(nr, "outputs/tables/Figure4_numbers_at_risk.csv")
cat("Figure 4 regenerated at 600 dpi\n")
invisible(dev.off())
