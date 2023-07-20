library(nortest)
library(extrafont)
font_import(prompt = FALSE, pattern = "calibri")
loadfonts(device = "win")
library(ggplot2)
library(viridis)
library(tidyverse)
library(RColorBrewer)
library(sysfonts)
library(showtext)
font_add("Work Sans", regular = "~/Downloads/WorkSans-Regular.ttf")
showtext_auto()
#install.packages("reshape2")
library(reshape2)
#install.packages("cowplot")
library(cowplot)
library(patchwork)
#library(hopach)
library(extrafont)
#font_import(prompt = FALSE, pattern = "calibri")
#loadfonts(device = "win")


################################################################################
##RSVA analytical sensitivity QC CV graphs
################################################################################


df1 = read.table('RSVA_analytical.txt', 
                 sep='\t', 
                 header = T, 
                 colClasses = "numeric")
View(df1)
df1 = data.frame(df1)
View(df1)
df1.m = melt(df1, 
             id.vars ="Log_Copies.uL")
View(df1.m)


theme_update(plot.title = element_text(hjust = 0.5), 
             text=element_text(size=11, family = "Calibri"))

label_names <- c(
  'A_Ct_Mean' = "RSVA RT-qPCR Mean Ct",
  'B..CV..Coverage' = "A: CV Coverage (genome)",
  'D_percCV_HQ_Filtered_Mapped_reads' = "C: CV HQ Filtered Mapped Reads (genome)",
  'C_percCV_Average_Depth' = "B: CV Mean Depth(x) (genome)",
  'E_percCV_Ambiguous_bases' = "D: CV Ambiguous Bases (genome)"
)


ggplot(df1.m, 
       aes(Log_Copies.uL, 
           value, 
           colour = variable)) + 
  geom_point(show.legend = FALSE) +
  geom_smooth(method=lm, 
              formula = y ~ poly(x, 3), 
              se=F, 
              alpha=0.2, 
              linewidth=0.5, 
              show.legend = FALSE) +
  facet_wrap(~ variable, 
             ncol = 1, 
             scales = "free_y", 
             labeller = as_labeller(label_names)) +
  scale_colour_discrete() +
  geom_vline(data = data.frame(xint=3.8), 
             aes(xintercept = xint), 
             size = 0.6, 
             linetype="dotted", 
             col="red") +
  geom_segment(data = data.frame(yint=31,variable="A_Ct_Mean"), 
               linetype="dotted",
               aes(x = 3, 
                   y = yint, 
                   xend = 4, 
                   yend = yint),
               size = 0.6, col="red") +
  scale_y_continuous(breaks = breaks_fun) +
  labs(y="CV     CV        CV      CV      Ct") +
  xlab(bquote("Log"[10]~"(Copies/mL)")) +
  labs(title = "RSVA analytical WGS recovery") +
  xlim(3,6.7)


################################################################################
##RSVB analytical sensitivity QC CV graphs
################################################################################

dfB = read.table('RSVB_analytical.txt', 
                 sep='\t', 
                 header = T, 
                 colClasses = "numeric")
View(dfB)
dfB = data.frame(dfB)
View(dfB)
dfB.m = melt(dfB, 
             id.vars ="Log_Copies.uL")
View(dfB.m)

theme_update(plot.title = element_text(hjust = 0.5), 
             text=element_text(size=11, family = "Calibri"))

label_names <- c(
  'A_Ct_Mean' = "RSVB RT-qPCR Mean Ct",
  'B..CV..Coverage' = "E: CV Coverage (genome)",
  'D_percCV_HQ_Filtered_Mapped_reads' = "G: CV HQ Filtered Mapped Reads (genome)",
  'C_percCV_Average_Depth' = "F: CV Mean Depth(x) (genome)",
  'E_percCV_Ambiguous_bases' = "H: CV Ambiguous Bases (genome)"
)

breaks_fun <- function(x) {
  if (max(x) > 20) {
    seq(19, 32, 4)
  } else if (max(x) > 1.2) {
    seq(0, 2, 0.5)
  } else if (max(x) > 0.6) {
    seq(0, 2, 0.1)
  } else {
    seq(0, 1, 0.05)
  }
}




ggplot(dfB.m, aes(Log_Copies.uL, value, colour = variable)) + 
  geom_point(show.legend = FALSE) +
  geom_smooth(method=lm, formula = y ~ poly(x, 3), se=F, alpha=0.2, 
              linewidth=0.5, show.legend = FALSE) +
  facet_wrap(~ variable, ncol = 1, scales = "free_y", 
             labeller = as_labeller(label_names)) +
  scale_colour_discrete() +
  geom_vline(data = data.frame(xint=3.625), aes(xintercept = xint), 
             size = 0.6, linetype="dotted", col="red") +
  geom_segment(data = data.frame(yint=31,variable="A_Ct_Mean"), 
               linetype="dotted",
               aes(x = 3, y = yint, xend = 4, yend = yint),
               size = 0.6, col="red") +
  geom_blank(data = data.frame(xint=0.46,variable="A_Ct_Mean"),
             aes(x = xint, y = 19, xend = xint, yend = 32)) +
  scale_y_continuous(breaks = breaks_fun) +
  labs(y="CV      CV       CV      CV     Ct") +
  xlab(bquote("Log"[10]~"(Copies/mL)")) +
  labs(title = "RSVB analytical WGS recovery") +
  xlim(3,6.7)

ggsave("test.pdf",dpi=600)

ggplot_build(p)$data






################################################################################
#Boxplots for QC: av depth and HQ filtered mapped reads
################################################################################





df = read.table('RSV_clinical_QC_meta.txt',sep='\t', header = T, 
                colClasses = "numeric")

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}

ggplot(data=as.data.frame(df$avDepth), aes(x="", y=df$avDepth))+
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_boxplot(width=0.6, color="navy", fill="navy", alpha=0.2)+
  scale_y_continuous(labels=fancy_scientific) +
  geom_jitter(position = position_jitter(0.15), alpha=0.2) +
  labs(y = "RSV WGS average depth (x/genome)") +
  labs(x = "Median depth coverage: 4029x (n=1037 genomes)") +
  labs(title = "RSV WGS depth coverage, n=1037")


ggplot(data=as.data.frame(df$filtered_mapped_reads), 
       aes(x="", 
           y=df$avDepth))+
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_boxplot(width=0.6, color="darkgreen", fill="darkgreen", alpha=0.2)+
  scale_y_continuous(labels=fancy_scientific) +
  geom_jitter(position = position_jitter(0.15), alpha=0.2) +
  labs(y = "HQ mapped filtered reads") +
  labs(x = "Median HQ mapped filtered reads: 946066 reads (n=1037 genomes)") +
  labs(title = "RSV WGS average HQ mapped filtered reads, n=1037")


################################################################################
#RSV WGS sequence quality: Amb bases and N calls
################################################################################





ggplot(data=as.data.frame(df$percN), aes(x=df$percN, y=1-..y.., 
                                         colour = "skyblue1"))+
  stat_ecdf(size = 2) +
  stat_ecdf(data=as.data.frame(df$percdegens), 
            aes(x=df$percdegens, 
                y=1-..y.., 
                colour = "dodgerblue3"), 
            size = 2) +
  geom_segment(data = data.frame(xint=0.5,variable="percN"), 
               linetype="dashed",
               aes(x = xint, y = 0, xend = xint, yend = 0.05),
               size = 1, col="red") +
  geom_segment(data = data.frame(xint=0.023,variable="percdegens"), 
               linetype="dashed",
               aes(x = xint, y = 0, xend = xint, yend = 0.05),
               size = 1, col="red") +
  geom_segment(data = data.frame(yint=0.05,variable="percN"), 
               linetype="dashed",
               aes(x = 0, y = yint, xend = 0.5, yend = yint),
               size = 1, col="red") +
  coord_cartesian(ylim=c(0, 0.2)) +
  labs(y = "Empirical probability of occurrence, (%)") +
  labs(x = "% of ambiguous bases in consensus sequence (genome)") +
  labs(title = "RSV WG sequencing quality, n=1037") +
  scale_y_continuous(labels = scales::percent) +
  scale_colour_manual(name = 'type_of_ambiguity', 
                      values =c("skyblue1"="skyblue1",
                                          "dodgerblue3"="dodgerblue3"), 
                                          labels = c('N%','mixed-bases%')) +
  theme(legend.position=c(.85,.9))



################################################################################
#WGS recovery comparison for 'under five over sixty':
################################################################################


under_five <- read.table("under_five_Ct.txt")
View(under_five)
median(under_five$V1)

theme_update(text=element_text(size=14, family = "Calibri"))

ggplot(data=as.data.frame(under_five), aes(x="", y=under_five$V1))+
  stat_boxplot(geom = "errorbar", width = 0.1, size = 1, color = "red") +
  geom_boxplot(width=0.6, color="red", size = 1, fill="yellow", alpha=0.2)+
  geom_jitter(position = position_jitter(0.15), alpha=0.4) +
  labs(y="RSV RT-qPCR Ct") +
  labs(x="Ct under 5 y.o. (n=546)") +
  ylim(15,31)

over_sixty <- read.table("over_sixty_Ct.txt")
View(over_sixty)
median(over_sixty$V1)

ggplot(data=as.data.frame(over_sixty), aes(x="", y=over_sixty$V1))+
  stat_boxplot(geom = "errorbar", width = 0.1, size = 1, color = "red") +
  geom_boxplot(width=0.6, color="red", size = 1, fill="orange", alpha=0.2)+
  geom_jitter(position = position_jitter(0.15), alpha=0.4) +
  labs(y="RSV RT-qPCR Ct") +
  labs(x="Ct over 60 y.o. (n=172)") +
  ylim(15,31)


ks.test(under_five$V1, over_sixty$V1)

median_differences <- NULL
for (i in 1:10000) {
  median_differences <- c(median_differences, 
                          median(sample(over_sixty$V1, 
                                        replace=TRUE))-median(sample(under_five$V1, 
                                                                                   replace=TRUE)))
}

sprintf("Bootstrap 95 percent confidence interval for differece in median change lies between %06.3f and %06.3f", 
        quantile(median_differences,0.025), 
        quantile(median_differences,0.975))

hist(median_differences, 
     col="lightcyan2", 
     xlab="median differences")

wilcox.test(over_sixty$V1, under_five$V1)



####################################################################
##RSVA co-infection with RSVB
#######################################################################


df1 = read.table('RSVA_analytical_spec1.txt', 
                 sep='\t', 
                 header = T, 
                 colClasses = "numeric")
View(df1)
df1 = data.frame(df1)
View(df1)
df1.m = melt(df1, id.vars ="Log_Copies.uL")
View(df1.m)

df2 = read.table('RSVA_analytical_spec2.txt', 
                 sep='\t', 
                 header = T, 
                 colClasses = "numeric")
View(df2)
df2 = data.frame(df2)
View(df2)
df2.m = melt(df2, id.vars ="Log_Copies.uL", 
             value.name = "value2")
View(df2.m)

theme_update(plot.title = element_text(hjust = 0.5), 
             text=element_text(size=11, 
                               family = "Calibri"))

label_names <- c(
  'A_Ct_Mean' = "RSVA RT-qPCR Mean Ct",
  'B..CV..Coverage' = "CV Coverage (genome)",
  'D_percCV_HQ_Filtered_Mapped_reads' = "CV HQ Filtered Mapped Reads (genome)",
  'C_percCV_Average_Depth' = "CV Mean Depth(x) (genome)",
  'E_percCV_Ambiguous_bases' = "CV Ambiguous Bases (genome)"
)

breaks_fun <- function(x) {
  if (max(x) > 20) {
    seq(19, 32, 4)
  } else if (max(x) > 0.9) {
    seq(0, 2, 0.5)
  } else if (max(x) > 0.5) {
    seq(0, 2, 0.25)
  } else {
    seq(0, 1, 0.1)
  }
}


vars=c("RSVA A/-"="#225ea8", "RSV A/B co-infection"="#d7191c")

ggplot() + 
  scale_color_manual("",values = vars) +
  geom_point(data=df1.m, 
             aes(Log_Copies.uL, 
                             value, 
                             color="RSVA A/-"), 
             shape=18, 
             linewidth=0.5) +
  geom_smooth(data=df1.m, 
              aes(Log_Copies.uL, 
                  value, 
                  color="RSVA A/-"), 
              method=loess, 
              alpha=0.2, 
              linewidth=0.6) +
  geom_point(data = df2.m, 
             aes(Log_Copies.uL, 
                 value2, 
                 color="RSV A/B co-infection"), 
             shape=18, 
             linewidth=0.5) +
  geom_smooth(data = df2.m, 
              aes(Log_Copies.uL, 
                  value2, 
                  color="RSV A/B co-infection"), 
              method=loess, 
              linewidth=0.5) +
  scale_y_continuous(breaks = breaks_fun) +
  theme(legend.position = "bottom", 
        legend.spacing = unit(0, "pt"), 
        legend.margin = margin(r = 0, l = 0)) +
  labs(y="CV") +
  xlab(bquote("RSV A/B respective concentration (Log"[10]~"Copies/mL)")) +
  labs(title = "RSVA analytical specificity: RSVA-RSVB co-infection") +
  xlim(3.25,6.7) +
  facet_wrap(~ variable, 
             ncol = 1, 
             scales = "free_y", 
             labeller = as_labeller(label_names)) +
  scale_x_continuous(breaks = c(3.52, 4.52, 5.52, 6.52), 
                     labels = c("3.52/6.62", "4.52/5.62", "5.52/4.62", "6.52/3.62"))


################################################################################
##RSVB co-infection with RSVA
################################################################################

dfB1 = read.table('RSVB_analytical_spec1.txt', 
                  sep='\t', 
                  header = T, 
                  colClasses = "numeric")
View(dfB1)
dfB1 = data.frame(dfB1)
View(dfB1)
dfB1.m = melt(dfB1, 
              id.vars ="Log_Copies.uL", 
              value.name = "value3")
View(dfB1.m)

dfB2 = read.table('RSVB_analytical_spec2.txt', 
                  sep='\t', 
                  header = T, 
                  colClasses = "numeric")
View(dfB2)
dfB2 = data.frame(dfB2)
View(dfB2)
dfB2.m = melt(dfB2, 
              id.vars ="Log_Copies.uL", 
              value.name = "value4")
View(dfB2.m)

theme_update(plot.title = element_text(hjust = 0.5), 
             text=element_text(size=11, family = "Calibri"))

label_names <- c(
  'A_Ct_Mean' = "RSVB RT-qPCR Mean Ct",
  'B..CV..Coverage' = "CV Coverage (genome)",
  'D_percCV_HQ_Filtered_Mapped_reads' = "CV HQ Filtered Mapped Reads (genome)",
  'C_percCV_Average_Depth' = "CV Mean Depth(x) (genome)",
  'E_percCV_Ambiguous_bases' = "CV Ambiguous Bases (genome)"
)

breaks_fun <- function(x) {
  if (max(x) > 20) {
    seq(19, 32, 4)
  } else if (max(x) > 0.9) {
    seq(0, 2, 0.5)
  } else if (max(x) > 0.6) {
    seq(0, 2, 0.25)
  } else {
    seq(0, 1, 0.05)
  }
}


vars2 = c("RSVB B/-"="springgreen4", 
          "RSV B/A co-infection"="mediumorchid")

ggplot() + 
  scale_color_manual("",values = vars2) +
  geom_point(data=dfB1.m, aes(Log_Copies.uL, 
                              value3, 
                              color="RSVB B/-"), 
             shape=18) +
  stat_smooth(data=dfB1.m, 
              aes(Log_Copies.uL, 
                  value3, 
                  color="RSVB B/-"), 
              method=loess, 
              alpha=0.2, 
              linewidth=0.6, 
              level=0.95) +
  geom_point(data = dfB2.m, 
             aes(Log_Copies.uL, 
                 value4, 
                 color="RSV B/A co-infection"), 
             shape=18) +
  geom_smooth(data = dfB2.m, 
              aes(Log_Copies.uL, 
                  value4, 
                  color="RSV B/A co-infection"), 
              method=loess, 
              span=0.9,
              linewidth=0.5) +
  scale_y_continuous(breaks = breaks_fun) +
  theme(legend.position = "bottom", 
        legend.spacing = unit(0, "pt"), 
        legend.margin = margin(r = 0, l = 0)) +
  labs(y="CV") +
  xlab(bquote("RSV B/A respective concentration (Log"[10]~"Copies/mL)")) +
  labs(title = "RSVB analytical specificity: RSVB-RSVA co-infection") +
  xlim(3.25,6.7) +
  facet_wrap(~ variable, 
             ncol = 1, 
             scales = "free_y", 
             labeller = as_labeller(label_names)) +
  scale_x_continuous(breaks = c(3.62, 4.62, 5.62, 6.62), 
                     labels = c("3.62/6.52", "4.62/5.52", "5.62/4.52", "6.62/3.52"))



################################################################################
## RSVA spiked with other pathogens
################################################################################

dfAS = read.table('RSVA_analytical_spiked1.txt', 
                  sep='\t', 
                  header = T, 
                  colClasses = "numeric")
View(dfAS)
dfAS = data.frame(dfAS)
View(dfAS)
dfAS.m = melt(dfAS, 
              id.vars ="Log_Copies.uL", 
              value.name = "valueAS")
View(dfAS.m)


theme_update(plot.title = element_text(hjust = 0.5), 
             text=element_text(size=11, 
                               family = "Calibri"))


breaks_fun <- function(x) {
  if (max(x) > 20) {
    seq(19, 32, 4)
  } else if (max(x) > 0.9) {
    seq(0, 2, 0.5)
  } else if (max(x) > 0.5) {
    seq(0, 2, 0.25)
  } else {
    seq(0, 1, 0.05)
  }
}


varsAS=c("RSVA"="green4", "RSVA spiked"="maroon3")

ggplot() + 
  scale_color_manual("",values = varsAS) +
  geom_point(data=df1.m, 
             aes(Log_Copies.uL, 
                 value, 
                 color="RSVA"), 
             shape=18) +
  stat_smooth(data=df1.m, 
              aes(Log_Copies.uL, 
                  value, 
                  color="RSVA"), 
              method=loess, 
              alpha=0.2, 
              linewidth=0.6, 
              level=0.95) +
  geom_point(data = dfAS.m, 
             aes(Log_Copies.uL, 
                 valueAS, 
                 color="RSVA spiked"), 
             shape=18) +
  geom_smooth(data = dfAS.m, 
              aes(Log_Copies.uL, 
                  valueAS, 
                  color="RSVA spiked"), 
              method=loess, 
              span=0.9,
              linewidth=0.5) +
  scale_y_continuous(breaks = breaks_fun) +
  theme(legend.position = "bottom", 
        legend.spacing = unit(0, "pt"), 
        legend.margin = margin(r = 0, l = 0)) +
  labs(y="CV") +
  xlab(bquote("RSVA (Log"[10]~"Copies/mL)")) +
  labs(title = "RSVA analytical specificity: RSVA spiked with other pathogens") +
  xlim(3.25,6.7) +
  facet_wrap(~ variable, 
             ncol = 1, 
             scales = "free_y", 
             labeller = as_labeller(label_names)) +
  scale_x_continuous(breaks = c(3.52, 4.52, 5.52, 6.52))


################################################################################
## RSVB spiked with other pathogens
################################################################################

dfBS = read.table('RSVB_analytical_spiked1.txt', 
                  sep='\t', 
                  header = T, 
                  colClasses = "numeric")
View(dfBS)
dfBS = data.frame(dfBS)
View(dfBS)
dfBS.m = melt(dfBS, 
              id.vars ="Log_Copies.uL", 
              value.name = "valueBS")
View(dfBS.m)


theme_update(plot.title = element_text(hjust = 0.5), 
             text=element_text(size=11, family = "Calibri"))


breaks_fun <- function(x) {
  if (max(x) > 20) {
    seq(19, 32, 4)
  } else if (max(x) > 0.9) {
    seq(0, 2, 0.5)
  } else if (max(x) > 0.5) {
    seq(0, 2, 0.25)
  } else {
    seq(0, 1, 0.05)
  }
}


varsBS=c("RSVB"="chartreuse4", "RSVB spiked"="hotpink2")

ggplot() + 
  scale_color_manual("",values = varsBS) +
  geom_point(data=dfB1.m, 
             aes(Log_Copies.uL, 
                 value3, 
                 color="RSVB"), 
             shape=18) +
  stat_smooth(data=dfB1.m, 
              aes(Log_Copies.uL, 
                  value3, 
                  color="RSVB"), 
              method=loess, 
              alpha=0.2, 
              linewidth=0.6, 
              level=0.95) +
  geom_point(data = dfBS.m, 
             aes(Log_Copies.uL, 
                 valueBS, 
                 color="RSVB spiked"), 
             shape=18) +
  geom_smooth(data = dfBS.m, 
              aes(Log_Copies.uL, 
                  valueBS,
                  color="RSVB spiked"), 
              method=loess, 
              span=0.9,
              linewidth=0.5) +
  scale_y_continuous(breaks = breaks_fun) +
  theme(legend.position = "bottom", 
        legend.spacing = unit(0, "pt"), 
        legend.margin = margin(r = 0, l = 0)) +
  labs(y="CV") +
  xlab(bquote("RSVB (Log"[10]~"Copies/mL)")) +
  labs(title = "RSVB analytical specificity: RSVB spiked with other pathogens") +
  xlim(3.25,6.7) +
  facet_wrap(~ variable, 
             ncol = 1, 
             scales = "free_y", 
             labeller = as_labeller(label_names)) +
  scale_x_continuous(breaks = c(3.62, 4.62, 5.62, 6.62))
