##### Cape Cod National Seashore coyote 
# Genetics, abundance/density, and occupancy R Script

#### Genetic summary stats

#packages 
library(adegenet)
library(tidyverse)

# 12 microsatellite loci (SRY is a sex marker, excluded)
loci <- c("AHTh171coy", "CPH18", "CXX-468", "CXX-602coy", "REN54P11",
          "AHT137", "AHT142", "CXX-279", "CXX-374", "INU055",
          "REN162C04", "REN169O18coyote")

#data
genotypes_22_23 <- read.csv("genotypes_22_23.csv")
genotypes_22    <- read.csv("genotypes_22.csv")
genotypes_23    <- read.csv("genotypes_23.csv")

#keep the 12 microsatellite loci (locusA-locusL)
old_names <- paste0("locus", LETTERS[1:12])
genotypes_22_23 <- genotypes_22_23[, old_names]
genotypes_22    <- genotypes_22[, old_names]
genotypes_23    <- genotypes_23[, old_names]

#rename to the real locus names (same order as columns A-L)
new_names <- c("AHTh171coy", "CPH18", "CXX-468", "CXX-602coy", "REN54P11",
               "AHT137", "AHT142", "CXX-279", "CXX-374", "INU055",
               "REN162C04", "REN169O18coyote")
names(genotypes_22_23) <- new_names
names(genotypes_22)    <- new_names
names(genotypes_23)    <- new_names

#genind
CCNS22_23_genind <- df2genind(genotypes_22_23, sep = ":", ploidy = 2)
CCNS22_genind    <- df2genind(genotypes_22,    sep = ":", ploidy = 2)
CCNS23_genind    <- df2genind(genotypes_23,    sep = ":", ploidy = 2)

#check
CCNS22_23_genind$loc.n.all

# Genetics summary statistics-----
#summary
summary(CCNS22_23_genind)
summary(CCNS22_genind)
summary(CCNS23_genind)

#average number of alleles 
mean(CCNS22_23_genind$loc.n.all) # 8.666667
min(CCNS22_23_genind$loc.n.all)
max(CCNS22_23_genind$loc.n.all)
mean(CCNS22_genind$loc.n.all) #7.666667
min(CCNS22_genind$loc.n.all)
max(CCNS22_genind$loc.n.all)
mean(CCNS23_genind$loc.n.all) #7.75
min(CCNS23_genind$loc.n.all)
max(CCNS23_genind$loc.n.all)

#expected heterozygosity for all individuals sampled in 2023 and 2023
mean(CCNS22_23_genind_summary$Hexp) # 0.7299793
min(CCNS22_23_genind_summary$Hexp)
max(CCNS22_23_genind_summary$Hexp)
mean(CCNS22_genind_summary$Hexp)  # 0.721876
min(CCNS22_genind_summary$Hexp)
max(CCNS22_genind_summary$Hexp)
max(CCNS22_23_genind_summary$Hexp)
mean(CCNS23_genind_summary$Hexp) # 0.7331315
min(CCNS23_genind_summary$Hexp)
max(CCNS23_genind_summary$Hexp)



#observed heterozygosity for all individuals sampled in 2023 and 2023
mean(CCNS22_23_genind_summary$Hobs) #0.6694139
min(CCNS22_23_genind_summary$Hobs)
max(CCNS22_23_genind_summary$Hobs)
mean(CCNS22_genind_summary$Hobs) #0.7090643
min(CCNS22_genind_summary$Hobs)
max(CCNS22_genind_summary$Hobs)
mean(CCNS23_genind_summary$Hobs) #0.6029412
min(CCNS23_genind_summary$Hobs)
max(CCNS23_genind_summary$Hobs)


# per locus 
# number of alleles per locus for each dataset

alleles_tab <- data.frame(
  Locus       = names(CCNS22_23_genind$loc.n.all),
  "2022"      = as.numeric(CCNS22_genind$loc.n.all),
  "2023"      = as.numeric(CCNS23_genind$loc.n.all),
  check.names = FALSE
)

alleles_tab
write.csv(alleles_tab, "Table_alleles_per_locus.csv", row.names = FALSE)

# PID 
datasets <- list("2022" = CCNS22_genind,
                 "2023" = CCNS23_genind)

N <- 1000   # number of simulations
set.seed(1)

PIDcombo <- data.frame()
report   <- data.frame()

for (yr in names(datasets)) {
  
  g <- datasets[[yr]]
  
  #allele frequencies for each locus -> list of vectors
  freqs      <- colMeans(tab(g, freq = TRUE), na.rm = TRUE)
  al.freq.ls <- list()
  for (i in 1:nLoc(g)) {
    al.freq.ls[[i]] <- as.numeric(freqs[g$loc.fac == levels(g$loc.fac)[i]])
  }
  names(al.freq.ls) <- levels(g$loc.fac)
  
  x <- length(al.freq.ls)   # max number of loci (12)
  
  PID.pop.tbl     <- matrix(data = NA, nrow = x, ncol = 3)
  PID.pop.sib.tbl <- matrix(data = NA, nrow = x, ncol = 3)
  
  for (j in 1:x) {                       # for each count of loci
    
    PIpop.loci     <- NULL
    PIpop.sib.loci <- NULL
    
    for (h in 1:N) {                     # do N simulations
      samp      <- sample(al.freq.ls, j, replace = FALSE)
      PIpop     <- NULL
      PIpop.sib <- NULL
      
      for (i in 1:length(samp)) {        # Waits et al. 2001
        p  <- samp[[i]]
        r1 <- sum(p^4)
        r2 <- 4 * outer(p, p)^2          # (2*p_i*p_j)^2
        r  <- sum(r2[lower.tri(r2)]) + r1
        s  <- 0.25 + (0.5 * sum(p^2)) + (0.5 * sum(p^2)^2) - (0.25 * sum(p^4))
        
        PIpop     <- c(PIpop, r)
        PIpop.sib <- c(PIpop.sib, s)
      }
      PIpop.loci     <- c(PIpop.loci, prod(PIpop))
      PIpop.sib.loci <- c(PIpop.sib.loci, prod(PIpop.sib))
    }
    
    PID.pop.tbl[j, ]     <- c(j, mean(PIpop.loci),     sd(PIpop.loci))
    PID.pop.sib.tbl[j, ] <- c(j, mean(PIpop.sib.loci), sd(PIpop.sib.loci))
  }
  
  PID.pop.tbl     <- as.data.frame(PID.pop.tbl)
  PID.pop.sib.tbl <- as.data.frame(PID.pop.sib.tbl)
  colnames(PID.pop.tbl)     <- c("Loci", "Probability", "SD")
  colnames(PID.pop.sib.tbl) <- c("Loci", "Probability", "SD")
  
  PID.pop.tbl$Type     <- "PID"
  PID.pop.sib.tbl$Type <- "PIDsib"
  PID.pop.tbl$Year     <- yr
  PID.pop.sib.tbl$Year <- yr
  
  PIDcombo <- rbind(PIDcombo, PID.pop.tbl, PID.pop.sib.tbl)
  
  # values to report (all 12 loci)
  report <- rbind(report,
                  data.frame(Year   = yr,
                             PID    = PID.pop.tbl$Probability[10],
                             PIDsib = PID.pop.sib.tbl$Probability[10]))
}

report   # PID and PIDsib for 10 loci

#### plot
ggplot(PIDcombo, aes(x = Loci, y = Probability, color = Type, linetype = Type)) +
  geom_line(linewidth = 1.5) +
  facet_wrap(~ Year) +
  scale_x_continuous(breaks = 1:12) +
  scale_y_log10() +
  scale_color_manual(values = c("gray68", "black")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  geom_hline(yintercept = 0.001, color = "red", linewidth = 1) +
  ylab("Probability (log10 scale)") +
  theme(text = element_text(size = 16, family = "sans"),
        panel.background = element_blank(),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(linewidth = 1),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.key.width = unit(2, "cm"),
        axis.text = element_text(size = 14),
        strip.background = element_blank(),
        strip.text = element_text(size = 16),
        legend.position = "top")
