## Calculations of individual PCB fractions in well-plates
# adipose cells + FBS (albumin) + air + freely dissolved
# 3 main caculations: (1) Cabinet mixture, (2) Aroclor 1016 (top 12)
# and (3) Aroclor 1254 (top 12)

# Install packages
install.packages("readxl")
install.packages("reshape2")
install.packages("ggplot2")

# load libraries
library(readxl)
library(reshape2)
library(ggplot2)

# Partitioning from:
# https://www.ufz.de/index.php?en=31698&contentonly=1&m=0&lserd_data[mvc]=Public/start

#Fixed parameters
# 24-well plate https://shop.gbo.com/en/usa/products/bioscience/cell-culture-products/cellstar-cell-culture-multiwell-plates/662160.html
# Multiwell dimensions
Vt <- 3.3/1000 # Total volume L
Vm <- 1/1000 # Medium volume L
Va <- Vt-Vm # Air volume L

# Fraction without adipose cells --------------------------------
# Function to calculate fractions
fraction = function(logKa.w, dUaw, logKpro.w, logKalb.w, R,
                    tst, texp) {
  
  # Albumin concentration from FBS
  # Concentration of albumin in FBS
  V.FBS.h <- 0.15 # mL 15%
  V.FBS.l <- 0.005 # mL 0.5%
  C.alb.initial <- 23 # mg/mL
  C.prot.initial <- (38 - C.alb.initial) # mg/mL
  C.alb.h <- C.alb.initial*V.FBS.h/Vm/1000/1000 # kg/L (15%)
  C.alb.l <- C.alb.initial*V.FBS.l/Vm/1000/1000 # kg/L (0.5%)
  dalb <- 1 # kg/L albumin density
  C.alb.h <- C.alb.h/dalb # Lalb/Lwater
  C.alb.l <- C.alb.l/dalb # Lalb/Lwater
  # Protein concentration from FBS
  C.prot.med.h <- C.prot.initial*V.FBS.h/Vm/1000/1000 # kg/L
  C.prot.med.l <- C.prot.initial*V.FBS.l/Vm/1000/1000 # kg/L
  dprot <- 1.43 # kg/L protein density. ref: https://pubmed.ncbi.nlm.nih.gov/10930825/
  C.prot.med.h <- C.prot.med.h/dprot # Lprot/Lwater
  C.prot.med.l <- C.prot.med.l/dprot # Lprot/Lwater
  
  # Temperature correction for Kaw
  Ka.w.t <- 10^(logKa.w)*exp(-dUaw/R*(1/texp-1/tst)) # Ka.w corrected by water and air temps
  logKa.w.t <- log10(Ka.w.t)
  
  # Fraction calculation
  # High FBS (15%)
  den.alb.h <- 1 + 10^(logKalb.w)*C.alb.h + 10^(logKpro.w)*C.prot.med.h +
    10^(logKa.w.t)*Va/Vm
  f.dis.alb.h <- 1/den.alb.h # freely dissolved fraction
  f.alb.alb.h <- 10^(logKalb.w)*C.alb.h/den.alb.h # albumin fraction from FBS
  f.prot.alb.h <- 10^(logKpro.w)*C.prot.med.h/den.alb.h # protein fraction from FBS
  f.air.alb.h <- 10^(logKa.w.t)*Va/den.alb.h/Vm # air fraction
  
  # Low FBS (0.5%)
  den.alb.l <- 1 + 10^(logKalb.w)*C.alb.l + 10^(logKpro.w)*C.prot.med.l +
    10^(logKa.w.t)*Va/Vm
  f.dis.alb.l <- 1/den.alb.l # freely dissolved fraction
  f.alb.alb.l <- 10^(logKalb.w)*C.alb.l/den.alb.l # albumin fraction from FBS
  f.prot.alb.l <- 10^(logKpro.w)*C.prot.med.l/den.alb.l # protein fraction from FBS
  f.air.alb.l <- 10^(logKa.w.t)*Va/den.alb.l/Vm # air fraction
  
  # No FBS
  den.alb.0 <- 1 + 10^(logKa.w.t)*Va/Vm
  f.dis.alb.0 <- 1/den.alb.0 # freely dissolved fraction
  f.air.alb.0 <- 10^(logKa.w.t)*Va/den.alb.0/Vm # air fraction
  
  frac <- c(f.dis.alb.h, f.alb.alb.h, f.prot.alb.h, f.air.alb.h,
            f.dis.alb.l, f.alb.alb.l, f.prot.alb.l, f.air.alb.l,
            f.dis.alb.0, f.air.alb.0)
}

# Cabinet mixture without adipose cells -----------------------------------
# Read data.xlsx
d <- data.frame(read_excel("DataAdi.xlsx", sheet = "cabinet",
                           col_names = TRUE, col_types = NULL))

# Name parameters
congener <- d$congener
logKa.w <- d$logKair.water
dUaw <- d$dUaw
logKlip.w <- d$logKlipid.water
logKpro.w <- d$logKprotein.water
logKalb.w <- d$logKalbumin.water
R <- d$R
tst <- d$tst
texp <- d$texp

num.congener <- length(congener)
result <- NULL
for (i in 1:num.congener) {
  result <- rbind(result, fraction(logKa.w[i], dUaw[i],
                                   logKpro.w[i], logKalb.w[i], R[i],
                                   tst[i], texp[i]))
}

final.result <- data.frame(congener, result)
names(final.result) <- c("congener", "frac.dis.FBS(15%)", "frac.alb.FBS(15%)",
                         "frac.prot.FBS(15%)", "frac.air.FBS(15%)",
                         "frac.dis.FBS(0.5%)", "frac.alb.FBS(0.5%)",
                         "frac.prot.FBS(0.5%)", "frac.air.FBS(0.5%)",
                         "frac.dis.FBS(0%)", "frac.air.FBS(0%)")

# Export results
write.csv(final.result, file = "CabinetNoAdipose.csv")
# To make plots, jump to section "Plots without adipose cells"

# Aroclor 1016 (top 12) without adipose cells -----------------------------
# Read data.xlsx
d.A1016 <- data.frame(read_excel("DataAdi.xlsx", sheet = "Aroclor1016",
                                 col_names = TRUE, col_types = NULL))

# Name parameters
congener <- d.A1016$congener
logKa.w <- d.A1016$logKair.water
dUaw <- d.A1016$dUaw
logKlip.w <- d.A1016$logKlipid.water
logKpro.w <- d.A1016$logKprotein.water
logKalb.w <- d.A1016$logKalbumin.water
R <- d.A1016$R
tst <- d.A1016$tst
texp <- d.A1016$texp

num.congener <- length(congener)
result <- NULL
for (i in 1:num.congener) {
  result <- rbind(result, fraction(logKa.w[i], dUaw[i],
                                   logKpro.w[i], logKalb.w[i], R[i],
                                   tst[i], texp[i]))
}

final.result <- data.frame(congener, result)
names(final.result) <- c("congener", "frac.dis.FBS(15%)", "frac.alb.FBS(15%)",
                         "frac.prot.FBS(15%)", "frac.air.FBS(15%)",
                         "frac.dis.FBS(0.5%)", "frac.alb.FBS(0.5%)",
                         "frac.prot.FBS(0.5%)", "frac.air.FBS(0.5%)",
                         "frac.dis.FBS(0%)", "frac.air.FBS(0%)")

# Export results
write.csv(final.result, file = "A1016NoAdiposeV3.csv")
# To make plots, jump to section "Plots without adipose cells"

# Aroclor 1254 (top 12) without adipose cells -----------------------------
# Read data.xlsx
d.A1254 <- data.frame(read_excel("DataAdi.xlsx", sheet = "Aroclor1254",
                                 col_names = TRUE, col_types = NULL))

# Name parameters
congener <- d.A1254$congener
logKa.w <- d.A1254$logKair.water
dUaw <- d.A1254$dUaw
logKlip.w <- d.A1254$logKlipid.water
logKpro.w <- d.A1254$logKprotein.water
logKalb.w <- d.A1254$logKalbumin.water
R <- d.A1254$R
tst <- d.A1254$tst
texp <- d.A1254$texp

num.congener <- length(congener)
result <- NULL
for (i in 1:num.congener) {
  result <- rbind(result, fraction(logKa.w[i], dUaw[i],
                                   logKpro.w[i], logKalb.w[i], R[i],
                                   tst[i], texp[i]))
}

final.result <- data.frame(congener, result)
names(final.result) <- c("congener", "frac.dis.FBS(15%)", "frac.alb.FBS(15%)",
                         "frac.prot.FBS(15%)", "frac.air.FBS(15%)",
                         "frac.dis.FBS(0.5%)", "frac.alb.FBS(0.5%)",
                         "frac.prot.FBS(0.5%)", "frac.air.FBS(0.5%)",
                         "frac.dis.FBS(0%)", "frac.air.FBS(0%)")

# Export results
write.csv(final.result, file = "A1254NoAdipose.csv")
# To make plots, jump to section "Plots without adipose cells"

# Plots without adipose cells ---------------------------------------------
# final.result needs to be changed, depending on the plot
# Create data.frame with needed fractions
# (1) FSB = 15%
p.FBS.h <- final.result[,!names(final.result) %in% c("frac.dis.FBS(0.5%)", "frac.alb.FBS(0.5%)",
                                                     "frac.prot.FBS(0.5%)", "frac.air.FBS(0.5%)",
                                                     "frac.dis.FBS(0%)", "frac.air.FBS(0%)")]
# Transform data.frame p.1 to 3 column data.frame
p.FBS.h <- melt(p.FBS.h, id.var = c("congener"),
                variable.name = "phase", value.name = "fraction")

# For A1016 and A1254, congeners need to be changed
# Cabinet mixture
p.FBS.h$congener <- factor(p.FBS.h$congener,
                           levels = c('PCB11', 'PCB47', 'PCB51',
                                      'PCB68'))
# A1016
p.FBS.h$congener <- factor(p.FBS.h$congener,
                           levels = c('PCB4', 'PCB8', 'PCB16',
                                      'PCB17', 'PCB18', 'PCB22',
                                      'PCB28', 'PCB31', 'PCB33',
                                      'PCB44', 'PCB49', 'PCB52'))

# A1254
p.FBS.h$congener <- factor(p.FBS.h$congener,
                           levels = c('PCB66', 'PCB70', 'PCB85',
                                      'PCB87', 'PCB97', 'PCB99',
                                      'PCB101', 'PCB105', 'PCB110',
                                      'PCB118', 'PCB138', 'PCB153'))

# Organize fraction to be displayed in plot
p.FBS.h$phase <- factor(p.FBS.h$phase,
                        levels = c('frac.air.FBS(15%)', 'frac.alb.FBS(15%)',
                                   'frac.prot.FBS(15%)', 'frac.dis.FBS(15%)'))
# Plot
ggplot(p.FBS.h, aes(x = congener, y = fraction, fill = phase)) + 
  geom_bar(stat = "identity", col = "white") +
  scale_fill_manual(labels = c("Air" , "Albumin-FBS (15%)", "Protein-FBS (15%)",
                               "Medium"),
                    values = c("deepskyblue", "lightgrey", "coral4", "red")) +
  theme_classic() +
  xlab(expression(bold(""))) +
  ylab(expression("Fraction in well")) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.title.x = element_text(face = "bold", size = 8))

# (2) FSB = 0.5%
p.FBS.l <- final.result[,!names(final.result) %in% c("frac.dis.FBS(15%)", "frac.alb.FBS(15%)",
                                                     "frac.prot.FBS(15%)", "frac.air.FBS(15%)",
                                                     "frac.dis.FBS(0%)", "frac.air.FBS(0%)")]
# Transform data.frame p.1 to 3 column data.frame
p.FBS.l <- melt(p.FBS.l, id.var = c("congener"),
                variable.name = "phase", value.name = "fraction")

# For A1016 and A1254, congeners need to be changed
# Cabinet mixture
p.FBS.l$congener <- factor(p.FBS.l$congener,
                           levels = c('PCB11', 'PCB47', 'PCB51',
                                      'PCB68'))

# A1016
p.FBS.l$congener <- factor(p.FBS.l$congener,
                           levels = c('PCB4', 'PCB8', 'PCB16',
                                      'PCB17', 'PCB18', 'PCB22',
                                      'PCB28', 'PCB31', 'PCB33',
                                      'PCB44', 'PCB49', 'PCB52'))

# A1254
p.FBS.l$congener <- factor(p.FBS.l$congener,
                           levels = c('PCB66', 'PCB70', 'PCB85',
                                      'PCB87', 'PCB97', 'PCB99',
                                      'PCB101', 'PCB105', 'PCB110',
                                      'PCB118', 'PCB138', 'PCB153'))

# Organize fraction to be displayed in plot
p.FBS.l$phase <- factor(p.FBS.l$phase,
                        levels = c('frac.air.FBS(0.5%)', 'frac.alb.FBS(0.5%)',
                                   'frac.prot.FBS(0.5%)', 'frac.dis.FBS(0.5%)'))

# Plot
ggplot(p.FBS.l, aes(x = congener, y = fraction, fill = phase)) + 
  geom_bar(stat = "identity", col = "white") +
  scale_fill_manual(labels = c("Air" , "Albumin-FBS (0.5%)", "Protein-FBS (0.5%)",
                               "Medium"),
                    values = c("deepskyblue", "lightgrey", "coral4", "red")) +
  theme_classic() +
  xlab(expression(bold(""))) +
  ylab(expression("Fraction in well")) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.title.x = element_text(face = "bold", size = 8))

# (3) FSB = 0.0%
p.FBS.0 <- final.result[,!names(final.result) %in% c("frac.dis.FBS(15%)", "frac.alb.FBS(15%)",
                                                     "frac.prot.FBS(15%)", "frac.air.FBS(15%)",
                                                     "frac.dis.FBS(0.5%)", "frac.alb.FBS(0.5%)",
                                                     "frac.prot.FBS(0.5%)", "frac.air.FBS(0.5%)")]

# Transform data.frame p.1 to 3 column data.frame
p.FBS.0 <- melt(p.FBS.0, id.var = c("congener"),
                variable.name = "phase", value.name = "fraction")

# For A1016 and A1254, congeners need to be changed
# Cabinet mixture
p.FBS.0$congener <- factor(p.FBS.0$congener,
                           levels = c('PCB11', 'PCB47', 'PCB51',
                                      'PCB68'))

# A1016
p.FBS.0$congener <- factor(p.FBS.0$congener,
                           levels = c('PCB4', 'PCB8', 'PCB16',
                                      'PCB17', 'PCB18', 'PCB22',
                                      'PCB28', 'PCB31', 'PCB33',
                                      'PCB44', 'PCB49', 'PCB52'))

# A1254
p.FBS.0$congener <- factor(p.FBS.0$congener,
                           levels = c('PCB66', 'PCB70', 'PCB85',
                                      'PCB87', 'PCB97', 'PCB99',
                                      'PCB101', 'PCB105', 'PCB110',
                                      'PCB118', 'PCB138', 'PCB153'))

# Organize fraction to be displayed in plot
p.FBS.0$phase <- factor(p.FBS.0$phase,
                        levels = c('frac.air.FBS(0%)', 'frac.dis.FBS(0%)'))

# Plot
ggplot(p.FBS.0, aes(x = congener, y = fraction, fill = phase)) + 
  geom_bar(stat = "identity", col = "white") +
  scale_fill_manual(labels = c("Air" , "Medium"),
                    values = c("deepskyblue", "red")) +
  theme_classic() +
  xlab(expression(bold(""))) +
  ylab(expression("Fraction in well")) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.title.x = element_text(face = "bold", size = 8))

# Fraction with adipose cells -----------------------------------
# Function to calculate fractions
fractionAdi = function(logKa.w, dUaw, logKlip.w, logKpro.w, logKalb.w,
                    R, tst, texp) {
  
  # Concentrations need to be in volume (L), not mass
  Adi <- 10^-7 # kg/well
  C.adi <- Adi/Vm # Concentration of cell per well #kg/L
  C.lip.adi <- C.adi*0.6 # lipid content 60% kg/L
  dlip <- 0.905 # kg/L lipid density. ref: https://pubmed.ncbi.nlm.nih.gov/8148928/
  C.lip.adi <- C.lip.adi/dlip # L/L
  C.prot.adi <- C.adi*0.004 # protein content 0.4% kg/L
  dprot <- 1.43 # kg/L protein density. ref: https://pubmed.ncbi.nlm.nih.gov/10930825/
  C.prot.adi <- C.prot.adi/dprot # L/L
  C.water.adi <- C.adi*0.4 # water content 40% kg/L
  dwater <- 0.99127 # kg/L water density at 37 C
  C.water.adi <- C.water.adi/dwater # L/L
  V.water.adi <- C.water.adi*Adi/10^6 # L water inside cell
  
  # Concentration of albumin in FBS
  V.FBS.h <- 0.15 # mL 15%
  V.FBS.l <- 0.005 # mL 0.5%
  C.alb.initial <- 23 # mg/mL
  C.prot.initial <- (38 - C.alb.initial) # mg/mL
  C.alb.h <- C.alb.initial*V.FBS.h/Vm/1000/1000 # kg/L (15%)
  C.alb.l <- C.alb.initial*V.FBS.l/Vm/1000/1000 # kg/L (0.5%)
  dalb <- 1 # kg/L albumin density
  C.alb.h <- C.alb.h/dalb # Lalb/Lwater
  C.alb.l <- C.alb.l/dalb # Lalb/Lwater
  # Protein concentration from FBS
  C.prot.med.h <- C.prot.initial*V.FBS.h/Vm/1000/1000 # kg/L
  C.prot.med.l <- C.prot.initial*V.FBS.l/Vm/1000/1000 # kg/L
  dprot <- 1.43 # kg/L protein density. ref: https://pubmed.ncbi.nlm.nih.gov/10930825/
  C.prot.med.h <- C.prot.med.h/dprot # Lprot/Lwater
  C.prot.med.l <- C.prot.med.l/dprot # Lprot/Lwater

  # Temperature correction for Kaw
  Ka.w.t <- 10^(logKa.w)*exp(-dUaw/R*(1/texp-1/tst)) # Ka.w corrected by water and air temps
  logKa.w.t <- log10(Ka.w.t)
  
  # Fraction calculation
  # High FBS (15%)
  den.alb.h.adi <- 1 + 10^(logKalb.w)*C.alb.h + 10^(logKpro.w)*C.prot.med.h +
    10^(logKlip.w)*C.lip.adi + 10^(logKpro.w)*C.prot.adi + 
    10^(logKa.w.t)*Va/(Vm-V.water.adi)
  f.dis.alb.h.adi <- 1/den.alb.h.adi # freely dissolved fraction
  f.alb.alb.h.adi <- 10^(logKalb.w)*C.alb.h/den.alb.h.adi # albumin fraction from FBS
  f.prot.alb.h.adi <- 10^(logKpro.w)*C.prot.med.h/den.alb.h.adi # protein fraction from FBS
  f.prot.adi.alb.h.adi <- 10^(logKpro.w)*C.prot.adi/den.alb.h.adi # protein fraction from adipose cells
  f.lip.adi.alb.h.adi <- 10^(logKlip.w)*C.lip.adi/den.alb.h.adi # lipid fraction from adipose cells
  f.air.alb.h.adi <- 10^(logKa.w.t)*Va/Vm/den.alb.h.adi # air fraction
  
  # Low FBS (0.5%)
  den.alb.l.adi <- 1 + 10^(logKalb.w)*C.alb.l + 10^(logKpro.w)*C.prot.med.l +
    10^(logKlip.w)*C.lip.adi + 10^(logKpro.w)*C.prot.adi + 
    10^(logKa.w.t)*Va/(Vm-V.water.adi)
  f.dis.alb.l.adi <- 1/den.alb.l.adi # freely dissolved fraction
  f.alb.alb.l.adi <- 10^(logKalb.w)*C.alb.l/den.alb.l.adi # albumin fraction from FBS
  f.prot.alb.l.adi <- 10^(logKpro.w)*C.prot.med.l/den.alb.l.adi # protein fraction from FBS
  f.prot.adi.alb.l.adi <- 10^(logKpro.w)*C.prot.adi/den.alb.l.adi # protein fraction from adipose cells
  f.lip.adi.alb.l.adi <- 10^(logKlip.w)*C.lip.adi/den.alb.l.adi # lipid fraction from adipose cells
  f.air.alb.l.adi <- 10^(logKa.w.t)*Va/Vm/den.alb.l.adi # air fraction
  
  # No FBS
  den.alb.0.adi <- 1 + 10^(logKlip.w)*C.lip.adi + 10^(logKpro.w)*C.prot.adi + 
    10^(logKa.w.t)*Va/(Vm-V.water.adi)
  f.dis.alb.0.adi <- 1/den.alb.0.adi # freely dissolved fraction
  f.prot.adi.alb.0.adi <- 10^(logKpro.w)*C.prot.adi/den.alb.0.adi # protein fraction from adipose cells
  f.lip.adi.alb.0.adi <- 10^(logKlip.w)*C.lip.adi/den.alb.0.adi # lipid fraction from adipose cells
  f.air.alb.0.adi <- 10^(logKa.w.t)*Va/Vm/den.alb.0.adi # air fraction
  
  frac.adi <- c(f.dis.alb.h.adi, f.alb.alb.h.adi, f.prot.alb.h.adi,
            f.prot.adi.alb.h.adi, f.lip.adi.alb.h.adi, f.air.alb.h.adi,
            f.dis.alb.l.adi, f.alb.alb.l.adi, f.prot.alb.l.adi,
            f.prot.adi.alb.l.adi, f.lip.adi.alb.l.adi, f.air.alb.l.adi,
            f.dis.alb.0.adi, f.prot.adi.alb.0.adi, f.lip.adi.alb.0.adi,
            f.air.alb.0.adi)
}

# Cabinet mixture with adipose cells --------------------------------------
# Read data.xlsx
d <- data.frame(read_excel("DataAdi.xlsx", sheet = "cabinet",
                           col_names = TRUE, col_types = NULL))

# Name parameters
congener <- d$congener
logKa.w <- d$logKair.water
dUaw <- d$dUaw
logKlip.w <- d$logKlipid.water
logKpro.w <- d$logKprotein.water
logKalb.w <- d$logKalbumin.water
R <- d$R
tst <- d$tst
texp <- d$texp

num.congener <- length(congener)
resultAdi <- NULL
for (i in 1:num.congener) {
  resultAdi <- rbind(resultAdi, fractionAdi(logKa.w[i], dUaw[i], logKlip.w[i],
                                   logKpro.w[i], logKalb.w[i], R[i],
                                   tst[i], texp[i]))
}

final.resultAdi <- data.frame(congener, resultAdi)
names(final.resultAdi) <- c("congener", "frac.dis.adi.FBS(15%)",
                            "frac.alb.adi.FBS(15%)", "frac.prot.adi.FBS(15%)",
                            "frac.prot-cell.adi.FBS(15%)", "frac.lip-cell.adi.FBS(15%)",
                            "frac.air.adi.FBS(15%)", "frac.dis.adi.FBS(0.5%)",
                            "frac.alb.adi.FBS(0.5%)", "frac.prot.adi.FBS(0.5%)",
                            "frac.prot-cell.adi.FBS(0.5%)", "frac.lip-cell.adi.FBS(0.5%)",
                            "frac.air.adi.FBS(0.5%)", "frac.dis.adi.FBS(0%)",
                            "frac.prot-cell.adi.FBS(0%)", "frac.lip-cell.adi.FBS(0%)",
                            "frac.air.adi.FBS(0%)")

# Export results
write.csv(final.resultAdi, file = "CabinetAdiposeV2.csv")
# To make plots, jump to section "Plots with adipose cells"

# Aroclor 1016 (top 12) with adipose cells --------------------------------
# Read data.xlsx
d.A1016 <- data.frame(read_excel("DataAdi.xlsx", sheet = "Aroclor1016",
                                 col_names = TRUE, col_types = NULL))

# Name parameters
congener <- d.A1016$congener
logKa.w <- d.A1016$logKair.water
dUaw <- d.A1016$dUaw
logKlip.w <- d.A1016$logKlipid.water
logKpro.w <- d.A1016$logKprotein.water
logKalb.w <- d.A1016$logKalbumin.water
R <- d.A1016$R
tst <- d.A1016$tst
texp <- d.A1016$texp

num.congener <- length(congener)
resultAdi <- NULL
for (i in 1:num.congener) {
  resultAdi <- rbind(resultAdi, fractionAdi(logKa.w[i], dUaw[i], logKlip.w[i],
                                            logKpro.w[i], logKalb.w[i], R[i],
                                            tst[i], texp[i]))
}

final.resultAdi <- data.frame(congener, resultAdi)
names(final.resultAdi) <- c("congener", "frac.dis.adi.FBS(15%)",
                            "frac.alb.adi.FBS(15%)", "frac.prot.adi.FBS(15%)",
                            "frac.prot-cell.adi.FBS(15%)", "frac.lip-cell.adi.FBS(15%)",
                            "frac.air.adi.FBS(15%)", "frac.dis.adi.FBS(0.5%)",
                            "frac.alb.adi.FBS(0.5%)", "frac.prot.adi.FBS(0.5%)",
                            "frac.prot-cell.adi.FBS(0.5%)", "frac.lip-cell.adi.FBS(0.5%)",
                            "frac.air.adi.FBS(0.5%)", "frac.dis.adi.FBS(0%)",
                            "frac.prot-cell.adi.FBS(0%)", "frac.lip-cell.adi.FBS(0%)",
                            "frac.air.adi.FBS(0%)")

# Export results
write.csv(final.resultAdi, file = "A1016AdiposeV2.csv")
# To make plots, jump to section "Plots with adipose cells"

# Aroclor 1254 (top 12) with adipose cells --------------------------------
# Read data.xlsx
d.A1254 <- data.frame(read_excel("DataAdi.xlsx", sheet = "Aroclor1254",
                                 col_names = TRUE, col_types = NULL))

# Name parameters
congener <- d.A1254$congener
logKa.w <- d.A1254$logKair.water
dUaw <- d.A1254$dUaw
logKlip.w <- d.A1254$logKlipid.water
logKpro.w <- d.A1254$logKprotein.water
logKalb.w <- d.A1254$logKalbumin.water
R <- d.A1254$R
tst <- d.A1254$tst
texp <- d.A1254$texp

num.congener <- length(congener)
resultAdi <- NULL
for (i in 1:num.congener) {
  resultAdi <- rbind(resultAdi, fractionAdi(logKa.w[i], dUaw[i], logKlip.w[i],
                                            logKpro.w[i], logKalb.w[i], R[i],
                                            tst[i], texp[i]))
}

final.resultAdi <- data.frame(congener, resultAdi)
names(final.resultAdi) <- c("congener", "frac.dis.adi.FBS(15%)",
                            "frac.alb.adi.FBS(15%)", "frac.prot.adi.FBS(15%)",
                            "frac.prot-cell.adi.FBS(15%)", "frac.lip-cell.adi.FBS(15%)",
                            "frac.air.adi.FBS(15%)", "frac.dis.adi.FBS(0.5%)",
                            "frac.alb.adi.FBS(0.5%)", "frac.prot.adi.FBS(0.5%)",
                            "frac.prot-cell.adi.FBS(0.5%)", "frac.lip-cell.adi.FBS(0.5%)",
                            "frac.air.adi.FBS(0.5%)", "frac.dis.adi.FBS(0%)",
                            "frac.prot-cell.adi.FBS(0%)", "frac.lip-cell.adi.FBS(0%)",
                            "frac.air.adi.FBS(0%)")

# Export results
write.csv(final.resultAdi, file = "A1254AdiposeV2.csv")
# To make plots, jump to section "Plots with adipose cells"

# Plots with adipose cells ------------------------------------------------
# final.resultAdi needs to be changed, depending on the plot
# Create data.frame with needed fractions
# (1) FBS = 15%
p.FBS.h.adi <- final.resultAdi[,!names(final.resultAdi) %in% c("frac.dis.adi.FBS(0.5%)",
                                                               "frac.alb.adi.FBS(0.5%)", "frac.prot.adi.FBS(0.5%)",
                                                               "frac.prot-cell.adi.FBS(0.5%)", "frac.lip-cell.adi.FBS(0.5%)",
                                                               "frac.air.adi.FBS(0.5%)", "frac.dis.adi.FBS(0%)",
                                                               "frac.prot-cell.adi.FBS(0%)", "frac.lip-cell.adi.FBS(0%)",
                                                               "frac.air.adi.FBS(0%)")]

# Transform data.frame p.1 to 3 column data.frame
p.FBS.h.adi <- melt(p.FBS.h.adi, id.var = c("congener"),
                 variable.name = "phase", value.name = "fraction")

# For A1016 and A1254, congeners need to be changed
# Cabinet mixture
p.FBS.h.adi$congener <- factor(p.FBS.h.adi$congener,
                            levels = c('PCB11', 'PCB47', 'PCB51',
                                       'PCB68'))

# A1016
p.FBS.h.adi$congener <- factor(p.FBS.h.adi$congener,
                           levels = c('PCB4', 'PCB8', 'PCB16',
                                      'PCB17', 'PCB18', 'PCB22',
                                      'PCB28', 'PCB31', 'PCB33',
                                      'PCB44', 'PCB49', 'PCB52'))

# A1254
p.FBS.h.adi$congener <- factor(p.FBS.h.adi$congener,
                           levels = c('PCB66', 'PCB70', 'PCB85',
                                      'PCB87', 'PCB97', 'PCB99',
                                      'PCB101', 'PCB105', 'PCB110',
                                      'PCB118', 'PCB138', 'PCB153'))

# Organize fraction to be displayed in plot
p.FBS.h.adi$phase <- factor(p.FBS.h.adi$phase,
                         levels = c('frac.dis.adi.FBS(15%)', 'frac.alb.adi.FBS(15%)',
                                    'frac.prot.adi.FBS(15%)', 'frac.prot-cell.adi.FBS(15%)',
                                    'frac.lip-cell.adi.FBS(15%)', 'frac.air.adi.FBS(15%)'))

# Plot
ggplot(p.FBS.h.adi, aes(x = congener, y = fraction, fill = phase)) + 
  geom_bar(stat = "identity", col = "white") +
  scale_fill_manual(labels = c("Air" , "Albumin-FBS (15%)", "Protein-FBS (15%)",
                               "Medium", "Lipid-cell", "Protein-cell"),
                    values = c("deepskyblue", "lightgrey", "coral4", "red", "blue",
                              "orange")) +
  theme_classic() +
  xlab(expression(bold(""))) +
  ylab(expression("Fraction in well")) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.title.x = element_text(face = "bold", size = 8))


# (2) FBS = 0.5%
p.FBS.l.adi <- final.resultAdi[,!names(final.resultAdi) %in% c("frac.dis.adi.FBS(15%)",
                                                     "frac.alb.adi.FBS(15%)", "frac.prot.adi.FBS(15%)",
                                                     "frac.prot-cell.adi.FBS(15%)", "frac.lip-cell.adi.FBS(15%)",
                                                     "frac.air.adi.FBS(15%)", "frac.dis.adi.FBS(0%)",
                                                     "frac.prot-cell.adi.FBS(0%)", "frac.lip-cell.adi.FBS(0%)",
                                                     "frac.air.adi.FBS(0%)")]

# Transform data.frame p.1 to 3 column data.frame
p.FBS.l.adi <- melt(p.FBS.l.adi, id.var = c("congener"),
                variable.name = "phase", value.name = "fraction")

# For A1016 and A1254, congeners need to be changed
# Cabinet mixture
p.FBS.l.adi$congener <- factor(p.FBS.l.adi$congener,
                           levels = c('PCB11', 'PCB47', 'PCB51',
                                      'PCB68'))

# A1016
p.FBS.l.adi$congener <- factor(p.FBS.l.adi$congener,
                               levels = c('PCB4', 'PCB8', 'PCB16',
                                          'PCB17', 'PCB18', 'PCB22',
                                          'PCB28', 'PCB31', 'PCB33',
                                          'PCB44', 'PCB49', 'PCB52'))

# A1254
p.FBS.l.adi$congener <- factor(p.FBS.l.adi$congener,
                               levels = c('PCB66', 'PCB70', 'PCB85',
                                          'PCB87', 'PCB97', 'PCB99',
                                          'PCB101', 'PCB105', 'PCB110',
                                          'PCB118', 'PCB138', 'PCB153'))

# Organize fraction to be displayed in plot
p.FBS.l.adi$phase <- factor(p.FBS.l.adi$phase,
                        levels = c('frac.dis.adi.FBS(0.5%)', 'frac.alb.adi.FBS(0.5%)',
                                   'frac.prot.adi.FBS(0.5%)', 'frac.prot-cell.adi.FBS(0.5%)',
                                   'frac.lip-cell.adi.FBS(0.5%)', 'frac.air.adi.FBS(0.5%)'))

# Plot
ggplot(p.FBS.l.adi, aes(x = congener, y = fraction, fill = phase)) + 
  geom_bar(stat = "identity", col = "white") +
  scale_fill_manual(labels = c("Air" , "Albumin-FBS (0.5%)", "Protein-FBS (0.5%)",
                               "Medium", "Lipid-cell", "Protein-cell"),
                    values = c("deepskyblue", "lightgrey", "coral4", "red", "blue",
                               "orange")) +
  theme_classic() +
  xlab(expression(bold(""))) +
  ylab(expression("Fraction in well")) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.title.x = element_text(face = "bold", size = 8))

# (3) FBS = 0.0%
p.FBS.0.adi <- final.resultAdi[,!names(final.resultAdi) %in% c("frac.dis.adi.FBS(15%)", "frac.alb.adi.FBS(15%)",
                                                               "frac.prot.adi.FBS(15%)", "frac.prot-cell.adi.FBS(15%)",
                                                               "frac.lip-cell.adi.FBS(15%)", "frac.air.adi.FBS(15%)",
                                                               "frac.dis.adi.FBS(0.5%)", "frac.alb.adi.FBS(0.5%)",
                                                               "frac.prot.adi.FBS(0.5%)", "frac.prot-cell.adi.FBS(0.5%)",
                                                               "frac.lip-cell.adi.FBS(0.5%)", "frac.air.adi.FBS(0.5%)")]

# Transform data.frame p.1 to 3 column data.frame
p.FBS.0.adi <- melt(p.FBS.0.adi, id.var = c("congener"),
                variable.name = "phase", value.name = "fraction")

# For A1016 and A1254, congeners need to be changed
# Cabinet mixture
p.FBS.0.adi$congener <- factor(p.FBS.0.adi$congener,
                           levels = c('PCB11', 'PCB47', 'PCB51',
                                      'PCB68'))

#A1016
p.FBS.0.adi$congener <- factor(p.FBS.0.adi$congener,
                               levels = c('PCB4', 'PCB8', 'PCB16',
                                          'PCB17', 'PCB18', 'PCB22',
                                          'PCB28', 'PCB31', 'PCB33',
                                          'PCB44', 'PCB49', 'PCB52'))

# A1254
p.FBS.0.adi$congener <- factor(p.FBS.0.adi$congener,
                               levels = c('PCB66', 'PCB70', 'PCB85',
                                          'PCB87', 'PCB97', 'PCB99',
                                          'PCB101', 'PCB105', 'PCB110',
                                          'PCB118', 'PCB138', 'PCB153'))

# Organize fraction to be displayed in plot
p.FBS.0.adi$phase <- factor(p.FBS.0.adi$phase,
                        levels = c('frac.air.adi.FBS(0%)', 'frac.dis.adi.FBS(0%)',
                                   'frac.lip-cell.adi.FBS(0%)',
                                   'frac.prot-cell.adi.FBS(0%)'))

# Plot
ggplot(p.FBS.0.adi, aes(x = congener, y = fraction, fill = phase)) + 
  geom_bar(stat = "identity", col = "white") +
  scale_fill_manual(labels = c("air", "medium", "Lipid-cell", "Protein-cell"),
                    values = c("deepskyblue", "red", "blue",
                               "orange")) +
  theme_classic() +
  xlab(expression(bold(""))) +
  ylab(expression("Fraction in well")) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.title.x = element_text(face = "bold", size = 8))

