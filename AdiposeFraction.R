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

# Cabinet mixture ---------------------------------------------------------
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

  # albumin concentration from FSB
  C.alb.h <- 5.145/1000 # kg/L (15%)
  C.alb.l <- 0.25/1000 # kg/L (0.5%)
  dalb <- 1 # kg/L ask!
  C.alb.h <- C.alb.h/dalb # Lalb/Lwater
  C.alb.l <- C.alb.l/dalb # Lalb/Lwater
  # protein concentration from FSB
  C.prot.med.h <- 5.145/1000 # kg/L
  dprot <- 1.43 # kg/L ask! ref: https://pubmed.ncbi.nlm.nih.gov/10930825/
  C.prot.med.h <- C.prot.med.h/dprot # Lprot/Lwater
  C.prot.med.l <- 0.25/1000 # kg/L
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

num.congener <- length(congener)
result <- NULL
for (i in 1:num.congener) {
  result <- rbind(result, fraction(logKa.w[i], dUaw[i],
                                   logKpro.w[i], logKalb.w[i], R[i],
                                   tst[i], texp[i]))
}

final.result <- data.frame(congener, result)
names(final.result) <- c("congener", "frac.dis.alb.h", "frac.alb.alb.h",
                         "frac.prot.alb.h", "frac.air.alb.h",
                         "frac.dis.alb.l", "frac.alb.alb.l",
                         "frac.prot.alb.l", "frac.air.alb.l",
                         "frac.dis.alb.0", "frac.air.alb.0")

# Plots
# create data.frame with needed fractions
# (1) FSB = 10%
p.FSB.h <- final.result[,!names(final.result) %in% c("frac.dis.alb.l", "frac.alb.alb.l",
                                                     "frac.prot.alb.l", "frac.air.alb.l",
                                                     "frac.dis.alb.0", "frac.air.alb.0")]

# Transform data.frame p.1 to 3 column data.frame
p.FSB.h <- melt(p.FSB.h, id.var = c("congener"),
                variable.name = "phase", value.name = "fraction")

# Name the compounds
p.FSB.h$congener <- factor(p.FSB.h$congener,
                           levels = c('PCB11', 'PCB47', 'PCB51',
                                      'PCB68'))

# Organize fraction to be displayed in plot
p.FSB.h$phase <- factor(p.FSB.h$phase,
                        levels = c('frac.air.alb.h', 'frac.alb.alb.h',
                                   'frac.prot.alb.h', 'frac.dis.alb.h'))

ggplot(p.FSB.h, aes(x = congener, y = fraction, fill = phase)) + 
  geom_bar(stat = "identity", col = "white") +
  scale_fill_manual(labels = c("air" , "FSB-albumin (10%)", "FSB-protein (10%)",
                               "medium"),
                    values = c("deepskyblue", "lightgrey", "coral4", "red")) +
  theme_classic() +
  xlab(expression(bold(""))) +
  ylab(expression("Fraction in well")) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.title.x = element_text(face = "bold", size = 8))

# (2) FSB = 0.5%
p.FSB.l <- final.result[,!names(final.result) %in% c("frac.dis.alb.h", "frac.alb.alb.h",
                                                     "frac.prot.alb.h", "frac.air.alb.h",
                                                     "frac.dis.alb.0", "frac.air.alb.0")]

# Transform data.frame p.1 to 3 column data.frame
p.FSB.l <- melt(p.FSB.l, id.var = c("congener"),
                variable.name = "phase", value.name = "fraction")

# Name the compounds
p.FSB.l$congener <- factor(p.FSB.l$congener,
                           levels = c('PCB11', 'PCB47', 'PCB51',
                                      'PCB68'))

# Organize fraction to be displayed in plot
p.FSB.l$phase <- factor(p.FSB.l$phase,
                        levels = c('frac.air.alb.l', 'frac.alb.alb.l',
                                   'frac.prot.alb.l', 'frac.dis.alb.l'))

ggplot(p.FSB.l, aes(x = congener, y = fraction, fill = phase)) + 
  geom_bar(stat = "identity", col = "white") +
  scale_fill_manual(labels = c("air" , "FSB-albumin (0.5%)", "FSB-protein (0.5%)",
                               "medium"),
                    values = c("deepskyblue", "lightgrey", "coral4", "red")) +
  theme_classic() +
  xlab(expression(bold(""))) +
  ylab(expression("Fraction in well")) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.title.x = element_text(face = "bold", size = 8))

# (3) FSB = 0.0%
p.FSB.0 <- final.result[,!names(final.result) %in% c("frac.dis.alb.h", "frac.alb.alb.h",
                                                     "frac.prot.alb.h", "frac.air.alb.h",
                                                     "frac.dis.alb.l", "frac.alb.alb.l",
                                                     "frac.prot.alb.l", "frac.air.alb.l")]

# Transform data.frame p.1 to 3 column data.frame
p.FSB.0 <- melt(p.FSB.0, id.var = c("congener"),
                variable.name = "phase", value.name = "fraction")

# Name the compounds
p.FSB.0$congener <- factor(p.FSB.0$congener,
                           levels = c('PCB11', 'PCB47', 'PCB51',
                                      'PCB68'))

# Organize fraction to be displayed in plot
p.FSB.0$phase <- factor(p.FSB.0$phase,
                        levels = c('frac.air.alb.0', 'frac.dis.alb.0'))

ggplot(p.FSB.0, aes(x = congener, y = fraction, fill = phase)) + 
  geom_bar(stat = "identity", col = "white") +
  scale_fill_manual(labels = c("air" , "medium"),
                    values = c("deepskyblue", "red")) +
  theme_classic() +
  xlab(expression(bold(""))) +
  ylab(expression("Fraction in well")) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.title.x = element_text(face = "bold", size = 8))

# Fraction with adipose cells -----------------------------------

fractionAdi = function(logKa.w, dUaw, logKlip.w, logKpro.w, logKalb.w,
                    R, tst, texp) {
  
  # Concentrations need to be in volume (L), not mass
  Adi <- 10^-7 # kg/well
  C.adi <- Adi/Vm # Concentration of cell per well #kg/L
  C.lip.adi <- C.adi*0.6 # lipid content kg/L
  dlip <- 0.905 # kg/L ask! ref: https://pubmed.ncbi.nlm.nih.gov/8148928/
  C.lip.adi <- C.lip.adi/dlip # L/L
  C.prot.adi <- C.adi*0.4 # protein content kg/L
  dprot <- 1.43 # kg/L ask! ref: https://pubmed.ncbi.nlm.nih.gov/10930825/
  C.prot.adi <- C.prot.adi/dprot # L/L
  C.water.adi <- C.adi*0.004 # water content kg/L
  dwater <- 0.99127 # kg/L at 37 C
  C.water.adi <- C.water.adi/dwater # L/L
  V.water.adi <- C.water.adi*Adi/10^6 # L water inside cell
  # albumin concentration from FSB
  C.alb.h <- 4.9/1000 # kg/L
  C.alb.l <- 0.25/1000 # kg/L
  dalb <- 1 # kg/L ask!
  C.alb.h <- C.alb.h/dalb # Lalb/Lwater
  C.alb.l <- C.alb.l/dalb # Lalb/Lwater
  # protein concentration from FSB
  C.prot.med.h <- 4.9/1000 # kg/L
  C.prot.med.h <- C.prot.med.h/dprot # Lprot/Lwater
  C.prot.med.l <- 0.25/1000 # kg/L
  C.prot.med.l <- C.prot.med.l/dprot # Lprot/Lwater
  
  # Temperature correction for Kaw
  Ka.w.t <- 10^(logKa.w)*exp(-dUaw/R*(1/texp-1/tst)) # Ka.w corrected by water and air temps
  logKa.w.t <- log10(Ka.w.t)
  
  # Fraction calculation
  # High FBS (10%)
  den.alb.h.adi <- 1 + 10^(logKalb.w)*C.alb.h + 10^(logKpro.w)*C.prot.med.h +
    10^(logKlip.w)*C.lip.adi + 10^(logKpro.w)*C.prot.adi + 
    10^(logKa.w.t)*Va/(Vm-V.water.adi)
  f.dis.alb.h.adi <- 1/den.alb.h.adi # freely dissolved fraction
  f.alb.alb.h.adi <- 10^(logKalb.w)*C.alb.h/den.alb.h.adi # albumin fraction from FSB
  f.prot.alb.h.adi <- 10^(logKpro.w)*C.prot.med.h/den.alb.h.adi # protein fraction from FSB
  f.prot.adi.alb.h.adi <- 10^(logKpro.w)*C.prot.adi/den.alb.h.adi # protein fraction from adipose cells
  f.lip.adi.alb.h.adi <- 10^(logKlip.w)*C.lip.adi/den.alb.h.adi # lipid fraction from adipose cells
  f.air.alb.h.adi <- 10^(logKa.w.t)*Va/Vm/den.alb.h.adi # air fraction
  
  # Low FBS (0.5%)
  den.alb.l.adi <- 1 + 10^(logKalb.w)*C.alb.l + 10^(logKpro.w)*C.prot.med.l +
    10^(logKlip.w)*C.lip.adi + 10^(logKpro.w)*C.prot.adi + 
    10^(logKa.w.t)*Va/(Vm-V.water.adi)
  f.dis.alb.l.adi <- 1/den.alb.l.adi # freely dissolved fraction
  f.alb.alb.l.adi <- 10^(logKalb.w)*C.alb.l/den.alb.l.adi # albumin fraction from FSB
  f.prot.alb.l.adi <- 10^(logKpro.w)*C.prot.med.l/den.alb.l.adi # protein fraction from FSB
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

num.congener <- length(congener)
resultAdi <- NULL
for (i in 1:num.congener) {
  resultAdi <- rbind(resultAdi, fractionAdi(logKa.w[i], dUaw[i], logKlip.w[i],
                                   logKpro.w[i], logKalb.w[i], R[i],
                                   tst[i], texp[i]))
}

final.resultAdi <- data.frame(congener, resultAdi)
names(final.resultAdi) <- c("congener", "frac.dis.alb.h.adi",
                            "frac.alb.alb.h.adi", "frac.prot.alb.h.adi",
                            "frac.prot.adi.alb.h.adi", "frac.lip.adi.alb.h.adi",
                            "frac.air.alb.h.adi", "frac.dis.alb.l.adi",
                            "frac.alb.alb.l.adi", "frac.prot.alb.l.adi",
                            "frac.prot.adi.alb.l.adi", "frac.lip.adi.alb.l.adi",
                            "frac.air.alb.l.adi", "frac.dis.alb.0.adi",
                            "frac.prot.adi.alb.0.adi", "frac.lip.adi.alb.0.adi",
                            "frac.air.alb.0.adi")

# Plots
# create data.frame with needed fractions
# (1) FSB = 10%
p.FSB.h.adi <- final.resultAdi[,!names(final.resultAdi) %in% c("frac.dis.alb.l.adi",
                                                               "frac.alb.alb.l.adi", "frac.prot.alb.l.adi",
                                                               "frac.prot.adi.alb.l.adi", "frac.lip.adi.alb.l.adi",
                                                               "frac.air.alb.l.adi", "frac.dis.alb.0.adi",
                                                               "frac.prot.adi.alb.0.adi", "frac.lip.adi.alb.0.adi",
                                                               "frac.air.alb.0.adi")]

# Transform data.frame p.1 to 3 column data.frame
p.FSB.h.adi <- melt(p.FSB.h.adi, id.var = c("congener"),
                 variable.name = "phase", value.name = "fraction")

# Name the compounds
p.FSB.h.adi$congener <- factor(p.FSB.h.adi$congener,
                            levels = c('PCB11', 'PCB47', 'PCB51',
                                       'PCB68'))

# Organize fraction to be displayed in plot
p.FSB.h.adi$phase <- factor(p.FSB.h.adi$phase,
                         levels = c('frac.air.alb.h.adi', 'frac.alb.alb.h.adi',
                                    'frac.prot.alb.h.adi', 'frac.dis.alb.h.adi',
                                    'frac.lip.adi.alb.h.adi', 'frac.prot.adi.alb.h.adi'))

ggplot(p.FSB.h.adi, aes(x = congener, y = fraction, fill = phase)) + 
  geom_bar(stat = "identity", col = "white") +
  scale_fill_manual(labels = c("air" , "FSB-albumin (10%)", "FSB-protein (10%)",
                               "medium", "lip-adipose", "prot-adipose"),
                    values = c("deepskyblue", "lightgrey", "coral4", "red", "blue",
                              "orange")) +
  theme_classic() +
  xlab(expression(bold(""))) +
  ylab(expression("Fraction in well")) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.title.x = element_text(face = "bold", size = 8))


# (2) FSB = 0.5%
p.FSB.l <- final.result[,!names(final.result) %in% c("frac.dis.alb.h", "frac.alb.alb.h",
                                                     "frac.prot.alb.h", "frac.air.alb.h",
                                                     "frac.dis.alb.0", "frac.air.alb.0")]

# Transform data.frame p.1 to 3 column data.frame
p.FSB.l <- melt(p.FSB.l, id.var = c("congener"),
                variable.name = "phase", value.name = "fraction")

# Name the compounds
p.FSB.l$congener <- factor(p.FSB.l$congener,
                           levels = c('PCB11', 'PCB47', 'PCB51',
                                      'PCB68'))

# Organize fraction to be displayed in plot
p.FSB.l$phase <- factor(p.FSB.l$phase,
                        levels = c('frac.air.alb.l', 'frac.alb.alb.l',
                                   'frac.prot.alb.l', 'frac.dis.alb.l'))

ggplot(p.FSB.l, aes(x = congener, y = fraction, fill = phase)) + 
  geom_bar(stat = "identity", col = "white") +
  scale_fill_manual(labels = c("air" , "FSB-albumin (0.5%)", "FSB-protein (0.5%)",
                               "medium"),
                    values = c("deepskyblue", "lightgrey", "coral4", "red")) +
  theme_classic() +
  xlab(expression(bold(""))) +
  ylab(expression("Fraction in well")) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.title.x = element_text(face = "bold", size = 8))

# (3) FSB = 0.0%
p.FSB.0 <- final.result[,!names(final.result) %in% c("frac.dis.alb.h", "frac.alb.alb.h",
                                                     "frac.prot.alb.h", "frac.air.alb.h",
                                                     "frac.dis.alb.l", "frac.alb.alb.l",
                                                     "frac.prot.alb.l", "frac.air.alb.l")]

# Transform data.frame p.1 to 3 column data.frame
p.FSB.0 <- melt(p.FSB.0, id.var = c("congener"),
                variable.name = "phase", value.name = "fraction")

# Name the compounds
p.FSB.0$congener <- factor(p.FSB.0$congener,
                           levels = c('PCB11', 'PCB47', 'PCB51',
                                      'PCB68'))

# Organize fraction to be displayed in plot
p.FSB.0$phase <- factor(p.FSB.0$phase,
                        levels = c('frac.air.alb.0', 'frac.dis.alb.0'))

ggplot(p.FSB.0, aes(x = congener, y = fraction, fill = phase)) + 
  geom_bar(stat = "identity", col = "white") +
  scale_fill_manual(labels = c("air" , "medium"),
                    values = c("deepskyblue", "red")) +
  theme_classic() +
  xlab(expression(bold(""))) +
  ylab(expression("Fraction in well")) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.title.x = element_text(face = "bold", size = 8))

# Aroclor 1016 (top 12) ---------------------------------------------------
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

# Aroclor 1016 fraction without adipose cells -----------------------------

# Function to calculate fractions
fractionA.1016 = function(logKa.w, dUaw, logKpro.w, logKalb.w, R,
                    tst, texp) {

  # albumin concentration from FSB
  C.alb.h <- 4.9/1000 # kg/L
  C.alb.l <- 0.25/1000 # kg/L
  dalb <- 1 # kg/L ask!
  C.alb.h <- C.alb.h/dalb # Lalb/Lwater
  C.alb.l <- C.alb.l/dalb # Lalb/Lwater
  # protein concentration from FSB
  C.prot.med.h <- 4.9/1000 # kg/L
  dprot <- 1.43 # kg/L ask! ref: https://pubmed.ncbi.nlm.nih.gov/10930825/
  C.prot.med.h <- C.prot.med.h/dprot # Lprot/Lwater
  C.prot.med.l <- 0.25/1000 # kg/L
  C.prot.med.l <- C.prot.med.l/dprot # Lprot/Lwater
  
  # Temperature correction for Kaw
  Ka.w.t <- 10^(logKa.w)*exp(-dUaw/R*(1/texp-1/tst)) # Ka.w corrected by water and air temps
  logKa.w.t <- log10(Ka.w.t)
  
  # Fraction calculation
  # High FBS (10%)
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
  den.alb.0 <- 1  + 10^(logKa.w.t)*Va/Vm
  f.dis.alb.0 <- 1/den.alb.0 # freely dissolved fraction
  f.air.alb.0 <- 10^(logKa.w.t)*Va/den.alb.0/Vm # air fraction
  
  frac <- c(f.dis.alb.h, f.alb.alb.h, f.prot.alb.h, f.air.alb.h,
            f.dis.alb.l, f.alb.alb.l, f.prot.alb.l, f.air.alb.l,
            f.dis.alb.0, f.air.alb.0)
}

num.congener <- length(congener)
resultA.1016 <- NULL
for (i in 1:num.congener) {
  resultA.1016 <- rbind(resultA.1016, fractionA.1016(logKa.w[i], dUaw[i],
                                   logKpro.w[i], logKalb.w[i], R[i],
                                   tst[i], texp[i]))
}

final.resultA.1016 <- data.frame(congener, resultA.1016)
names(final.resultA.1016) <- c("congener", "frac.dis.alb.h.A1016", "frac.alb.alb.h.A1016",
                         "frac.prot.alb.h.A1016", "frac.air.alb.h.A1016",
                         "frac.dis.alb.l.A1016", "frac.alb.alb.l.A1016",
                         "frac.prot.alb.l.A1016", "frac.air.alb.l.A1016",
                         "frac.dis.alb.0.A1016", "frac.air.alb.0.A1016")

# Plots
# create data.frame with needed fractions
# (1) FSB = 10%
p.FSB.h <- final.result[,!names(final.result) %in% c("frac.dis.alb.l", "frac.alb.alb.l",
                                                     "frac.prot.alb.l", "frac.air.alb.l",
                                                     "frac.dis.alb.0", "frac.air.alb.0")]

# Transform data.frame p.1 to 3 column data.frame
p.FSB.h <- melt(p.FSB.h, id.var = c("congener"),
                variable.name = "phase", value.name = "fraction")

# Name the compounds
p.FSB.h$congener <- factor(p.FSB.h$congener,
                           levels = c('PCB11', 'PCB47', 'PCB51',
                                      'PCB68'))

# Organize fraction to be displayed in plot
p.FSB.h$phase <- factor(p.FSB.h$phase,
                        levels = c('frac.air.alb.h', 'frac.alb.alb.h',
                                   'frac.prot.alb.h', 'frac.dis.alb.h'))

ggplot(p.FSB.h, aes(x = congener, y = fraction, fill = phase)) + 
  geom_bar(stat = "identity", col = "white") +
  scale_fill_manual(labels = c("air" , "FSB-albumin (10%)", "FSB-protein (10%)",
                               "medium"),
                    values = c("deepskyblue", "lightgrey", "coral4", "red")) +
  theme_classic() +
  xlab(expression(bold(""))) +
  ylab(expression("Fraction in well")) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.title.x = element_text(face = "bold", size = 8))

# (2) FSB = 0.5%
p.FSB.l <- final.result[,!names(final.result) %in% c("frac.dis.alb.h", "frac.alb.alb.h",
                                                     "frac.prot.alb.h", "frac.air.alb.h",
                                                     "frac.dis.alb.0", "frac.air.alb.0")]

# Transform data.frame p.1 to 3 column data.frame
p.FSB.l <- melt(p.FSB.l, id.var = c("congener"),
                variable.name = "phase", value.name = "fraction")

# Name the compounds
p.FSB.l$congener <- factor(p.FSB.l$congener,
                           levels = c('PCB11', 'PCB47', 'PCB51',
                                      'PCB68'))

# Organize fraction to be displayed in plot
p.FSB.l$phase <- factor(p.FSB.l$phase,
                        levels = c('frac.air.alb.l', 'frac.alb.alb.l',
                                   'frac.prot.alb.l', 'frac.dis.alb.l'))

ggplot(p.FSB.l, aes(x = congener, y = fraction, fill = phase)) + 
  geom_bar(stat = "identity", col = "white") +
  scale_fill_manual(labels = c("air" , "FSB-albumin (0.5%)", "FSB-protein (0.5%)",
                               "medium"),
                    values = c("deepskyblue", "lightgrey", "coral4", "red")) +
  theme_classic() +
  xlab(expression(bold(""))) +
  ylab(expression("Fraction in well")) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.title.x = element_text(face = "bold", size = 8))

# (3) FSB = 0.0%
p.FSB.0 <- final.result[,!names(final.result) %in% c("frac.dis.alb.h", "frac.alb.alb.h",
                                                     "frac.prot.alb.h", "frac.air.alb.h",
                                                     "frac.dis.alb.l", "frac.alb.alb.l",
                                                     "frac.prot.alb.l", "frac.air.alb.l")]

# Transform data.frame p.1 to 3 column data.frame
p.FSB.0 <- melt(p.FSB.0, id.var = c("congener"),
                variable.name = "phase", value.name = "fraction")

# Name the compounds
p.FSB.0$congener <- factor(p.FSB.0$congener,
                           levels = c('PCB11', 'PCB47', 'PCB51',
                                      'PCB68'))

# Organize fraction to be displayed in plot
p.FSB.0$phase <- factor(p.FSB.0$phase,
                        levels = c('frac.air.alb.0', 'frac.dis.alb.0'))

ggplot(p.FSB.0, aes(x = congener, y = fraction, fill = phase)) + 
  geom_bar(stat = "identity", col = "white") +
  scale_fill_manual(labels = c("air" , "medium"),
                    values = c("deepskyblue", "red")) +
  theme_classic() +
  xlab(expression(bold(""))) +
  ylab(expression("Fraction in well")) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.title.x = element_text(face = "bold", size = 8))


