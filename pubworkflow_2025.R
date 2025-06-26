# Script for publication 2025 #
# created by Michelle Lee
# may 2025

# The purpose of this code is to create the network analyses, figures, model outputs of the manuscript sent for review at Diversity and Distributions. The code will do the following:

# 1 - List the libraries used
# 2 - FIGURE 2 (example networks)
# 3 - Network level models
# 4 - Trophic level models
# 4.3 - FIGURE 3
# 4.6 - FIGURE 4
# 5 - Species level models
# 5.3 - FIGURE 5
# 6 - Environmental models


# 1 - List the libraries used =======
require(tidyverse)
require(vegan)
require(reshape2)
require(bipartite)
require(glmmTMB)
require(DHARMa)
require(MuMIn)
require(ggeffects)
require(car)
setwd("./pubworkflowFILES_2025")

# 1.1 - code for network matrix ====
# code adapted from bipartite function frame2webs author: Jochen Fruend
df2intmatrix <- function(dframe, varnames = c("lower", "higher", "freq"), type.out = "list", emptylist = TRUE) {
  if (length(varnames)==3) {
    if (any(is.na(dframe[,varnames[3]]))) warning(paste("NAs in", varnames[3], "converted to 0"))
    webarray <- tapply(dframe[,varnames[3]],dframe[,varnames[1:2]], sum)
  }
  if (length(varnames)==2) webarray <- tapply(rep(1,nrow(dframe)),dframe[,varnames[1:2]], sum)
  webarray[is.na(webarray)] <- 0   # needs to be done when using tapply: unobserved combinations always get a zero, even with na.rm=T
  if (type.out=="array") return(webarray)
  if (type.out=="list") {
    weblist <- list()
    #for (i in dimnames(webarray)[[3]]) weblist[[i]] <- webarray[,,i]
    #if (emptylist) weblist <- lapply(weblist,empty)
    #return(weblist)
  }
}

# 1.2 - FIGURE 1 =======

figdat <- readRDS("explan_data_03may22.rds") %>% 
  mutate(log_area = log(area))

relationships <- ggplot(figdat, aes(x = log_area, y = min_distance, color = prop_no3n)) +
  xlab("Log islet area") +
  ylab("Islet proximity (m)") +
  geom_point(size = 3)+
  scale_color_gradient(low="blue", high="red", name = "Nitrate (ppm)") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15))
relationships


# 2 - FIGURE 2 =======

# 2.0 - interaction data ====
bothyears_full <- read_csv("bothyears_full.csv")

# 2.1 - color scheme ------

intcol <- "#ABA481"

poll_tax <- read_csv("pollinator_taxonomy.csv")
# 9 order of invertebrate
orderinverts_color <- data.frame(
  poll_order = c("Araneae", "Orthoptera", "Blattodae", "Thysanoptera", "Hemiptera", "Hymenoptera", "Coleoptera", "Lepidoptera", "Diptera"),
  color = c("#DFC23F","#FFEA8F","#E0B700","#E0B700","#BFAD5F","#AFA36F","#FFE25F","#B89600","#FFD61F"),
  poll_num = c(1:9))
orderinverts <- data.frame("poll_order" = c("Araneae", "Orthoptera", "Blattodae", "Thysanoptera", "Hemiptera", "Hymenoptera", "Coleoptera", "Lepidoptera", "Diptera")) %>% 
  right_join(.,
             poll_tax,
             by = "poll_order") %>% 
  left_join(., orderinverts_color, by = "poll_order")
# 9 shades of yellow
# fix with a leftjoin

ggplot(orderinverts_color, aes(x = poll_order, fill = poll_order)) +
  scale_fill_manual("Pollinator order", values = orderinverts_color$color) +
  geom_bar()


# take plant list and filter it based on the plants that could be in any of the networks
# the full plant list has 28 species, but the list with interactions shouldn't have more than 13
plant_tax <- read_csv("tnrs_result.csv") %>% 
  select(plant_species, plant_family, plant_order)
plant_tax$plant_species[ plant_tax$plant_species == "Coccolobo uvifera" ] <- "Coccoloba uvifera"

# join to plant phylogeny to order (numbered from order of angiosperm phylogeny)
plant_phylog <- read_csv("plant_phylog.csv") %>% 
  right_join(plant_tax, by = "plant_order")

# unique plant orders (11):
#Laurales"       "Pandanales"     "Arecales"       "Malpighiales"   "Rosales"        "Myrtales"      
#"Malvales"       "Caryophyllales" "Gentianales"    "Boraginales"    "Asterales"


orderplants_color <- data.frame(
  plant_order = c("Laurales","Pandanales","Arecales","Malpighiales","Rosales","Myrtales","Malvales","Caryophyllales","Gentianales","Boraginales","Asterales"),
  color = c("#5B851C", # med
            "#556142", #lightest
            "#6A9B21", #med dark
            "#59772C", # med light
            "#4F7318", #darkest
            "#5B851C", #med
            "#556142", #lightest
            "#6A9B21", #med dark
            "#59772C", # med light
            "#4F7318",
            "#5B851C" #med
  ))
ggplot(orderplants_color, aes(x = plant_order, fill = plant_order)) +
  scale_fill_manual("Plant order", values = orderplants_color$color) +
  geom_bar()
orderplants <- right_join(plant_phylog,
                          orderplants_color,
                          by = "plant_order")


# 2.2 - print cooper -----

# filtered interactions list
cooper <- bothyears_full %>% filter(islet == "Cooper")

# make interaction web based on the filtered interaction list
coopmatrix <- as.data.frame(df2intmatrix(cooper, varnames = c("plant_species", "poll_species"), type.out = 'array'))

# filter colors for pollinators
coop_pollcols <- data.frame(
  poll_species = c(colnames(coopmatrix))
) %>% 
  left_join(orderinverts, by = "poll_species")

# filter colors for plants
coop_plantcols <- data.frame(
  plant_species = c(rownames(coopmatrix))
) %>% 
  left_join(orderplants, by = "plant_species")

# make sequence for pollinators
coop_pollseq <- filter(orderinverts, orderinverts$poll_species %in% cooper$poll_species)

# make sequence for plants
coop_plantseq <- filter(orderplants, orderplants$plant_species %in% cooper$plant_species)

# make sequence argument for both plants and pollinators
by_order <- list(
  seq.high = coop_pollseq$poll_species,
  seq.low = coop_plantseq$plant_species
)


# web with colors and fading
cooperweb <- plotweb(as.data.frame(coopmatrix),
                     method = "normal", empty = TRUE, arrow="no", 
                     col.interaction= adjustcolor(intcol, alpha.f = 0.75), 
                     col.high = adjustcolor(as.character(coop_pollcols$color)), 
                     col.low= adjustcolor(as.character(coop_plantcols$color)),
                     bor.col.interaction =NA, bor.col.high=NA, 
                     bor.col.low=NA, high.lablength = 5, low.lablength = 5,
                     sequence=by_order,
                     low.abun=NULL, low.abun.col="green", 
                     bor.low.abun.col ="black", high.abun=NULL, high.abun.col="red", 
                     bor.high.abun.col="black", text.rot=90, text.high.col="black", 
                     text.low.col="black", adj.high=NULL, adj.low=NULL, plot.axes = FALSE,
                     add=FALSE, y.lim=c(0,2), x.lim=NULL, low.plot=TRUE, 
                     high.plot=TRUE, high.lab.dis = NULL, 
                     low.lab.dis = NULL, abuns.type="additional")



# 2.3 - print kaula -----

kauladf <- bothyears_full %>% filter(islet == "Kaula")

# make interaction web based on the filtered interaction list
kaumatrix <- as.data.frame(df2intmatrix(kauladf, varnames = c("plant_species", "poll_species"), type.out = 'array'))

# filter colors for pollinators
kau_pollcols <- data.frame(
  poll_species = c(colnames(kaumatrix))
) %>% 
  left_join(orderinverts, by = "poll_species")

# filter colors for plants
kau_plantcols <- data.frame(
  plant_species = c(rownames(kaumatrix))
) %>% 
  left_join(orderplants, by = "plant_species")

# make sequence for pollinators
kau_pollseq <- filter(orderinverts, orderinverts$poll_species %in% kauladf$poll_species)

# make sequence for plants
kau_plantseq <- filter(orderplants, orderplants$plant_species %in% kauladf$plant_species)

# make sequence argument for both plants and pollinators
kau.by_order <- list(
  seq.high = kau_pollseq$poll_species,
  seq.low = kau_plantseq$plant_species
)


# web with colors and fading
kauweb <- plotweb(as.data.frame(kaumatrix),
                  method = "normal", empty = TRUE, arrow="no", 
                  col.interaction= adjustcolor(intcol, alpha.f = 0.75), 
                  col.high = adjustcolor(as.character(kau_pollcols$color)), 
                  col.low= adjustcolor(as.character(kau_plantcols$color)),
                  bor.col.interaction =NA, bor.col.high=NA, 
                  bor.col.low=NA, high.lablength = 5, low.lablength = 5,
                  sequence=kau.by_order,
                  low.abun=NULL, low.abun.col="green", 
                  bor.low.abun.col ="black", high.abun=NULL, high.abun.col="red", 
                  bor.high.abun.col="black", text.rot=90, text.high.col="black", 
                  text.low.col="black", adj.high=NULL, adj.low=NULL, plot.axes = FALSE,
                  add=FALSE, y.lim=c(0,2), x.lim=NULL, low.plot=TRUE, 
                  high.plot=TRUE, high.lab.dis = NULL, 
                  low.lab.dis = NULL, abuns.type="additional")


# 2.4 - print portsmouth -----
port <- bothyears_full %>% filter(islet == "Portsmouth")

# make interaction web based on the filtered interaction list
pormatrix <- as.data.frame(df2intmatrix(port, varnames = c("plant_species", "poll_species"), type.out = 'array'))

# filter colors for pollinators
por_pollcols <- data.frame(
  poll_species = c(colnames(pormatrix))
) %>% 
  left_join(orderinverts, by = "poll_species")

# filter colors for plants
por_plantcols <- data.frame(
  plant_species = c(rownames(pormatrix))
) %>% 
  left_join(orderplants, by = "plant_species")

# make sequence for pollinators
por_pollseq <- filter(orderinverts, orderinverts$poll_species %in% port$poll_species)

# make sequence for plants
por_plantseq <- filter(orderplants, orderplants$plant_species %in% port$plant_species)

# make sequence argument for both plants and pollinators
por.by_order <- list(
  seq.high = por_pollseq$poll_species,
  seq.low = por_plantseq$plant_species
)


# web with colors and fading
porweb <- plotweb(as.data.frame(pormatrix), #------------CHANGE
                  method = "normal", empty = TRUE, arrow="no", 
                  col.interaction= adjustcolor(intcol, alpha.f = 0.75), 
                  col.high = adjustcolor(as.character(por_pollcols$color)), # ------- CHANGE
                  col.low= adjustcolor(as.character(por_plantcols$color)), # -------- CHANGE
                  bor.col.interaction =NA, bor.col.high=NA, 
                  bor.col.low=NA, high.lablength = 5, low.lablength = 5,
                  sequence=por.by_order, # ------------------ CHANGE
                  low.abun=NULL, low.abun.col="green", 
                  bor.low.abun.col ="black", high.abun=NULL, high.abun.col="red", 
                  bor.high.abun.col="black", text.rot=90, text.high.col="black", 
                  text.low.col="black", adj.high=NULL, adj.low=NULL, plot.axes = FALSE,
                  add=FALSE, y.lim=c(0,2), x.lim=NULL, low.plot=TRUE, 
                  high.plot=TRUE, high.lab.dis = NULL, 
                  low.lab.dis = NULL, abuns.type="additional")



# 2.5 - print eastern ------
eastern <- bothyears_full %>% filter(islet == "Eastern")

# make interaction web based on the filtered interaction list
easmatrix <- as.data.frame(df2intmatrix(eastern, varnames = c("plant_species", "poll_species"), type.out = 'array'))

# filter colors for pollinators
eas_pollcols <- data.frame(
  poll_species = c(colnames(easmatrix))
) %>% 
  left_join(orderinverts, by = "poll_species")

# filter colors for plants
eas_plantcols <- data.frame(
  plant_species = c(rownames(easmatrix))
) %>% 
  left_join(orderplants, by = "plant_species")

# make sequence for pollinators
eas_pollseq <- filter(orderinverts, orderinverts$poll_species %in% eastern$poll_species)

# make sequence for plants
eas_plantseq <- filter(orderplants, orderplants$plant_species %in% eastern$plant_species)

# make sequence argument for both plants and pollinators
eas.by_order <- list(
  seq.high = eas_pollseq$poll_species,
  seq.low = eas_plantseq$plant_species
)


# web with colors and fading
easweb <- plotweb(as.data.frame(easmatrix),
                  method = "normal", empty = TRUE, arrow="no", 
                  col.interaction= adjustcolor(intcol, alpha.f = 0.75), 
                  col.high = adjustcolor(as.character(eas_pollcols$color)), 
                  col.low= adjustcolor(as.character(eas_plantcols$color)),
                  bor.col.interaction =NA, bor.col.high=NA, 
                  bor.col.low=NA, high.lablength = 5, low.lablength = 5,
                  sequence=eas.by_order,
                  low.abun=NULL, low.abun.col="green", 
                  bor.low.abun.col ="black", high.abun=NULL, high.abun.col="red", 
                  bor.high.abun.col="black", text.rot=90, text.high.col="black", 
                  text.low.col="black", adj.high=NULL, adj.low=NULL, plot.axes = FALSE,
                  add=FALSE, y.lim=c(0,2), x.lim=NULL, low.plot=TRUE, 
                  high.plot=TRUE, high.lab.dis = NULL, 
                  low.lab.dis = NULL, abuns.type="additional")


# 3 - Network models =====
# 3.1 - Data for network + trophic level ======

# combined networks (data for both 2017 and 2019 included in one network), n = 19 islets
# removes south fighter because there are no soil samples from that islet
dat <- readRDS("full_modelinput_may24.rds") %>% 
  # remove environmental variables that Hillary used in my 2010 papers
  dplyr::select(-c(35:39)) %>%
  na.omit() %>% 
  mutate(log_area = log(area))
dat$fc_LL <- as.numeric(dat$fc_LL)

# create file with south fighter that can be used if soil is dropped, n = 20 islets
dat.wSF <- readRDS("full_modelinput_may24.rds") %>% 
  # remove environmental variables that Hillary used in my 2010 papers
  dplyr::select(-c(35:39)) %>% 
  filter(islet != "Bird", islet != "Pollux") %>% 
  mutate(log_area = log(area))
dat.wSF$fc_LL <- as.numeric(dat.wSF$fc_LL)


# 3.2 - H2 -------

h2.m1 <- lm(
  H2 ~ log(area) + min_distance + prop_no3n + log(area):min_distance,
  data = dat,
  family = gaussian,
  na.action = na.fail
)
h2.dredge.1 <- dredge(h2.m1)
# no model great here
# 2nd best model is area


h2.m2 <- lm(
  H2 ~ log(area),
  data = dat.wSF,
  family = gaussian()
)
# n = 21
summary(h2.m2)
# nothing strong relationships here
simulateResiduals(h2.m2, plot=T)
# no significant problems detected

# 3.3 - average links ======

links.m1 <- lm(
  links_per_sp ~ log_area + min_distance + prop_no3n + log_area:min_distance,
  data = dat,
  family = gaussian,
  na.action = na.fail
)
summary(links.m1)
links.dredge.1 <- dredge(links.m1)
# no model great here
# 2nd best model is area



# 4 - Trophic level =====

vars <- c("min_distance",
          "prop_no3n",
          "log_area"
)
dat.stand <- dat %>% mutate_at(vars, ~(scale(.) %>% as.vector))
str(dat.stand)
dat.stand2 <- dat.wSF %>% mutate_at(vars, ~(scale(.) %>% as.vector))


# 4.1 - pollinator functional complementarity ---------------

fcHL.m1 <- lm(
  fc_HL ~ log_area + min_distance + prop_no3n + log_area:min_distance,
  data = dat.stand,
  family = gaussian(),
  na.action = na.fail
)
summary(fcHL.m1)
fcHL.dredge <- dredge(fcHL.m1)
vif(fcHL.m1)
# area*minimum distance
# note: running area versus log area does make a difference in the model output
# 2nd best model is area
# save model output

fcHL.m2 <- lm(
  fc_HL ~ log_area * min_distance,
  data = dat.stand2,
)
# n = 20
summary(fcHL.m2)
#modelsummary(fcHL.m2)
AIC(fcHL.m2)
# strong relationships here, but what do they mean?
# bigger islets -> higher fc
# larger distance -> higher fc
# biggest islets, farther away -> lower fc
simulateResiduals(fcHL.m2, plot=T)
# no significant problems detected


# poll fc ggeffects figure
pred_fcpoll_area <- ggpredict(fcHL.m2, c("log_area", "min_distance [10:300 by=60]"))
plot(pred_fcpoll_area)

# make raw data its own object to map onto ggpredict
raw <- attr(pred_fcpoll_area, "rawdata")


pollfc_predict <- ggplot(pred_fcpoll_area, aes(x = x, y = predicted, colour = group)) +
  geom_smooth(se = F) +
  ylab("Pollinator functional complementarity") +
  xlab("Log islet area") +
  # used colors seven "blocks" away
  scale_color_manual(name = "Proximity (m)",
                     values = c("#BAC2E9",
                                "#ACB5E5",
                                "#9EA9E0",
                                "#7584D3",
                                "#6072CD",
                                "#4B60C6",
                                "#3B50BB",
                                "#2E3E91",
                                "#27367D",
                                "#1A2453")) +
  geom_point(data = raw, mapping = aes(x = x, y = response), color = "black") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.position = "none"
  )
pollfc_predict


# need to create gradient from "#BAC2E9" to "#1A2453" and take the legend from another figure

pollfc_legend <- ggplot(dat.wSF, aes(x = min_distance, y = fc_LL, color = min_distance)) +
  geom_point() +
  scale_color_gradient(low = "#BAC2E9",
                       high = "#1A2453", 
                       name = "Proximity (m)") +
  theme(legend.title=element_text(size=11), 
        legend.text=element_text(size=12))
pollfc_legend



# 4.2 - plant functional complementarity ---------------

fcLL.m1 <- lm(
  fc_LL ~ log_area + min_distance + prop_no3n + log_area:min_distance,
  data = dat,
  family = gaussian(),
  na.action = na.fail
)
fcLL.dredge <- dredge(fcLL.m1)
# area


fcLL.m2 <- lm(
  fc_LL ~ log_area,
  data = dat.wSF,
  family = gaussian()
)
# n = 20

summary(fcLL.m2)
AIC(fcLL.m2)


# strong relationships here, but what do they mean?
# bigger islets -> higher fc
simulateResiduals(fcLL.m2, plot=T)
# no significant problems detected


# plant fc ggeffects figure
pred_fcplant_area <- ggpredict(fcLL.m2, terms = "log_area", ci.lvl = 0.95)
plot(pred_fcplant_area)



# make raw data its own object to map onto ggpredict
raw_plant <- attr(pred_fcplant_area, "rawdata")


plantfc_predict <- ggplot(pred_fcplant_area, aes(x = x, y = predicted)) +
  geom_smooth(colour = "#CCA65C") +
  ylab("Plant functional complementarity") +
  xlab("Log islet area") +
  geom_point(data = raw_plant, mapping = aes(x = x, y = response)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#CCA65C", colour = "NA") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    axis.title.x=element_blank()
  )
plantfc_predict

#ggsave("./plantfc.png", plantfc_predict, width = 4, height = 4, units = "in")


# 4.3 - FIGURE 3 =======

# make cowplot of two fc figures
fc_vert <- plot_grid(pollfc_predict + xlab(NULL),
                     plantfc_predict,
                     ncol = 1, align = "v")
fc_vert





# 4.4 - pollinator links ======

linksHL.m1 <- lm(
  links_HL ~ log_area + min_distance + prop_no3n + log_area:min_distance,
  data = dat,
  family = gaussian,
  na.action = na.fail
)
summary(linksHL.m1)
linksHL.dredge.1 <- dredge(linksHL.m1)

linksHL.m2 <- lm(
  links_HL ~ log_area,
  data = dat.wSF,
  family = gaussian,
  na.action = na.fail
)
summary(linksHL.m2)


# 4.5 - plant links ========

linksLL.m1 <- lm(
  links_LL ~ log_area + min_distance + prop_no3n + log_area:min_distance,
  data = dat,
  family = gaussian,
  na.action = na.fail
)
summary(linksLL.m1)
linksLL.dredge.1 <- dredge(linksLL.m1)

linksLL.m2 <- lm(
  links_LL ~ 1,
  data = dat.wSF,
  family = gaussian,
  na.action = na.fail
)
summary(linksLL.m2)

# 4.6 - FIGURE 4 ======

# plantRichfig <- ggplot(dat, aes(x = log_area, y = plant_rich)) +
#   geom_point(alpha = 0.5, size = 5) +
#   geom_smooth(method = "lm", linewidth = 3, 
#               color = sizecol, fill = sizecol) +
#   ylab("Plant species richness") +
#   xlab("Log islet area") +
#   theme_bw() +
#   theme(
#     axis.text = element_text(size = 14),
#     axis.title = element_text(size = 18)
#   )
# plantRichfig

plantRichfig <- ggplot(dat, aes(x = log_area, y = plant_rich)) +
  geom_point() +
  ylab("Plant species richness") +
  xlab("Log islet area") +
  geom_smooth(method = "lm", color = "#CCA65C", fill = "#CCA65C") +
  theme_bw() +
  theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 15),axis.title.x=element_blank()
    )
ggsave("./plantsprich.png", plantRichfig, width = 4, height = 4, units = "in")

pollRichfig <- ggplot(dat, aes(x = log_area, y = poll_gen_rich)) +
  geom_point() +
  ylab("Pollinator species richness") +
  xlab("Log islet area") +
  geom_smooth(method = "lm", color = "#CCA65C", fill = "#CCA65C") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),axis.title.x=element_blank()
  )
ggsave("./pollsprich.png", pollRichfig, width = 4, height = 4, units = "in")

pollLinksfig <- ggplot(dat, aes(x = log_area, y = links_HL)) +
  geom_point() +
  ylab("Pollinator niche breadth") +
  xlab("Log islet area") +
  geom_smooth(method = "lm", color = "#CCA65C", fill = "#CCA65C") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),axis.title.x=element_blank()
  )
ggsave("./polllinks.png", pollLinksfig, width = 4, height = 4, units = "in")

netLinksfig <- ggplot(dat, aes(x = log_area, y = links_per_sp)) +
  geom_point() +
  ylab("Average links per species") +
  xlab("Log islet area") +
  geom_smooth(method = "lm", color = "#CCA65C", fill = "#CCA65C") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),axis.title.x=element_blank()
  )
ggsave("./netlinks.png", netLinksfig, width = 4, height = 4, units = "in")


extrafig.comb <- ggarrange(netLinksfig + rremove("xlab") + rremove("legend"),
                           pollLinksfig + rremove("xlab") + rremove("legend"),
                           pollRichfig + rremove("xlab") + rremove("legend"),
                           plantRichfig+ rremove("xlab") + rremove("legend"),
                           labels = "AUTO",
                           font.label = list(size = 18),
                           #legend.grob = ggpubr::get_legend(ggmeansbare),
                           ncol = 2,
                           nrow = 2,
                           align = "hv"
                           #legend = "none"
                           
)
# sometimes issues with this line
# use 
extrafig1 <- annotate_figure(extrafig.comb,
                             bottom = grid::textGrob("Log islet area", gp = grid::gpar(cex = 1.5)),
)
extrafig1


# 5 - Species level models =====

# 5.1 - Data =====
datspec1 <- readRDS("speclev01_combinedfull_dbluthpoll.rds")
datspec2 <- readRDS("speclev01_combinedfull_degpoll.rds")
env.dat <- readRDS("explan_data_03may22.rds") %>% 
  dplyr::select(-c(21:25)) %>% 
  mutate(log_area = log(area)) %>% 
  filter(islet != "Pollux", islet != "Bird", islet != "South Fighter")

# take spaces out of column names
dat_specd <- datspec1 %>% 
  rename_all(function(x) gsub(" ", "_", x)) %>% 
  dplyr::select(islet, Chrysosoma_sp, Megachile_fullawayi, Pachodynerus_nasidens, Piletocera_signiferalis, Scholastes_palmyra, Siphunculina_striolata) %>% 
  melt(.) %>% 
  na.omit() %>% 
  left_join(env.dat, by = "islet") %>% 
  rename(d = value, poll_species = variable) %>% 
  mutate(log_area = log(area)) %>% 
  filter(islet != "Pollux", islet != "Bird", islet != "South Fighter")

dat_specdeg <- datspec2 %>% 
  rename_all(function(x) gsub(" ", "_", x)) %>% 
  dplyr::select(islet, Chrysosoma_sp, Megachile_fullawayi, Pachodynerus_nasidens, Piletocera_signiferalis, Scholastes_palmyra, Siphunculina_striolata) %>% 
  melt(.) %>% 
  na.omit() %>% 
  left_join(env.dat, by = "islet") %>% 
  rename(deg = value, poll_species = variable) %>% 
  mutate(log_area = log(area)) %>% 
  filter(islet != "Pollux", islet != "Bird", islet != "South Fighter")


# 5.2 - Pollinator  model ------

m1.d <- lm(d ~ (log_area * min_distance * prop_no3n):poll_species,
           data = dat_specd,
           na.action = na.fail)
summary(m1.d)
m1sum <- summary(m1.d)

m1.2.d <- lm(d ~ (log_area + min_distance + prop_no3n):poll_species,
             data = dat_specd,
             na.action = na.fail)
summary(m1.2.d)

m1.dredge1 <- dredge(m1.d)
m1.2.dredge1 <- dredge(m1.2.d)
# from the interaction dredge: nothing, then islet size and no3n
# from the no interaction dredge: nothing, then islet size and no3n

m2.deg <- lm(deg ~ log_area:poll_species,
             data = dat_specdeg)
summary(m2.deg)

polldegmeans <- ggemmeans(m2.deg, terms = c("log_area", "poll_species"))
rawpolldeg <- attr(polldegmeans, "rawdata")

plotspecdeg <- ggplot(polldegmeans, aes(x = x, y = predicted, fill = group, color = group)) +
  ylab("Pollinator species interactions") +
  xlab("Log islet area") +
  geom_jitter(data = rawpolldeg, 
             mapping = aes(x = x, y = response), 
             alpha = 0.8, size = 1) +
  #geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill = group), alpha=0.2, color = NA) +
  geom_smooth(
    linewidth = 1
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18)
  )
plotspecdeg





# 5.4 - Plant data =====

# combined
datspec1 <- readRDS("speclevFLOWER_combinedfull_dbluthpoll.rds")
datspec2 <- readRDS("speclevFLOWER_combinedfull_degpoll.rds")
env.dat <- readRDS("explan_data_03may22.rds") %>% 
  dplyr::select(-c(21:25)) %>% 
  mutate(log_area = log(area)) %>% 
  filter(islet != "Pollux", islet != "Bird", islet != "South Fighter")

# take spaces out of column names
dat_specd <- datspec1 %>% 
  rename_all(function(x) gsub(" ", "_", x)) %>% 
  select(islet, Heliotropium_foertherianum, Scaevola_taccada, Cocos_nucifera) %>% 
  melt(.) %>% 
  na.omit() %>% 
  left_join(env.dat, by = "islet") %>% 
  rename(d = value, plant_species = variable) %>% 
  filter(islet != "Pollux", islet != "Bird", islet != "South Fighter")

dat_specdeg <- datspec2 %>% 
  rename_all(function(x) gsub(" ", "_", x))%>% 
  select(islet, Heliotropium_foertherianum, Scaevola_taccada, Cocos_nucifera) %>% 
  melt(.) %>% 
  na.omit() %>% 
  left_join(env.dat, by = "islet") %>% 
  rename(deg = value, plant_species = variable) %>% 
  filter(islet != "Pollux", islet != "Bird", islet != "South Fighter")


# 5.5 - Plant model =====

# Plant species specialization
m1.d <- lm(d ~ (log_area * min_distance * prop_no3n):plant_species,
           data = dat_specd,
           na.action = na.fail)
summary(m1.d)
m1sum <- summary(m1.d)

m1.2.d <- lm(d ~ (log_area + min_distance + prop_no3n):plant_species,
             data = dat_specd,
             na.action = na.fail)
summary(m1.2.d)

m2.d <- lm(d ~ (log_area * min_distance * prop_no3n),
           data = dat_specd,
           na.action = na.fail)
summary(m2.d)
m2.dredge <- dredge(m2.d)
m2.2d <- lm(d ~ min_distance,
            data = dat_specd,
            na.action = na.fail)
summary(m2.2d)


m1.dredge1 <- dredge(m1.d)
m1.2.dredge1 <- dredge(m1.2.d)
# from the interaction dredge: nothing, then minimum distance
# from the no interaction dredge: nothing, then minimum distance

# Plant species degree
m1.deg <- lm(deg ~ (log_area + 
                      min_distance + 
                      prop_no3n + 
                      log_area:min_distance):plant_species,
             data = dat_specdeg, 
             na.action  = na.fail)
summary(m1.deg)
m1.deg.dredge <- dredge(m1.deg)


m2.2deg <- lm(deg ~ (log_area):plant_species,
              data = dat_specdeg,
              na.action = na.fail)
summary(m2.2deg)





# 5.6 - FIGURE 5 =====

# 1 for plants
summary(m2.2deg)
plantdegmeans <- ggemmeans(m2.2deg, terms = c("log_area", "plant_species"))
rawplantdeg <- attr(plantdegmeans, "rawdata")

plotplantspecdeg <- ggplot(plantdegmeans, aes(x = x, y = predicted, fill = group, color = group)) +
  ylab("Plant species interactions") +
  xlab("Log islet area") +
  geom_point(data = rawplantdeg, 
             mapping = aes(x = x, y = response), 
             alpha = 0.2, size = 3) +
  geom_smooth(
    linewidth = 1
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    legend.text = element_text(face ="italic"),axis.title.x=element_blank(),
    legend.position = "none"
  )  +
  scale_fill_discrete("Species", labels = ~ gsub("_", " ", .x)) +
  scale_color_discrete("Species", labels = ~ gsub("_", " ", .x))
plotplantspecdeg
#ggsave("./plantspecdeg.png", plotplantspecdeg, width = 4, height = 4, units = "in")

# 1 for pollinators
polldegmeans <- ggemmeans(m2.deg, terms = c("log_area", "poll_species"))
rawpolldeg <- attr(polldegmeans, "rawdata")

plotspecdeg <- ggplot(polldegmeans, aes(x = x, y = predicted, fill = group, color = group)) +
  ylab("Pollinator species interactions") +
  xlab("Log islet area") +
  geom_point(data = rawpolldeg, 
             mapping = aes(x = x, y = response), 
             alpha = 0.2, size = 3) +
  geom_smooth(
    linewidth = 1
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    axis.title.x=element_blank(),
    legend.text = element_text(face ="italic"),
    legend.position = "none"
  )  +
  scale_fill_discrete("Species", labels = ~ gsub("_", " ", .x)) +
  scale_color_discrete("Species", labels = ~ gsub("_", " ", .x))
plotspecdeg
#ggsave("./pollspecdeg.png", plotspecdeg, width = 4, height = 4, units = "in")


extrafig.comb_specdeg <- ggarrange(plotspecdeg + rremove("xlab"),
                                   plotplantspecdeg + rremove("xlab"),
                                   labels = "AUTO",
                                   font.label = list(size = 18),
                                   #legend.grob = ggpubr::get_legend(ggmeansbare),
                                   ncol = 2,
                                   nrow = 1,
                                   align = "hv"
                                   #legend = "none"
                                   
)
# sometimes issues with this line
# use 
extrafig2 <- annotate_figure(extrafig.comb_specdeg,
                             bottom = grid::textGrob("Log islet area", gp = grid::gpar(cex = 1.5)),
)
extrafig2






# 6 - Environmental variables on species richness =====
# 6.1 - Plant richness ====

plant.div.m1 <- lm(plant_rich ~ log_area + prop_no3n + min_distance, dat)
summary(plant.div.m1)
plant.div.m2 <- lm(plant_rich ~ log_area + prop_no3n + min_distance + log_area:min_distance, 
                   df,
                   na.action = na.fail)
summary(plant.div.m2)
dredge.plantdiv.m1 <- dredge(plant.div.m2)
# log area best predictor here


plantm1 <- lm(plant_rich ~ log_area, dat)
summary(plantm1)


# 6.2 Pollinator richness ===== 

poll.div.m1 <- lm(poll_gen_rich ~ log_area + prop_no3n + min_distance, dat)
summary(poll.div.m1)
poll.div.m2 <- lm(poll_gen_rich ~ log_area + prop_no3n + min_distance + log_area:min_distance, 
                  df,
                  na.action = na.fail)
summary(poll.div.m2)
dredge.polldiv.m1 <- dredge(poll.div.m2)
# log area best predictor here


# 6.3 - Network size ======

dim.div.m1 <- lm(dimensions ~ log_area + prop_no3n + min_distance, dat)
summary(dim.div.m1)
dim.div.m2 <- lm(dimensions ~ log_area * prop_no3n * min_distance, 
                 dat,
                 na.action = na.fail)
summary(dim.div.m2)
dredge.dim.div.m1 <- dredge(dim.div.m2)
# log area best predictor here, but closely followed by log area + no3n + logarea:no3n

dim.div.m3 <- lm(dimensions ~ log_area + prop_no3n + log_area:prop_no3n, dat)
summary(dim.div.m3)
