# Install and upload packages

library(dplyr)
library(ggplot2)
library(ggpubr)
install.packages("BiocManager")

#lefser and SummarizedExperiment packages were installed with BiocManager
BiocManager::install("SummarizedExperiment")
BiocManager::install("lefser")
library(lefser)
library(SummarizedExperiment)

#upload data
data_lefse <- read.csv2("normalized_genus_Fem_C_B_lefse.csv")
View(data_lefse)


data_one <- data_lefse%>%dplyr::select(starts_with("Group"),starts_with("Control"),
                                       starts_with("Binge"))

data_main<- data_one
data_main$Group <- make.unique(as.character(data_main$Group))
data_main<- data_main[, -1]

View(data_main)
View(data_one)

# Function to rename columns based on pattern
rename_columns <- function(df, pattern, new_name) {
  # Find columns that are related to the pattern
  cols_to_rename <- grep(pattern, colnames(df), value = FALSE)
  
  # Update columns' names
  colnames(df)[cols_to_rename] <- new_name
  
  return(df)
}

# Rename columns that start with "Control" to "Control_fem"
data_main <- rename_columns(data_main, "^Control", "Control")

# Rename columns that start with "Binge" to "Binge_fem"
data_main <- rename_columns(data_main, "^Binge", "Binge")

# Check new columns' names
print(colnames(data_main))


duplicated_rows <- data_one$Group[duplicated(data_one$Group)]
print(duplicated_rows)
data_one$Group <- make.unique(as.character(data_one$Group))


#aplly numeric
data_main <- as.data.frame(lapply(data_main, as.numeric))

rownames(data_main) <- data_one$Group
View(data_main)
data_one<- data_one[, -1]


# modification of data

relative_abundances <- apply(data_main, 2, function(x) (x / sum(x)) * 1e6)

coldata <- data.frame(
  SampleName = colnames(relative_abundances),
  Grupo = gl(2, 5, labels = c("Control", "Binge"))
)

coldata$Grupo <- factor(coldata$Grupo, levels = c("Binge", "Control"))


se <- SummarizedExperiment(
  assays =  list(counts = relative_abundances),
  colData = coldata
)


res <- lefser(se, groupCol = "Grupo", lda.threshold = 2.0)

res<-res%>%
  mutate(Grupinho = ifelse(scores< 0.0, "Binge", "Control"))%>%
  mutate(Names = factor(Names, levels = res %>%
                          arrange(Grupinho, scores) %>%
                          pull(Names))) 

C_B_fem  <- ggplot(res, aes(x = Names, y = scores, fill = Grupinho)) +
  geom_bar(stat = "identity", color = "black", size = 0.4, width = 0.8) +
  scale_fill_manual(values = c("Binge" = "red", "Control" = "forestgreen")) +
  ylab("LDA SCORE (log 10)") +
  theme_void() + scale_y_continuous(limits = c(-5.0,5.0))+
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 20, face = "bold", family = "Times", vjust = -2.5),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 18, family = "Times", vjust = -1.5)
  ) + 
  theme(    
    legend.position = "top",
    legend.title = element_blank(),
    legend.key.height = unit(0.07, 'cm'),
    legend.key.width = unit(0.6, 'cm'),
    legend.text = element_text(size = 18, family = "Times"),
    legend.margin=margin(-20,0,2,0),
    legend.box.margin=margin(10,10,10,10),
    plot.margin = unit(c(10,20,20,10), "mm")
  ) +
  geom_text(
    aes(y = 0, label = Names), 
    hjust = ifelse(res$scores < 0, 0, 1),
    nudge_y = ifelse(res$scores < 0, 0.1, -0.1),
    color = "black",
    size = 6,
    family = "Times"
  ) +    
  theme(    
    panel.grid.major.x = element_line(
      color = "grey", linewidth =  0.5, linetype = "dotted"),
    panel.grid.minor.x = element_line(
      color = "grey", linewidth =  0.5, linetype = "dotted")
  )  +
  coord_flip() 

print(C_B_fem)



############################################################################################
############################################################################################
############################################################################################

#Second plot

getwd()
#upload data
data_lefse_2_EB_fem <- read.csv2("normalized_genus_Female_C_versus_EB.csv")
View(data_lefse_2_EB_fem)


data_two <- data_lefse_2_EB_fem%>%dplyr::select(starts_with("Group"),starts_with("Control"),
                                       starts_with("Chronic"))

data_main_2<- data_two
data_main_2$Group <- make.unique(as.character(data_main_2$Group))
data_main_2<- data_main_2[, -1]

View(data_main_2)
View(data_two)


# Rename columns that start with "Control" to "Control_fem"
data_main_2 <- rename_columns(data_main_2, "^Control", "Control")

# Rename columns that start with "Binge" to "Binge_fem"
data_main_2 <- rename_columns(data_main_2, "^Chronic", "Chronic+Binge")

# Check new columns' names
print(colnames(data_main_2))


duplicated_rows_1 <- data_two$Group[duplicated(data_two$Group)]
print(duplicated_rows_1)
data_two$Group <- make.unique(as.character(data_two$Group))


#aplly numeric
data_main_2 <- as.data.frame(lapply(data_main_2, as.numeric))

rownames(data_main_2) <- data_two$Group
View(data_main_2)
data_two<- data_two[, -1]


# modification of data

relative_abundances_EB <- apply(data_main_2, 2, function(x) (x / sum(x)) * 1e6)

coldata_EB <- data.frame(
  SampleName = colnames(relative_abundances_EB),
  Grupo = gl(2, 5, labels = c("Control", "Chronic+Binge"))
)


coldata_EB$Grupo <- factor(coldata_EB$Grupo, levels = c("Chronic+Binge", "Control"))

se_EB <- SummarizedExperiment(
  assays =  list(counts = relative_abundances_EB),
  colData = coldata_EB
)


res_EB <- lefser(se_EB, groupCol = "Grupo", lda.threshold = 2.0)


res_EB<-res_EB%>%
  mutate(Grupinho = ifelse(scores< 0.0, "Chronic+Binge", "Control"))%>%
  mutate(Names = factor(Names, levels = res_EB %>%
                          arrange(Grupinho, scores) %>%
                          pull(Names))) 

EB_fem   <- ggplot(res_EB, aes(x = Names, y = scores, fill = Grupinho)) +
  geom_bar(stat = "identity", color = "black", size = 0.4, width = 0.9) +
  scale_fill_manual(values = c("Chronic+Binge" = "red", "Control" = "forestgreen")) +
  ylab("LDA SCORE (log 10)") +
  theme_void() + scale_y_continuous(limits = c(-5.0,5.0))+
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 20, face = "bold", family = "Times", vjust = -2.5),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 18, family = "Times", vjust = -1.5)
  ) + 
  theme(    
    legend.position = "top",
    legend.title = element_blank(),
    legend.key.height = unit(0.07, 'cm'),
    legend.key.width = unit(0.6, 'cm'),
    legend.text = element_text(size = 18, family = "Times"),
    legend.margin=margin(-20,0,2,0),
    legend.box.margin=margin(10,10,10,10),
    plot.margin = unit(c(10,20,10,10), "mm")
  ) +
  geom_text(
    aes(y = 0, label = Names), 
    hjust = ifelse(res_EB$scores < 0, 0, 1),
    nudge_y = ifelse(res_EB$scores < 0, 0.1, -0.1),
    color = "black",
    size = 6,
    family = "Times"
  ) +    
  theme(    
    panel.grid.major.x = element_line(
      color = "grey", linewidth =  0.5, linetype = "dotted"),
    panel.grid.minor.x = element_line(
      color = "grey", linewidth =  0.5, linetype = "dotted")
  )  +
  coord_flip() 



print(EB_fem)

##############################################################################################################
##############################################################################################################

#Third plot

getwd()
#upload data
data_lefse_Cro_fem <- read.csv2("normalized_genus_female_C_versus_ECronico.csv")
View(data_lefse_Cro_fem)


data_three <- data_lefse_Cro_fem%>%dplyr::select(starts_with("Group"),starts_with("Control"),
                                                starts_with("Chronic"))

data_main_3<- data_three
data_main_3$Group <- make.unique(as.character(data_main_3$Group))
data_main_3 <- data_main_3[, -1]

View(data_main_3)
View(data_three)


# Rename columns that start with "Control" to "Control_fem"
data_main_3 <- rename_columns(data_main_3, "^Control", "Control")

# Rename columns that start with "Binge" to "Binge_fem"
data_main_3 <- rename_columns(data_main_3, "^Chronic", "Chronic")

# Check new columns' names
print(colnames(data_main_3))


duplicated_rows_2 <- data_three$Group[duplicated(data_three$Group)]
print(duplicated_rows_2)
data_three$Group <- make.unique(as.character(data_three$Group))


#aplly numeric
data_main_3 <- as.data.frame(lapply(data_main_3, as.numeric))

rownames(data_main_3) <- data_three$Group
View(data_main_3)
data_three<- data_three[, -1]


# modification of data

relative_abundances_Cro_fem <- apply(data_main_3, 2, function(x) (x / sum(x)) * 1e6)

coldata_Cro <- data.frame(
  SampleName = colnames(relative_abundances_Cro_fem),
  Grupo = gl(2, 5, labels = c("Control", "Chronic"))
)


coldata_Cro$Grupo <- factor(coldata_Cro$Grupo, levels = c("Chronic", "Control"))

se_Cro <- SummarizedExperiment(
  assays =  list(counts = relative_abundances_Cro_fem),
  colData = coldata_Cro
)


res_Cro <- lefser(se_Cro, groupCol = "Grupo", lda.threshold = 2.0)

res_Cro <- lefser(se_Cro, groupCol = "Grupo", lda.threshold = 2.0)%>%
  mutate(Grupinho = ifelse(scores< 0.0, "Chronic", "Control"))%>%
  mutate(Names = factor(Names, levels = res_Cro %>%
                          arrange(Grupinho, scores) %>%
                          pull(Names))) 

Cro_fem <- ggplot(res_Cro, aes(x = Names, y = scores, fill = Grupinho)) +
  geom_bar(stat = "identity", color = "black", size = 0.4, width = 0.8) +
  scale_fill_manual(values = c("Chronic" = "red", "Control" = "forestgreen")) +
  ylab("LDA SCORE (log 10)") +
  theme_void() + scale_y_continuous(limits = c(-5.0,5.0))+
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 20, face = "bold", family = "Times", vjust = -2.5),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 18, family = "Times", vjust = -1.5)
  ) + 
  theme(    
    legend.position = "top",
    legend.title = element_blank(),
    legend.key.height = unit(0.07, 'cm'),
    legend.key.width = unit(0.6, 'cm'),
    legend.text = element_text(size = 18, family = "Times"),
    legend.margin=margin(-20,0,2,0),
    legend.box.margin=margin(10,10,10,10),
    plot.margin = unit(c(10,20,20,10), "mm")
  ) +
  geom_text(
    aes(y = 0, label = Names), 
    hjust = ifelse(res_Cro$scores < 0, 0, 1),
    nudge_y = ifelse(res_Cro$scores < 0, 0.1, -0.1),
    color = "black",
    size = 6,
    family = "Times"
  ) +    
  theme(    
    panel.grid.major.x = element_line(
      color = "grey", linewidth =  0.5, linetype = "dotted"),
    panel.grid.minor.x = element_line(
      color = "grey", linewidth =  0.5, linetype = "dotted")
  )  +
  coord_flip() 

print(Cro_fem)

#############################################################################

#Join the plots 

library(ggpubr)

plot_lefse_fem <-ggarrange(Cro_fem, C_B_fem,EB_fem,  nrow =3, ncol = 1, labels = c("D", "E", "F"), 
                           font.label = list(size = 16, color = "black", face = "bold", family = "Times"))

hi <- annotate_figure(plot_lefse_fem, top = text_grob("FEMALE", 
                                                color = "black", face = "bold", 
                                                size = 24, family = "Times"))


############################################################################################################
############################################################################################################
############################################################################################################


########################################################################

#Male plot



#first plot
getwd()
#upload data
data_lefse_Cro_male <- read.csv2("normalized_genus Male C versus EtCr.csv")
View(data_lefse_Cro_male)


data_four <- data_lefse_Cro_male%>%dplyr::select(starts_with("Group"),starts_with("Control"),
                                                 starts_with("Chronic"))

data_main_4<- data_four
data_main_4$Group <- make.unique(as.character(data_main_4$Group))
data_main_4 <- data_main_4[, -1]

View(data_main_4)
View(data_four)


# Rename columns that start with "Control" to "Control_fem"
data_main_4 <- rename_columns(data_main_4, "^Control", "Control")

# Rename columns that start with "Binge" to "Binge_fem"
data_main_4 <- rename_columns(data_main_4, "^Chronic", "Chronic")

# Check new columns' names
print(colnames(data_main_4))


duplicated_rows_4 <- data_four$Group[duplicated(data_four$Group)]
print(duplicated_rows_4)
data_four$Group <- make.unique(as.character(data_four$Group))


#aplly numeric
data_main_4 <- as.data.frame(lapply(data_main_4, as.numeric))

rownames(data_main_4) <- data_four$Group
View(data_main_4)
data_four <- data_four[, -1]


# modification of data

relative_abundances_Cro_male <- apply(data_main_4, 2, function(x) (x / sum(x)) * 1e6)

coldata_Cro_male <- data.frame(
  SampleName = colnames(relative_abundances_Cro_male),
  Grupo = gl(2, 5, labels = c("Control", "Chronic"))
)


coldata_Cro_male$Grupo <- factor(coldata_Cro_male$Grupo, levels = c("Chronic", "Control"))

se_Cro_male <- SummarizedExperiment(
  assays =  list(counts = relative_abundances_Cro_male),
  colData = coldata_Cro_male
)


res_Cro_male <- lefser(se_Cro_male, groupCol = "Grupo", lda.threshold = 2.0)



Cro_male <-lefserPlot(res_Cro_male,  label.font.size = 6, trim.names = TRUE)

print(Cro_male)


Cro_male <- Cro_male + theme(legend.margin=margin(-20,0,0,0),
                           legend.box.margin=margin(10,10,10,10),
                           axis.text.x = element_text(size = 18, family = "Times" , vjust = -1.5),
                           axis.title.x = element_text(size = 20, face="bold",family = "Times", vjust = -2.5),
                           legend.text=element_text(size=18, family = "Times"),  plot.margin = unit(c(6,30,20,10), "mm")) 
print(Cro_male)


#########################################################################################################################
#########################################################################################################################
#########################################################################################################################

#second plot
getwd()
#upload data
data_lefse_Binge_male <- read.csv2("normalized_genus_Male_C_versus_B.csv")
View(data_lefse_Binge_male)


data_five <- data_lefse_Binge_male%>%dplyr::select(starts_with("Group"),starts_with("Control"),
                                                 starts_with("Binge"))

data_main_5<- data_five
data_main_5$Group <- make.unique(as.character(data_main_5$Group))
data_main_5 <- data_main_5[, -1]

View(data_main_5)
View(data_five)



# Rename columns that start withm "Control" to "Control_fem"
data_main_5 <- rename_columns(data_main_5, "^Control", "Control")

# Rename columns that start with "Binge" to "Binge_fem"
data_main_5 <- rename_columns(data_main_5, "^Binge", "Binge")

# Check new columns' names
print(colnames(data_main_5))


duplicated_rows_5 <- data_five$Group[duplicated(data_five$Group)]
print(duplicated_rows_5)
data_five$Group <- make.unique(as.character(data_five$Group))


#aplly numeric
data_main_5 <- as.data.frame(lapply(data_main_5, as.numeric))

rownames(data_main_5) <- data_five$Group
View(data_main_5)
data_five <- data_five[, -1]


# modification of data

relative_abundances_Binge_male <- apply(data_main_5, 2, function(x) (x / sum(x)) * 1e6)

coldata_Binge_male <- data.frame(
  SampleName = colnames(relative_abundances_Binge_male),
  Grupo = gl(2, 5, labels = c("Control", "Binge"))
)


coldata_Binge_male$Grupo <- factor(coldata_Binge_male$Grupo, levels = c("Binge", "Control"))

se_Binge_male <- SummarizedExperiment(
  assays =  list(counts = relative_abundances_Binge_male),
  colData = coldata_Binge_male
)


res_Binge_male <- lefser(se_Binge_male, groupCol = "Grupo", lda.threshold = 2.0)

Binge_male <-lefserPlot(res_Binge_male,  label.font.size = 6, trim.names = TRUE)

print(Binge_male)


Binge_male <- Binge_male + theme(legend.margin=margin(-20,0,0,0),
                                 legend.box.margin=margin(10,10,10,10),
                             axis.text.x = element_text(size = 18, family = "Times" , vjust = -1.5),
                             axis.title.x = element_text(size = 20, face="bold",family = "Times", vjust = -2.5),
                             legend.text=element_text(size=18, family = "Times"),  plot.margin = unit(c(6,30,20,10), "mm")) 
print(Binge_male)


##########################################################################################################################
##########################################################################################################################

#third plot
getwd()
#upload data
data_lefse_EB_male <- read.csv2("normalized_genus_Male_C_versus_EB.csv")
View(data_lefse_EB_male)


data_six <- data_lefse_EB_male%>%dplyr::select(starts_with("Group"),starts_with("Control"),
                                                   starts_with("Chronic"))

data_main_6<- data_six
data_main_6$Group <- make.unique(as.character(data_main_6$Group))
data_main_6 <- data_main_6[, -1]

View(data_main_6)
View(data_six)


# Rename columns that start with "Control" to "Control_fem"
data_main_6 <- rename_columns(data_main_6, "^Control", "Control")

# Rename columns that start with "Binge" to "Binge_fem"
data_main_6 <- rename_columns(data_main_6, "^Chronic", "Chronic+Binge")

# Check new columns' names
print(colnames(data_main_6))


duplicated_rows_6 <- data_six$Group[duplicated(data_six$Group)]
print(duplicated_rows_6)
data_six$Group <- make.unique(as.character(data_six$Group))


#aplly numeric
data_main_6 <- as.data.frame(lapply(data_main_6, as.numeric))

rownames(data_main_6) <- data_six$Group
View(data_main_6)
data_six <- data_six[, -1]


# modification of data

relative_abundances_EB_male <- apply(data_main_6, 2, function(x) (x / sum(x)) * 1e6)

coldata_EB_male <- data.frame(
  SampleName = colnames(relative_abundances_Binge_male),
  Grupo = gl(2, 5, labels = c("Control", "Chronic+Binge"))
)

coldata_EB_male$Grupo <- factor(coldata_EB_male$Grupo, levels = c("Chronic+Binge", "Control"))


se_EB_male <- SummarizedExperiment(
  assays =  list(counts = relative_abundances_EB_male),
  colData = coldata_EB_male
)



res_EB_male <- lefser(se_EB_male, groupCol = "Grupo", lda.threshold = 2.0)%>%
  mutate(Grupinho = ifelse(scores< 0.0, "Chronic+Binge", "Control"))%>%
  mutate(Names = factor(Names, levels = res_EB_male %>%
                          arrange(Grupinho, scores) %>%
                          pull(Names))) 

EB_male <- ggplot(res_EB_male, aes(x = Names, y = scores, fill = Grupinho)) +
  geom_bar(stat = "identity", color = "black", size = 0.4, width = 0.8) +
  scale_fill_manual(values = c("Chronic+Binge" = "red", "Control" = "forestgreen")) +
  ylab("LDA SCORE (log 10)") +
  theme_void() + scale_y_continuous(limits = c(-5.0,5.0))+
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 20, face = "bold", family = "Times", vjust = -2.5),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 18, family = "Times", vjust = -1.5)
  ) + 
  theme(    
    legend.position = "top",
    legend.title = element_blank(),
    legend.key.height = unit(0.07, 'cm'),
    legend.key.width = unit(0.6, 'cm'),
    legend.text = element_text(size = 18, family = "Times"),
    legend.margin=margin(-20,0,2,0),
    legend.box.margin=margin(10,10,10,10),
    plot.margin = unit(c(10,30,10,10), "mm")
  ) +
  geom_text(
    aes(y = 0, label = Names), 
    hjust = ifelse(res_EB_male$scores < 0, 0, 1),
    nudge_y = ifelse(res_EB_male$scores < 0, 0.1, -0.1),
    color = "black",
    size = 6,
    family = "Times"
  ) +    
  theme(    
    panel.grid.major.x = element_line(
      color = "grey", linewidth =  0.5, linetype = "dotted"),
    panel.grid.minor.x = element_line(
      color = "grey", linewidth =  0.5, linetype = "dotted")
  )  +
  coord_flip()

print(EB_male)


#Join the plots 

plot_lefse_male <-ggarrange(Cro_male, Binge_male,EB_male, nrow =3, ncol = 1, labels = c("A", "B", "C"), 
                            font.label = list(size = 16, color = "black", face = "bold", family = "Times"))

oi <-annotate_figure(plot_lefse_male, top = text_grob("MALE", 
                                                color = "black", face = "bold", 
                                                size = 24, family = "Times"))

##########################################
all_lefse_plot <-ggarrange(oi,hi)



print(all_lefse_plot)

