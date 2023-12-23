#Gonçalo Fonseca Martins Garcia Costa
#Tese: Bioinformatics in Complex Biological Processes: Data Integration techniques to study Translation Regulation
#Universidade de Aveiro | 2023 


############    Packages    ###############

library(pheatmap)
library(readxl)

set.seed(11)

library(tidyr)
library(dplyr)
library(broom)
library(ez)
library(FSA)
library(dunn.test)


#################################################
#                                               #
#                                               #
#             Criação dos heatmaps              #
#                                               #
#                                               #
#################################################


############ RNA-Seq data todos os codoes ##################

# input dos dados
tab1 = read_excel("C:/Users/gonca/Downloads/Resultados_tese/RSCU_rnaseq.xlsx", sheet = "Matrix_R", col_names = T, range = cell_cols("B:BJ") )
tab1 = as.matrix(tab1)

# estabelecer nomes das condições
row.names(tab1) = c("2 hpi - up DEGs", "4 hpi - up DEGs","8 hpi - up DEGs", "4 hpi - down DEGs", "8 hpi - down DEGs")
#View(tab1)

# transpor para ficar com os codoes nas linhas
tab1_inverted = t(tab1)
tab1_inverted
#View(tab1_inverted)

# criação e personalização do heatmap
pheatmap(tab1_inverted, cluster_cols = F, cluster_rows = T,angle_col = "0", display_numbers = F,fontsize = 10, fontsize_col = 13)


############ RNA-Seq data so para os codoes que codificam a leucina ##################

# input dos dados
tab1 = read_excel("C:/Users/gonca/Downloads/Resultados_tese/RSCU_rnaseq.xlsx", sheet = "Leucine_Matrix_R", col_names = T, range = cell_cols("B:G") )
tab1 = as.matrix(tab1)

# estabelecer nomes das condições
row.names(tab1) = c("2 hpi - up DEGs", "4 hpi - up DEGs","8 hpi - up DEGs", "4 hpi - down DEGs", "8 hpi - down DEGs")
#View(tab1)

# transpor para ficar com os codoes nas linhas
tab1_inverted = t(tab1)
tab1_inverted
#View(tab1_inverted)

# criação e personalização do heatmap
pheatmap(tab1_inverted, cluster_cols = F, cluster_rows = T,angle_col = "0",fontsize = 10, fontsize_col = 10)#, display_numbers = T)



############ Ribo-Seq data todos os codoes ##################

# input dos dados
tab2 = read_excel("C:/Users/gonca/Downloads/Resultados_tese/RSCU_Riboseq.xlsx", sheet = "Matrix_R", col_names = T, range = cell_cols("B:BJ") )
tab2 = as.matrix(tab2)

# estabelecer nomes das condições
row.names(tab2) = c("2 hpi - up DEGs", "4 hpi - up DEGs","8 hpi - up DEGs", "4 hpi - down DEGs", "8 hpi - down DEGs")
#View(tab2)

# transpor para ficar com os codoes nas linhas
tab2_inverted = t(tab2)
tab2_inverted
#View(tab2_inverted)

# criação e personalização do heatmap
pheatmap(tab2_inverted, cluster_cols = F,angle_col = "0",fontsize = 10, fontsize_col = 13, display_numbers = F)



############ Ribo-Seq data so para os codoes que codificam a leucina ##################

# input dos dados
tab2 = read_excel("C:/Users/gonca/Downloads/Resultados_tese/RSCU_Riboseq.xlsx", sheet = "Leucine_Matrix_R", col_names = T, range = cell_cols("B:G") )
tab2 = as.matrix(tab2)

# estabelecer nomes das condições
row.names(tab2) = c("2 hpi - up DEGs", "4 hpi - up DEGs","8 hpi - up DEGs", "4 hpi - down DEGs", "8 hpi - down DEGs")
#View(tab2)

# transpor para ficar com os codoes nas linhas
tab2_inverted = t(tab2)
tab2_inverted
#View(tab2_inverted)

# criação e personalização do heatmap
pheatmap(tab2_inverted, cluster_cols = F, cluster_rows = T,angle_col = "0",fontsize = 10, fontsize_col = 10)





#################################################
#                                               #
#                                               #
#    Teste de independencia do Qui-quadrado     #
#                 RNA-Seq data                  #
#                                               #
#################################################




### calcular RSCUs
calculate_RSCU = function(table_with_data, line_number_of_dataframe, organism) {
  new_df = table_with_data[line_number_of_dataframe, ]
  new_rscu = as.data.frame(table_with_data[line_number_of_dataframe, ])
  row.names(new_rscu) = organism
  for (i in 1:ncol(new_df)) {
    new_rscu[1, i] = ncol(new_df) * new_df[1, i] / rowSums(new_df[1, ])
  }
  return(new_rscu)
}


## Aplicar o teste por aminoacido
test_Amino = function(table_with_data, amino) {
  test_amino = chisq.test(table_with_data)#,simulate.p.value = T) # O simulate.p.value é para os casos que possuem frequencias abaixo de 5.
  p_value = test_amino$p.value
  if ( p_value < 0.05 ) {
    message(paste("########### Amino:", amino))
    message(paste("pValue:", p_value))
    print(as.data.frame(table_with_data))
    rscu_test = calculate_RSCU(table_with_data, 1, "Test")
    rscu_human = calculate_RSCU(table_with_data, 2, "Human")
    print(rbind(rscu_test, rscu_human))
  }
  else {
    message(paste("########### Amino:", amino))
    message(paste("Test fail"))
    message(paste("pValue:", p_value))
  }
  return(test_amino)
}

#Input dos dados para o teste estatístico de 2 hr vs referencia
setwd("C:/Users/gonca/Downloads/Resultados_tese/")
table1 = read_excel("RSCU_rnaseq.xlsx", sheet = "Stat_test_rnaseq_2hr_up-human", col_names = T, range = cell_cols("B:BJ") )
#row.names(table1) = c("2hr_up","human_ref")
#table1 = t(table1)
View(table1)

#Tabelas individuais (para isolar aminoacidos)

table_Lys = table1[,1:2]
table_Asn = table1[,3:4]
table_Thr = table1[,5:8]
table_Arg = table1[,9:14]
table_Ser = table1[,15:20]
table_Ile = table1[,21:23]
##table_Met = table1[,24]
table_Phe = table1[,25:26]
table_Tyr = table1[,27:28]
table_Cys = table1[,29:30]
#table_Trp = table1[,31]
table_Leu = table1[,32:37]
table_Pro = table1[,38:41]
table_His = table1[,42:43]
table_Gln = table1[,44:45]
table_Val = table1[,46:49]
table_Ala = table1[,50:53]
table_Gly = table1[,54:57]
table_Asp = table1[,58:59]
table_Glu = table1[,60:61]

p_values_list = list()

#Testes de independência do qui-quadrado para as tabelas de cada aminoácido

test_Lys = test_Amino(table_Lys, "Lys")
p_values_list = append(p_values_list, test_Lys$p.value)
test_Asn = test_Amino(table_Asn, "Asn")
p_values_list = append(p_values_list, test_Asn$p.value)
test_Thr = test_Amino(table_Thr, "Thr")
p_values_list = append(p_values_list, test_Thr$p.value)
test_Arg = test_Amino(table_Arg, "Arg")
p_values_list = append(p_values_list, test_Arg$p.value)
test_Ser = test_Amino(table_Ser, "Ser")
p_values_list = append(p_values_list, test_Ser$p.value)
test_Ile = test_Amino(table_Ile, "Ile")
p_values_list = append(p_values_list, test_Ile$p.value)
test_Phe = test_Amino(table_Phe, "Phe")
p_values_list = append(p_values_list, test_Phe$p.value)
test_Tyr = test_Amino(table_Tyr, "Tyr")
p_values_list = append(p_values_list, test_Tyr$p.value)
test_Cys = test_Amino(table_Cys, "Cys")
p_values_list = append(p_values_list, test_Cys$p.value)
test_Leu = test_Amino(table_Leu, "Leu")
p_values_list = append(p_values_list, test_Leu$p.value)
test_Pro = test_Amino(table_Pro, "Pro")
p_values_list = append(p_values_list, test_Pro$p.value)
test_His = test_Amino(table_His, "His")
p_values_list = append(p_values_list, test_His$p.value)
test_Gln = test_Amino(table_Gln, "Gln")
p_values_list = append(p_values_list, test_Gln$p.value)
test_Val = test_Amino(table_Val, "Val")
p_values_list = append(p_values_list, test_Val$p.value)
test_Ala = test_Amino(table_Ala, "Ala")
p_values_list = append(p_values_list, test_Ala$p.value)
test_Gly = test_Amino(table_Gly, "Gly")
p_values_list = append(p_values_list, test_Gly$p.value)
test_Asp = test_Amino(table_Asp, "Asp")
p_values_list = append(p_values_list, test_Asp$p.value)
test_Glu = test_Amino(table_Glu, "Glu")
p_values_list = append(p_values_list, test_Glu$p.value)


p_values_list

######## Exportar p values para txt
file_path = "C:/Users/gonca/Downloads/Resultados_tese/p_values_RSCU.txt"
file_conn = file(file_path, "w")

for (i in 1:18) {
  line = paste(p_values_list[i])
  writeLines(line, file_conn)
}
close(file_conn)

####################



ajusted_p_values_list = p.adjust(p_values_list,method = "bonferroni", n = 18 )





######## Exportar adjusted p values para txt
file_path = "C:/Users/gonca/Downloads/Resultados_tese/adjusted_p_values_RSCU.txt"
file_conn = file(file_path, "w")

for (i in 1:18) {
  line = paste(ajusted_p_values_list[i])
  writeLines(line, file_conn)
}
close(file_conn)





#################################################
#                                               #
#                                               #
#    Teste de independencia do Qui-quadrado     #
#                 Ribo-Seq data                 #
#                                               #
#################################################




### calcular RSCUs
calculate_RSCU = function(table_with_data, line_number_of_dataframe, organism) {
  new_df = table_with_data[line_number_of_dataframe, ]
  new_rscu = as.data.frame(table_with_data[line_number_of_dataframe, ])
  row.names(new_rscu) = organism
  for (i in 1:ncol(new_df)) {
    new_rscu[1, i] = ncol(new_df) * new_df[1, i] / rowSums(new_df[1, ])
  }
  return(new_rscu)
}


## Aplicar o teste por aminoacido
test_Amino = function(table_with_data, amino) {
  test_amino = chisq.test(table_with_data)#,simulate.p.value = T) # O simulate.p.value é para os casos que possuem frequencias abaixo de 5.
  p_value = test_amino$p.value
  if ( p_value < 0.05 ) {
    message(paste("########### Amino:", amino))
    message(paste("pValue:", p_value))
    print(as.data.frame(table_with_data))
    rscu_test = calculate_RSCU(table_with_data, 1, "Test")
    rscu_human = calculate_RSCU(table_with_data, 2, "Human")
    print(rbind(rscu_test, rscu_human))
  }
  else {
    message(paste("########### Amino:", amino))
    message(paste("Test fail"))
    message(paste("pValue:", p_value))
  }
  return(test_amino)
}

#Input dos dados para o teste estatístico de 2 hr vs referencia
setwd("C:/Users/gonca/Downloads/Resultados_tese/")
#setwd("/home/projects/ua/master_2/2022/master/GoncaloCosta/rscu")
table1 = read_excel("RSCU_riboseq.xlsx", sheet = "Stat_test_ribo_2hr_up-human", col_names = T, range = cell_cols("B:BJ") )
#row.names(table1) = c("2hr_up","human_ref")
#table1 = t(table1)
View(table1)

#Tabelas individuais (para isolar aminoacidos)

table_Lys = table1[,1:2]
table_Asn = table1[,3:4]
table_Thr = table1[,5:8]
table_Arg = table1[,9:14]
table_Ser = table1[,15:20]
table_Ile = table1[,21:23]
##table_Met = table1[,24]
table_Phe = table1[,25:26]
table_Tyr = table1[,27:28]
table_Cys = table1[,29:30]
#table_Trp = table1[,31]
table_Leu = table1[,32:37]
table_Pro = table1[,38:41]
table_His = table1[,42:43]
table_Gln = table1[,44:45]
table_Val = table1[,46:49]
table_Ala = table1[,50:53]
table_Gly = table1[,54:57]
table_Asp = table1[,58:59]
table_Glu = table1[,60:61]

p_values_list = list()

#Testes de independência do qui quadrado para as tabelas de cada aminoácido

test_Lys = test_Amino(table_Lys, "Lys")
p_values_list = append(p_values_list, test_Lys$p.value)
test_Asn = test_Amino(table_Asn, "Asn")
p_values_list = append(p_values_list, test_Asn$p.value)
test_Thr = test_Amino(table_Thr, "Thr")
p_values_list = append(p_values_list, test_Thr$p.value)
test_Arg = test_Amino(table_Arg, "Arg")
p_values_list = append(p_values_list, test_Arg$p.value)
test_Ser = test_Amino(table_Ser, "Ser")
p_values_list = append(p_values_list, test_Ser$p.value)
test_Ile = test_Amino(table_Ile, "Ile")
p_values_list = append(p_values_list, test_Ile$p.value)
test_Phe = test_Amino(table_Phe, "Phe")
p_values_list = append(p_values_list, test_Phe$p.value)
test_Tyr = test_Amino(table_Tyr, "Tyr")
p_values_list = append(p_values_list, test_Tyr$p.value)
test_Cys = test_Amino(table_Cys, "Cys")
p_values_list = append(p_values_list, test_Cys$p.value)
test_Leu = test_Amino(table_Leu, "Leu")
p_values_list = append(p_values_list, test_Leu$p.value)
test_Pro = test_Amino(table_Pro, "Pro")
p_values_list = append(p_values_list, test_Pro$p.value)
test_His = test_Amino(table_His, "His")
p_values_list = append(p_values_list, test_His$p.value)
test_Gln = test_Amino(table_Gln, "Gln")
p_values_list = append(p_values_list, test_Gln$p.value)
test_Val = test_Amino(table_Val, "Val")
p_values_list = append(p_values_list, test_Val$p.value)
test_Ala = test_Amino(table_Ala, "Ala")
p_values_list = append(p_values_list, test_Ala$p.value)
test_Gly = test_Amino(table_Gly, "Gly")
p_values_list = append(p_values_list, test_Gly$p.value)
test_Asp = test_Amino(table_Asp, "Asp")
p_values_list = append(p_values_list, test_Asp$p.value)
test_Glu = test_Amino(table_Glu, "Glu")
p_values_list = append(p_values_list, test_Glu$p.value)


p_values_list

######## Exportar p values para txt
file_path = "C:/Users/gonca/Downloads/Resultados_tese/p_values_RSCU.txt"
file_conn = file(file_path, "w")

for (i in 1:18) {
  line = paste(p_values_list[i])
  writeLines(line, file_conn)
}
close(file_conn)

####################



ajusted_p_values_list = p.adjust(p_values_list,method = "bonferroni", n = 18 )



######## Exportar adjusted p values para txt
file_path = "C:/Users/gonca/Downloads/Resultados_tese/adjusted_p_values_RSCU.txt"
file_conn = file(file_path, "w")

for (i in 1:18) {
  line = paste(ajusted_p_values_list[i])
  writeLines(line, file_conn)
}
close(file_conn)




#################################################
#                                               #
#                                               #
#      Testes estatisticos aos tRFs             # 
#              sncRNA-Seq data                  #
#                                               #
#################################################



setwd("C:/Users/gonca/Downloads/Resultados_tese")
table2 = read_excel("IAV_readcounts.xlsx", sheet = "only_partialpre_iav_stat_test", col_names = T)
table2 = as.data.frame(table2)

#View(table2)




table2$tRNA_fragment = factor(table2$tRNA_fragment)
long_data = gather(table2, timepoint, counts, "Mock":"8 hours", factor_key=T)
long_data$timepoint = factor(long_data$timepoint, levels = unique(long_data$timepoint), ordered = TRUE)
#View(long_data)

# p-values list
p_values_list = list()

# unique tRNA fragments
tRNA_fragments_list = unique(long_data$tRNA_fragment)


#################################################
#                                               #
#                                               #
#          Kruskal test and dunn test           #
#                                               #
#                                               #
#################################################




##############  Kruskal test  ###################

for (tRNA_fragment in tRNA_fragments_list) {
  # Aplicar o teste a 1 tRNA fragment de cada vez
  tRNA_data = long_data[long_data$tRNA_fragment == tRNA_fragment, ]
  
  kruskal_result = kruskal.test(counts ~ timepoint, data = tRNA_data)
  
  p_value = kruskal_result$p.value
  p_values_list[[tRNA_fragment]] = p_value
}


#p_values_list



#Exportar p_values_list

file_path = "C:/Users/gonca/Downloads/Resultados_tese/p_values.txt"
file_conn = file(file_path, "w")

for (i in seq_along(p_values_list)) {
  tRNA_fragment = names(p_values_list)[i]
  p_value = p_values_list[[i]]
  line = paste(tRNA_fragment, p_value, sep = "\t")
  writeLines(line, file_conn)
}
close(file_conn)







##############  Dunn test  ###################

for (tRNA_fragment in tRNA_fragments_list) {
  # Subset the data for the current tRNA fragment
  tRNA_data = long_data[long_data$tRNA_fragment == tRNA_fragment, ]
  
  # Perform ANOVA test
  dunn_result = tryCatch({
    dunnTest(tRNA_data$counts, tRNA_data$timepoint, method="bonferroni")
  }, error = function(e) {
    return(NULL)  # Return NULL if there's an error
  })
  
  # Check if dunn_test_result is NULL (all zero values) and skip the iteration
  if (is.null(dunn_result)) {
    next
  }
  
  
  # Store the p-value in the list
  
  
  p_values_list = append(p_values_list, dunn_result$res[4][1,1])
  p_values_list = append(p_values_list, dunn_result$res[4][2,1])
  p_values_list = append(p_values_list, dunn_result$res[4][3,1])
  p_values_list = append(p_values_list, dunn_result$res[4][4,1])
  p_values_list = append(p_values_list, dunn_result$res[4][5,1])
  p_values_list = append(p_values_list, dunn_result$res[4][6,1])
  
  
}


p_values_list_names = list()
#####Excluir pq todas as contagens para os tRFs são 0 em todas as replicas
excluded_tRNA_fragments_five_prime= c("Val-CAC-14","Thr-TGT-3","Thr-CGT-2","Thr-CGT-3","Met-CAT-4","Met-CAT-6")
excluded_tRNA_fragments_three_prime= c("Cys-GCA-15","Cys-GCA-19","Cys-GCA-4","Gln-CTG-5","Gln-CTG-7","Gln-TTG-2","Gly-CCC-3","Gly-GCC-4","Gly-GCC-5","Gly-TCC-4","Ile-TAT-3","Lys-CTT-11","Lys-CTT-5","Lys-CTT-7","Lys-TTT-11","Lys-TTT-5","Pro-CGG-2","Ser-GCT-2","Ser-GCT-4","Ser-GCT-6","Thr-TGT-6","Tyr-GTA-9","Val-AAC-3","Val-AAC-4","Val-CAC-14","Val-CAC-4","Val-CAC-6","iMet-CAT-2")
excluded_tRNA_fragments_partialprecounts_prime = "Ser-GCT-6-1"

for (tRNA_fragment in tRNA_fragments_list) {
  if (tRNA_fragment %in% excluded_tRNA_fragments_partialprecounts_prime) { 
    #trocar no if entre os excluded conforme os dados utilizados
    next }
  else {
    p_values_list_names = append(p_values_list_names, rep(tRNA_fragment,6))
  }
}


length(p_values_list)
length(p_values_list_names)


names(p_values_list) = p_values_list_names


#Store p_values_list

file_path = "C:/Users/gonca/Downloads/Resultados_tese/p_values.txt"
file_conn = file(file_path, "w")

for (i in seq_along(p_values_list)) {
  tRNA_fragment = names(p_values_list)[i]
  p_value = p_values_list[[i]]
  line = paste(tRNA_fragment, p_value, sep = "\t")
  writeLines(line, file_conn)
}
close(file_conn)








#################################################
#                                               #
#                                               #
# Mann-Whitney U test (Wilcoxon rank-sum test)  #
#                                               #
#                                               #
#################################################


adjusted_p_values_list = list()
mannwhitney_results_p_value_adjusted = list()
mannwhitney_results = list()
tRNA_fragment_results = list()
adjusted_p_values_list = list()

#timepoints for comparison
timepoint_names = c("Mock", "2 hours", "6 hours", "8 hours")

for (tRNA_fragment in tRNA_fragments_list) {
  # Filter the data for the current tRNA fragment
  tRNA_data = long_data[long_data$tRNA_fragment == tRNA_fragment, ]
  
  
  
  for (i in 1:(length(timepoint_names) - 1)) {
    for (j in (i + 1):length(timepoint_names)) {
      timepoint1 = timepoint_names[i]
      timepoint2 = timepoint_names[j]
      
      # Filter the data for the two timepoints
      counts_timepoint1 = tRNA_data$counts[tRNA_data$timepoint == timepoint1]
      counts_timepoint2 = tRNA_data$counts[tRNA_data$timepoint == timepoint2]
      
      
      mwu_result = wilcox.test(counts_timepoint1, counts_timepoint2, exact = F)
      
      
      p_value = mwu_result$p.value
      
      # Adjust the p-value using the Bonferroni correction
      adjusted_p_value = p.adjust(p_value, method = "bonferroni")
      
      
      test_name = paste(timepoint1, "vs", timepoint2)
      tRNA_fragment_results[[test_name]] = mwu_result
      adjusted_p_values_list[[test_name]] = adjusted_p_value
    }
  }
  
  
  mannwhitney_results[[tRNA_fragment]] = tRNA_fragment_results
  mannwhitney_results_p_value_adjusted[[tRNA_fragment]] = adjusted_p_values_list
}




comparisons_list = list()
comparisons_list = append(comparisons_list, c("Mock vs 2 hours","Mock vs 6 hours","Mock vs 8 hours","2 hours vs 6 hours","2 hours vs 8 hours","6 hours vs 8 hours"))


file_path = "C:/Users/gonca/Downloads/Resultados_tese/p_values_mann_whitney.txt"
file_conn = file(file_path, "w")

for (i in tRNA_fragments_list) {
  tRNA_fragment = i
  for (comparisons in comparisons_list) {
    p_value = mannwhitney_results[[tRNA_fragment]][[comparisons]]$p.value
    p_value_adjusted = mannwhitney_results_p_value_adjusted[[tRNA_fragment]][[comparisons]]
    line = paste(tRNA_fragment, comparisons , p_value, p_value_adjusted, sep = "\t")
    writeLines(line, file_conn)
  }
}
close(file_conn)














