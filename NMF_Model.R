#install.packages("ggalluvial")
library(NMF)
library(tidyverse)
library(ggalluvial)
set.seed(123456789)  # Impostiamo un seed per ridurre la casualità 

# Requisiti della matrice:
# - deve essere non negativa
# - i geni devono essere sulle righe
# - i pazienti devono essere sulle colonne


# ===============================
# 1. Caricamento dei dati
# ===============================

# Carichiamo la matrice RNA-Seq da file CSV
#matrice_espressione <- read.csv("C:\\Users\\saran\\Downloads\\rna_norm.csv", header = TRUE, stringsAsFactors = FALSE)
matrice_espressione <-tpm_matrix_top_mad
# Visualizziamo le prime righe
head(matrice_espressione)

# Stampiamo le dimensioni (geni x pazienti)
print(dim(matrice_espressione))  #  #[1] 2953  123

# Visualizzazione completa 
View(matrice_espressione)

# Verifichiamo la presenza di valori negativi
any(matrice_espressione < 0) #FALSE

# Rimuoviamo la prima colonna (presumibilmente nomi dei geni)
#matrice_espressione <- matrice_espressione[, -1]

# Controlliamo ancora valori negativi
#any(matrice_espressione < 0)

# ===============================
# 2. Parametri iniziali
# ===============================


# 3) Definiamo alcune variabili utili
numero_cluster_iniziale <- 2
numero_cluster_finale <- 6
sequenza_cluster <- seq(numero_cluster_iniziale, numero_cluster_finale, 1)

# 4) Approccio 1: validazione personalizzata

# Calcoliamo il valore cophentic per ogni K su tutte le osservazioni
indice <- 1
valori_cophenetic_totali <- vector(mode = "logical", length = length(sequenza_cluster))
matrici_consenso <- list()
risultati_nmf <- list()
matrici_h <- list()
matrici_w <- list()

#seed dell'articolo 123 #nmero di run 10
for (indice in seq_along(sequenza_cluster)) {
  cluster_corrente <- sequenza_cluster[indice]
  print(paste("Numero di cluster", cluster_corrente, "in NMF su tutte le osservazioni"))
  risultati_nmf[[indice]] <- nmf(matrice_espressione, cluster_corrente, method = "brunet", seed = 123, nrun = 100)
  valori_cophenetic_totali[[indice]] <- cophcor(risultati_nmf[[indice]])
  matrici_consenso[[indice]] <- consensus(risultati_nmf[[indice]])
  matrici_h[[indice]] <- coef(risultati_nmf[[indice]])
  matrici_w[[indice]] <- basis(risultati_nmf[[indice]])
}

# Salviamo i risultati
saveRDS(valori_cophenetic_totali, "C:\\Users\\saran\\Downloads\\ValoriCophentic_Totali_BRCA.rds")
saveRDS(matrici_consenso, "C:\\Users\\saran\\Downloads\\MatriciConsenso_Totali_BRCA.rds")
saveRDS(matrici_w, "C:\\Users\\saran\\Downloads\\MatriciW_Totali_BRCA.rds")
saveRDS(matrici_h, "C:\\Users\\saran\\Downloads\\MatriciH_Totali_BRCA.rds")
saveRDS(risultati_nmf, "C:\\Users\\saran\\Downloads\\RisultatiNMF_Totali_BRCA.rds")

# Carichiamo i risultati
valori_cophenetic_totali <- readRDS("C:\\Users\\saran\\Downloads\\ValoriCophentic_Totali_BRCA.rds")
matrici_consenso <- readRDS("C:\\Users\\saran\\Downloads\\MatriciConsenso_Totali_BRCA.rds")
matrici_w <- readRDS("C:\\Users\\saran\\Downloads\\MatriciW_Totali_BRCA.rds")
matrici_h <- readRDS("C:\\Users\\saran\\Downloads\\MatriciH_Totali_BRCA.rds")
risultati_nmf <- readRDS("C:\\Users\\saran\\Downloads\\RisultatiNMF_Totali_BRCA.rds")

#--------------------VERIFICA PAZIENTI NEI CLUSTER------------------#
#-------------------------------------------------------------------#
# Supponiamo tu abbia scelto k = 2, che corrisponde al primo elemento (indice 1)
H_matrix <- matrici_h[[1]]  # matrice H per K = 2

# Cluster assignment: ogni paziente è assegnato al cluster con valore H più alto
cluster_pazienti <- apply(H_matrix, 2, which.max)

# Assegna i nomi dei pazienti (se non già presenti)
names(cluster_pazienti) <- colnames(matrice_espressione)

# Pazienti per ogni cluster
split(names(cluster_pazienti), cluster_pazienti)

# Crea data.frame con pazienti e cluster
df_cluster <- data.frame(
  paziente = names(cluster_pazienti),
  cluster = cluster_pazienti
)
view(df_cluster)
# Salva in CSV
write.csv(df_cluster, "C:\\Users\\saran\\Downloads\\cluster_pazienti_k2.csv", row.names = FALSE)


#---------------------------VERIFICA GENI ---------------------------#
#--------------------------------------------------------------------#

W_matrix <- matrici_w[[1]]  # matrice W per K = 2 (o altro K scelto)

# Cluster assignment: ogni gene è assegnato al cluster con valore W più alto
cluster_geni <- apply(W_matrix, 1, which.max)

# Assegna i nomi dei geni (se non già presenti)
names(cluster_geni) <- rownames(W_matrix)

# Geni per ogni cluster
split(names(cluster_geni), cluster_geni)

#--------------------VERIFICA PAZIENTI NEI CLUSTER (k = 3)------------------#
#-------------------------------------------------------------------------#
H_matrix <- matrici_h[[2]]  # matrice H per K = 3

# Cluster assignment: ogni paziente è assegnato al cluster con valore H più alto
cluster_pazienti <- apply(H_matrix, 2, which.max)

# Assegna i nomi dei pazienti (se non già presenti)
names(cluster_pazienti) <- colnames(matrice_espressione)

# Pazienti per ogni cluster
pazienti_per_cluster <- split(names(cluster_pazienti), cluster_pazienti)

print(pazienti_per_cluster)

# Crea data.frame con pazienti e cluster
df_cluster <- data.frame(
  paziente = names(cluster_pazienti),
  cluster = cluster_pazienti
)

# Salva in CSV
write.csv(df_cluster, "C:\\Users\\saran\\Downloads\\cluster_pazienti_k3.csv", row.names = FALSE)

#---------------------------VERIFICA GENI (k = 3) -------------------------#
#-------------------------------------------------------------------------#
W_matrix <- matrici_w[[2]]  # matrice W per K = 3

# Cluster assignment: ogni gene è assegnato al cluster con valore W più alto
cluster_geni <- apply(W_matrix, 1, which.max)

# Assegna i nomi dei geni (se non già presenti)
names(cluster_geni) <- rownames(W_matrix)

# Geni per ogni cluster
geni_per_cluster <- split(names(cluster_geni), cluster_geni)

print(geni_per_cluster)

#--------------------Alluvion PLOT GENI K2 K3 -------------------------#
#----------------------------------------------------------------------#
# Cluster k=2 per i geni
W_k2 <- matrici_w[[1]]
cluster_geni_k2 <- apply(W_k2, 1, which.max)
names(cluster_geni_k2) <- rownames(W_k2)

# Cluster k=3 per i geni
W_k3 <- matrici_w[[2]]
cluster_geni_k3 <- apply(W_k3, 1, which.max)
names(cluster_geni_k3) <- rownames(W_k3)

# Crea dataframe con le due assegnazioni per ogni gene
df_clusters_geni <- data.frame(
  gene = names(cluster_geni_k2),
  cluster_k2 = as.factor(cluster_geni_k2),
  cluster_k3 = as.factor(cluster_geni_k3)
)

# Controlla la struttura (opzionale)
head(df_clusters_geni)

# Alluvial plot per i geni
ggplot(df_clusters_geni,
       aes(y = 1, axis1 = cluster_k2, axis2 = cluster_k3)) +
  geom_alluvium(aes(fill = cluster_k2), width = 1/12, alpha = 0.8) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, color = "white") +
  scale_x_discrete(limits = c("k=2", "k=3"), expand = c(.1, .1)) +
  labs(title = "Alluvial plot: geni da cluster k=2 a cluster k=3",
       y = "Numero geni",
       x = "Numero cluster",
       fill = "Cluster k=2") +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  )

#--------------------Alluvion PLOT CLUSTER K2-K3 -------------------------#
#--------------------------------------------------------------------#
# Ottieni cluster k=2
H_k2 <- matrici_h[[1]]
cluster_pazienti_k2 <- apply(H_k2, 2, which.max)
names(cluster_pazienti_k2) <- colnames(matrice_espressione)

# Ottieni cluster k=3
H_k3 <- matrici_h[[2]]
cluster_pazienti_k3 <- apply(H_k3, 2, which.max)
names(cluster_pazienti_k3) <- colnames(matrice_espressione)

# Ottieni cluster k=4
H_k4 <- matrici_h[[3]]
cluster_pazienti_k4 <- apply(H_k4, 2, which.max)
names(cluster_pazienti_k4) <- colnames(matrice_espressione)


# Ottieni cluster k=5
H_k5 <- matrici_h[[4]]
cluster_pazienti_k5 <- apply(H_k5, 2, which.max)
names(cluster_pazienti_k5) <- colnames(matrice_espressione)


# Crea data.frame con le due assegnazioni per ogni paziente
df_clusters <- data.frame(
  paziente = names(cluster_pazienti_k2),
  cluster_k2 = as.factor(cluster_pazienti_k2),
  cluster_k3 = as.factor(cluster_pazienti_k3)
)

# Crea data.frame con le due assegnazioni per ogni paziente
df_clusters4 <- data.frame(
  paziente = names(cluster_pazienti_k2),
  cluster_k2 = as.factor(cluster_pazienti_k2),
  cluster_k4 = as.factor(cluster_pazienti_k4)
)

# Crea data.frame con le due assegnazioni per ogni paziente
df_clusters5 <- data.frame(
  paziente = names(cluster_pazienti_k2),
  cluster_k2 = as.factor(cluster_pazienti_k2),
  cluster_k5 = as.factor(cluster_pazienti_k5)
)

# Controlla la struttura (opzionale)
head(df_clusters)

# Crea alluvial plot
ggplot(df_clusters,
       aes(y = 1, axis1 = cluster_k2, axis2 = cluster_k3)) +
  geom_alluvium(aes(fill = cluster_k2), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("k=2", "k=3"), expand = c(.1, .1)) +
  labs(title = "Alluvial plot: pazienti da cluster k=2 a cluster k=3",
       y = "Numero pazienti",
       x = "Numero cluster") +
  theme_minimal()


# Crea alluvial plot
ggplot(df_clusters4,
       aes(y = 1, axis1 = cluster_k2, axis2 = cluster_k4)) +
  geom_alluvium(aes(fill = cluster_k2), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("k=2", "k=4"), expand = c(.1, .1)) +
  labs(title = "Alluvial plot: pazienti da cluster k=2 a cluster k=4",
       y = "Numero pazienti",
       x = "Numero cluster") +
  theme_minimal()


# Crea alluvial plot
ggplot(df_clusters5,
       aes(y = 1, axis1 = cluster_k2, axis2 = cluster_k5)) +
  geom_alluvium(aes(fill = cluster_k2), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("k=2", "k=5"), expand = c(.1, .1)) +
  labs(title = "Alluvial plot: pazienti da cluster k=2 a cluster k=5",
       y = "Numero pazienti",
       x = "Numero cluster") +
  theme_minimal()
#-----------------------GRAFICO COEFFICENTE--------------------------#
#--------------------------------------------------------------------#
# Grafico valori cophentic totali
valori_cophenetic_totali <- round(valori_cophenetic_totali, digits = 3)
dati_cophenetic_totali <- data.frame(NumeroCluster = seq(2, 6, 1), Valori = valori_cophenetic_totali)

grafico_totale <- ggplot(data = dati_cophenetic_totali, aes(x = NumeroCluster, y = Valori)) +
  geom_line(color = "#D0E6A5", size = 2, alpha = 0.5) +
  geom_point(color = "#D0E6A5", size = 1.5) +
  geom_text(aes(label = Valori)) +
  scale_x_continuous(breaks = seq(2, 6, 1)) +
  labs(title = "Indice Cophenetico – Tutti i campioni")

grafico_totale