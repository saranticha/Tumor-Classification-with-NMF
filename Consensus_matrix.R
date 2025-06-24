# Carichiamo il pacchetto necessario per l'NMF (Non-negative Matrix Factorization)
library(NMF)

# Carichiamo i risultati dell'NMF salvati in precedenza da un file RDS
risultatiNMF <- readRDS("C:\\Users\\saran\\Downloads\\RisultatiNMF_Totali_BRCA.rds")
#risultatiNMF <- risultati_nmf
# Estraiamo i risultati per K = 2, 3 e 4 cluster
K2 <- risultatiNMF[[1]]  # Risultato per 2 cluster
K3 <- risultatiNMF[[2]]  # Risultato per 3 cluster
K4 <- risultatiNMF[[3]]  # Risultato per 4 cluster
K5 <- risultatiNMF[[4]]  # Risultato per 5 cluster
K6 <- risultatiNMF[[5]]  # Risultato per 6 cluster
#K7 <- risultatiNMF[[6]]  # Risultato per 7 cluster
#K8 <- risultatiNMF[[7]]  # Risultato per 8 cluster
#K9 <- risultatiNMF[[8]]  # Risultato per 9 cluster
#K10 <- risultatiNMF[[9]]  # Risultato per 10 cluster

# Per ciascuna matrice di consenso, assegniamo i nomi dei campioni come nomi di righe e colonne
colnames(K2@consensus) <- sampleNames(K2)
rownames(K2@consensus) <- sampleNames(K2)

# Salviamo la mappa di consenso per K=2 in formato PDF
pdf("C:\\Users\\saran\\Downloads\\K2.pdf")
consensusmap(K2, hclustfun="average")  # Mappa di consenso con clustering gerarchico (funzione average)
dev.off()

#  K=3
colnames(K3@consensus) <- sampleNames(K3)
rownames(K3@consensus) <- sampleNames(K3)
pdf("C:\\Users\\saran\\Downloads\\K3.pdf")
consensusmap(K3, hclustfun="average")
dev.off()

#  K=4
colnames(K4@consensus) <- sampleNames(K4)
rownames(K4@consensus) <- sampleNames(K4)
pdf("C:\\Users\\saran\\Downloads\\K4.pdf")
consensusmap(K4, hclustfun="average")
dev.off()

#  K=5
colnames(K5@consensus) <- sampleNames(K5)
rownames(K5@consensus) <- sampleNames(K5)
pdf("C:\\Users\\saran\\Downloads\\K5.pdf")
consensusmap(K5, hclustfun="average")
dev.off()

#  K=6
colnames(K6@consensus) <- sampleNames(K6)
rownames(K6@consensus) <- sampleNames(K6)
pdf("C:\\Users\\saran\\Downloads\\K6.pdf")
consensusmap(K6, hclustfun="average")
dev.off()

#  K=7
colnames(K7@consensus) <- sampleNames(K7)
rownames(K7@consensus) <- sampleNames(K7)
pdf("C:\\Users\\saran\\Downloads\\K7.pdf")
consensusmap(K7, hclustfun="average")
dev.off()

#  K=8
colnames(K8@consensus) <- sampleNames(K8)
rownames(K8@consensus) <- sampleNames(K8)
pdf("C:\\Users\\saran\\Downloads\\K8.pdf")
consensusmap(K8, hclustfun="average")
dev.off()

#  K=9
colnames(K9@consensus) <- sampleNames(K9)
rownames(K9@consensus) <- sampleNames(K9)
pdf("C:\\Users\\saran\\Downloads\\9.pdf")
consensusmap(K9, hclustfun="average")
dev.off()

#  K=10
colnames(K10@consensus) <- sampleNames(K10)
rownames(K10@consensus) <- sampleNames(K10)
pdf("C:\\Users\\saran\\Downloads\\K10.pdf")
consensusmap(K10, hclustfun="average")
dev.off()