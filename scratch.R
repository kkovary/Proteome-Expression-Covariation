

hist(rnorm(10000, 100, 10) / 100)
hist(log2(rnorm(10000)))


a <- c("linear","log")

x <- rnorm(1000, mean = 1, sd = 0.1)
sd(x) / mean(x)
sd(x*100) / mean(x*100)

x = rnorm(100,100,10)
y = x * rnorm(100,100,10)



#mat <- cor(matrix(rnorm(50), nrow = 5, ncol = 10))

prot_names <- filter(allGO, GOterms == 'GO:0005840')$Uniprot %>% unlist()
mat <- cor[rownames(cor) %in% prot_names,colnames(cor) %in% prot_names]

colnames(mat) <- NULL
rownames(mat) <- NULL
hc = hclust(as.dist(1 - mat))
mat = mat[hc$order, hc$order]
mat[lower.tri(mat)] = NA
diag(mat) <- NA

pheatmap(mat, 
         cluster_col = F, 
         cluster_row = F, 
         color = rev(colorRampPalette(brewer.pal(11, 'RdBu'))(201)), 
         border_color = NA,
         na_col = "white",
         legend = FALSE,
         breaks = seq(-1,1, by = 0.01))

protImputeTidy %>% filter(Uniprot %in% prot_names) %>%
  ggplot(., aes(x = abundance, group = Uniprot, alpha = 0.05)) + geom_density(color = "transparent", fill = "grey")

protVar %>% filter(Uniprot %in% prot_names) %>%
  ggplot(., aes(x = 1, y = Uniprot, fill = sd)) + geom_tile() +
  coord_equal()


pca = prcomp(as.matrix(glyc[,2:ncol(glyc)]))

pcaDF = tibble(PC1 = pca$rotation[,1], PC2 = pca$rotation[,2])
ggplot(pcaDF, aes(x = PC1, y = PC2)) + geom_point()

int <- matrix(data = c(0,1,0,
                       0,0,1,
                       1,0,0),
              nrow = 3)
mat <- cor[1:3,1:3]
