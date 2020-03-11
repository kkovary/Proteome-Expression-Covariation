library(mvtnorm)

fun = function(cor_matrix, list_distributions, L)
{
  n = length(list_distributions)
  # Correlated Gaussian variables
  Gauss = rmvnorm(n=L, mean = rep(0,n), sig=cor_matrix)
  # convert them to uniform distribution.
  Unif = pnorm(Gauss) 
  # Convert them to whatever I want
  vars = sapply(1:n, FUN = function(i) list_distributions[[i]](Unif[,i]))
  return(vars)
}

# Create a matrix of values with a given R, mean, and standard deviation
corProts <- function(nProt, nCell, R, mean, sd) {
  
  cor_matrix =  matrix(data = R,
                       ncol = nProt,
                       nrow = nProt)
  
  diag(cor_matrix) <- 1
  list_distributions <- list()
  for (i in 1:nProt) {
    list_distributions <-
      c(list_distributions, function(nCell)
        round(qnorm(nCell, mean, sd)))
  }
  
  mat = fun(cor_matrix, list_distributions, nCell)
  return(mat)
}

# Rug Plots
cor_rug_plot <- function(x, color = '#e08214'){
  temp <- filter(allGroups, accession %in% x) %>%
    dplyr::select(name, cors) %>% unnest(cols = cors)
  
  temp <- allGroups %>% select(name, cors) %>%
    unnest(cols = cors) %>% mutate(name = "all") %>%
    filter(!duplicated(cors)) %>%
    rbind(temp,.)
  
  p1 <- ggplot(temp, aes(x = cors, fill = name)) + 
    geom_density(alpha = 0.75, color = "transparent") +
    scale_fill_manual(values = c('#bababa',color), name = "") +
    geom_rug(data = filter(temp, name != "all"), color = color) +
    xlab("") + ylab("Density") +
    theme_bw() +
    # theme(legend.position = c(.2,.9),
    #       legend.background = element_rect(fill="transparent")) +
    theme(legend.position = "none",
          text = element_text(size = 6)) +
    xlim(-1,1) +
    ggtitle(filter(temp, name != "all") %>% pull(name) %>% unique())
  
  p2 <- ggplot(temp, aes(x = name, y = cors, fill = name)) + 
    geom_boxplot(notch = TRUE, alpha = 0.75) +
    scale_fill_manual(values = c('#bababa',color)) +
    ylab("Correlation Coefficients") +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = 6),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    ylim(-1,1) +
    coord_flip()
  
  # p1/p2 + 
  #   #plot_annotation(tag_levels = "A") +
  #   plot_layout(heights = c(2,1))
  plot_grid(p1, p2, ncol = 1, rel_heights = c(2,1), align = "v")
}
