# read in all dir name with all qc files
set.seed(2024)
args <- commandArgs(trailingOnly = TRUE)
rlen_plot <- function(fn_read_length_txt){
    # input file
    rl_d <- read.csv(fn_read_length_txt, header = F)
    rl_d <- rl_d %>% dplyr::rename(rl = V1)

    # plot the distribution of the read length
    # Create a data frame from the vector
    
    ggplot(rl_d %>% sample_frac(0.01, replace = FALSE), aes(x = rl)) +
    geom_histogram(fill = "#367DB0",color = NA, alpha = 1, bins=100) +
    labs(title = glue("Read length distribution of the {sample_name} data"), x = "Read length", y = "Density") + 
    scale_x_continuous(limits = c(0, 5000))

    set.seed(seed)
    rl_d <- rl_d %>% sample_frac(0.01, replace = FALSE)
    p1 <- ggplot(rl_d, aes(x = rl)) +
    geom_histogram(aes(y = ..density..), fill = "#367DB0",color = NA, alpha = 1, bins=100) +
    geom_density(fill = NA,color = "red", alpha = 1, bins=100) + 
    labs(title = glue("Read length distribution of the {sample_name} data"), y = "Density") + 
    scale_x_continuous(limits = c(0, 5000)) + theme_minimal() + theme(axis.text.y = element_blank(), 
                                                                        axis.text.x = element_blank(),
                                                                        axis.title.x= element_blank())
    p2 <- 
    ggplot(rl_d, aes(y = rl)) +
    geom_boxplot(width = 0.2, fill = "#367DB0", color = "grey") +  # Adjust width as needed
    theme_minimal()  + scale_y_continuous(limits = c(0, 5000)) + coord_flip() + 
    labs(y = "Read length", x = "Density") + 
    theme(axis.text.y = element_blank() )

    p_rl_dist <- arrangeGrob(
    p1,
    p2,
    ncol = 1,
    heights= c(2,1)
    )

    rm(rl_d, p1,p2)
    p_rl_dist %>% grid.arrange
}
