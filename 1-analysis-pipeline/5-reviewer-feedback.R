library(dplyr)
library(ggplot2)

load("1-analysis-pipeline/1-output/sd_breast_set2.Rdata")
load("1-analysis-pipeline/1-output/sd_breast_all.Rdata")

head(sd_breast)
head(sd_breast_set2)

sd_breast <- data.frame(cg = names(sd_breast),
                        sd_breast = sd_breast,
                        array = "epic")
sd_breast_set2 <- data.frame(cg = names(sd_breast_set2),
                             sd_breast = sd_breast_set2,
                             array = "k450")
sd <- rbind(sd_breast, sd_breast_set2)

sdwide <- sd |> 
  tidyr::pivot_wider(id_cols = "cg",names_from = "array",values_from = sd_breast) |> 
  dplyr::filter(!is.na(epic) & !is.na(k450))
head(sdwide)

library(ggtext)
plot <- sdwide |> 
  ggplot(aes(x = epic,
             y = k450)) +
  geom_hex(bins = 250,
           aes(colour = log(after_stat(count)),
               fill = log(after_stat(count)))) +
  geom_abline(intercept = 0, slope = 1,
              colour = "gray30",
              linetype = "dashed") +
  ggpubr::stat_cor() +
  viridis::scale_colour_viridis(aesthetics = c("fill", "colour"),
                                name = "log<sub>count</sub>",
                                option = "plasma") +
  theme_bw() +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),legend.title = element_markdown(),
        legend.position = "top",
        aspect.ratio = 1) +
  labs(x = "<b>Variability set [EPIC]</b><br>n=50",
       y = "<b>Variability set [450K]</b><br>n=146")

ggsave(plot, file = "~/Dropbox/index-dev/3C/1-BC/23-mqtl/8-manuscript/3-manuscript/5-revision-jul23/plot_corr.pdf", width = 5, height = 5)


# Top 20 % overlaps
top_epic_overlap2 <- sdwide |> 
  dplyr::arrange(desc(epic)) |> 
  slice(1:(nrow(sdwide)*0.1))
top_450_overlap2 <- sdwide |> 
  dplyr::arrange(desc(k450)) |> 
  slice(1:(nrow(sdwide)*0.1))

list <- list(`Top 10% CpGs \n(Variability\nset EPIC)` = top_epic_overlap2$cg,
             `Top 10% CpGs\n(Variability\nset 450k)` = top_450_overlap2$cg)

library(ggvenn)
plot <- ggvenn::ggvenn(list,
               fill_color = viridis::viridis(2),
               stroke_size = 0.5)
ggsave(plot, file = "~/Dropbox/index-dev/3C/1-BC/23-mqtl/8-manuscript/3-manuscript/5-revision-jul23/plot_venn.pdf", width = 5, height = 5)
