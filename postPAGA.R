library(zellkonverter)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(tibble)
library(wesanderson)

sub_name = 'Monocytes_remove_dead'
file_path = file.path('/Users/tan/Kancera', 'PAGA_result_data', 
                      paste0('all_', sub_name, '_sample1000.h5ad'))

datah5 = readH5AD(file = file_path)
data = t(assay(datah5)) %>% as_tibble()
meta = colData(datah5) %>% 
  as_tibble() %>% 
  select(timepoint, Sample.ID, Subject.ID, group, leiden_new) %>%
  rename(randtrt = group) %>%
  rename(leiden = leiden_new)

# subpop freq change
data_freq = meta %>% 
  group_by(timepoint, randtrt, leiden) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(timepoint, randtrt) %>% # for each timepoint x randtrt, frequency add up to 1
  mutate(freq = n / sum(n))

g_list_0 <- lapply(sort(unique(meta$leiden)), function(x){
  ggplot(data = data_freq %>% filter(leiden == x), 
         aes(x = timepoint, y = freq, color = randtrt, group = randtrt)) +
    geom_line() +
    geom_point() +
    labs(title = paste0(sub_name, '_cluster_', x))
})
ggarrange(plotlist = g_list_0, ncol = 4, nrow = 3) %>% 
  ggexport(filename = file.path(getwd(), '../figures', paste0(sub_name, '_cluster_changes_timepoint.pdf')),
           width = 20,
           height = 10)
################################

# subpop marker expression
data_anno <- data %>% add_column(leiden = meta$leiden,
                                 randtrt = meta$randtrt)

g_list_1 <- lapply(colnames(data), function(x){
  d <- data_anno %>% select(x, leiden) %>% rename(marker_expression = x)
  ggplot(d, 
         aes(x = marker_expression, y = leiden, fill = leiden)) +
    geom_density_ridges(scale = 4, rel_min_height=.01) +
    scale_fill_manual(values = wes_palette("Rushmore1", nlevels(meta$leiden), type = "continuous")) +
    theme_ridges() + 
    theme(legend.position = "none") +
    xlim(quantile(d$marker_expression, 0.01),
         quantile(d$marker_expression, 0.99)) +
    labs(title = x)
})
ggarrange(plotlist = g_list_1, ncol = 5, nrow = 9) %>% 
  ggexport(filename = file.path(getwd(), '../figures', paste0(sub_name, '_cluster_marker_expressions.pdf')),
           width = 20,
           height = 40)
################################
# subpop marker expression separate placebo/treatment

g_list_1_1 <- lapply(colnames(data), function(x){
  d <- data_anno %>% 
    filter(randtrt == 'PLACEBO') %>%
    select(all_of(x), leiden) %>% rename(marker_expression = x)
  ggplot(d, 
         aes(x = marker_expression, y = leiden, fill = leiden)) +
    geom_density_ridges(scale = 4, rel_min_height=.01) +
    scale_fill_manual(values = wes_palette("Rushmore1", nlevels(meta$leiden), type = "continuous")) +
    theme_ridges() + 
    theme(legend.position = "none") +
    xlim(quantile(d$marker_expression, 0.01),
         quantile(d$marker_expression, 0.99)) +
    labs(title = x)
})
ggarrange(plotlist = g_list_1_1, ncol = 5, nrow = 9) %>% 
  ggexport(filename = file.path(getwd(), '../figures', paste0(sub_name, '_cluster_marker_expressions_PLACEBO.pdf')),
           width = 20,
           height = 40)
g_list_1_2 <- lapply(colnames(data), function(x){
  d <- data_anno %>% 
    filter(randtrt == 'KAND567') %>%
    select(all_of(x), leiden) %>% rename(marker_expression = x)
  ggplot(d, 
         aes(x = marker_expression, y = leiden, fill = leiden)) +
    geom_density_ridges(scale = 4, rel_min_height=.01) +
    scale_fill_manual(values = wes_palette("Rushmore1", nlevels(meta$leiden), type = "continuous")) +
    theme_ridges() + 
    theme(legend.position = "none") +
    xlim(quantile(d$marker_expression, 0.01),
         quantile(d$marker_expression, 0.99)) +
    labs(title = x)
})
ggarrange(plotlist = g_list_1_2, ncol = 5, nrow = 9) %>% 
  ggexport(filename = file.path(getwd(), '../figures', paste0(sub_name, '_cluster_marker_expressions_KAND567.pdf')),
           width = 20,
           height = 40)


################################
# subpop freq change by subject

data_freq = meta %>% 
  group_by(timepoint, Subject.ID, leiden) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(Subject.ID) %>% # for each timepoint x subject, frequency add up to 1
  mutate(freq = n / sum(n)) %>%
  left_join(meta %>% select(Subject.ID, randtrt) %>% distinct(), by='Subject.ID')

g_list_2 <- lapply(sort(unique(meta$leiden)), function(x){
  ggplot(data = data_freq %>% filter(leiden == x), 
         aes(x = timepoint, y = freq, color = randtrt)) +
    geom_violin(scale="width") +
    geom_jitter(aes(group=randtrt), position=position_jitterdodge()) + 
    stat_summary(fun.y="mean",geom="crossbar", mapping=aes(ymin=..y.., ymax=..y..), 
                 width=1, position=position_dodge(),show.legend = FALSE) + 
    labs(title = paste0(sub_name, '_cluster_', x))
})
ggarrange(plotlist = g_list_2, ncol = 4, nrow = 3) %>% 
  ggexport(filename = file.path(getwd(), '../figures', paste0(sub_name, '_cluster_changes_timepoint_by_subject_2.pdf')),
           width = 20,
           height = 10)


################################

# subject distribution of each cluster
# g_list_3 <- lapply(sort(unique(meta$leiden)), function(cluster_name){
#   df <- meta %>% 
#     filter(leiden == cluster_name) %>%
#     group_by(Subject.ID) %>%
#     summarise(count=n()) %>%
#     left_join(meta %>% select(Subject.ID, randtrt) %>% distinct(), by='Subject.ID') %>%
#     arrange(randtrt) %>%
#     mutate(
#       fraction = count/sum(count),
#       ymax = cumsum(fraction),
#       ymin = c(0, head(ymax, n=-1)),
#     )
#   breaks <- df %>% group_by(randtrt) %>% summarise(label_frac=sum(fraction)) %>% 
#     mutate(max = cumsum(label_frac), 
#            min = c(0,head(max,n=-1)),
#            pos = (max+min)/2
#     )
#   # make a named vector for color
#   cols = wes_palette("Rushmore1")[c(1,4)]
#   names(cols) = levels(breaks$randtrt)
#   g <- ggplot(df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=2.7, fill=randtrt)) +
#     geom_rect(color = 'white', size=0.5) +
#     annotate("text", x=2, y=0, label = sum(df$count), size = 8) +
#     coord_polar(theta="y") +
#     xlim(c(2, 4)) +
#     scale_y_continuous(breaks = breaks$pos, labels = breaks$randtrt) +
#     scale_fill_manual(values = cols) +
#     theme(axis.ticks = element_blank(),
#           axis.title = element_blank(),
#           axis.text.y = element_blank(),
#           axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
#           axis.text = element_text(size = 12), 
#           legend.position = "none",
#           panel.background = element_rect(fill = "white"),
#           plot.title = element_text(hjust = 0.5)) +
#     labs(title = paste0(sub_name, '_cluster ', cluster_name))
#   return(g)
# })
# ggarrange(plotlist = g_list_3, ncol = 3, nrow = 4) %>% 
#   ggexport(filename = file.path(getwd(), '../figures', paste0(sub_name, '_cluster_subject_distribution.pdf')),
#            width = 20,
#            height = 30)
