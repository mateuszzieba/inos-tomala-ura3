###########################
# require needed packages #
###########################
install.packages("magrittr")
install.packages("dplyr")
install.packages("purrr")
install.packages("tibble")
install.packages("tidyr")
install.packages("ggplot2")
install.packages("stringr")


require(magrittr)
require(dplyr)
require(purrr)
require(tibble)
require(tidyr)
require(ggplot2)
require(stringr)

experimental_sequences <- read.csv(file = "data/result_experimental_sequences.tsv", sep = "\t")%>%
  mutate(mutation = str_split_fixed(file, "[.]", 2)[,1],
         group = "experimental") 

missense_sequences <- read.csv(file = "data/result_missense_possible_sequences.tsv", sep = "\t") %>% 
  mutate(mutation = "missense",
         group = "missense")



rbind(experimental_sequences, missense_sequences) %>%
  mutate(diff_possition = MINCEN - MIN,
         pair = paste(position1, position2, sep = "_")) %>% 
  ggplot(aes(x = mutation, y = diff_possition)) +
  geom_boxplot() 


rbind({missense_sequences %>% mutate(group = "missense")},
      {experimental_sequences %>% mutate(group = "experimental")}) %>%
  mutate(diff_possition = MINCEN - MIN,
         pair = paste(position1, position2, sep = "_")) %>%
  group_by(pair) %>%
  nest() %>%
  mutate(t_test = map(data, ~t.test(.x[.$group == "missense",]$diff_possition,
                                    .x[.$group == "experimental",]$diff_possition,
                                    var.equal = TRUE)[[3]])) %>%
  unnest(t_test)


rbind({missense_sequences %>% mutate(group = "missense")},
      {experimental_sequences %>% mutate(group = "experimental")}) %>%
  mutate(diff_possition = MINCEN - MIN) %>%
  ggplot(aes(x = diff_possition, color = group)) +
  geom_histogram(position = "dodge")
