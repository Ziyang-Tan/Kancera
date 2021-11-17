library(mmR)
tab1 = mm.fastread('/Users/tan/Kancera/data/EXP-21-DG3638_Sample info.xlsx') %>%
  select(`Subject ID`, `Sample info`, `Sample ID`)
tab2 = mm.fastread('/Users/tan/Kancera/data/EXP-21-DG3653_Sample info.xlsx') %>%
  select(`Subject ID`, `Sample info`, `Sample ID`)
tab3 = mm.fastread('/Users/tan/Kancera/data/EXP-21-DG3654_Sample info.xlsx') %>%
  select(`Subject ID`, `Sample info`, `Sample ID`)
randinfo = mm.fastread('/Users/tan/Kancera/data/Randomiseringskod.xlsx') %>%
  mutate(RANDNO = as.character(RANDNO)) %>%
  rename(`Subject ID` = RANDNO)
  


sample_info = do.call(rbind, list(tab1, tab2, tab3)) %>%
  left_join(randinfo, by='Subject ID') %>%
  filter(`Subject ID` != 'Empty')


sample_info %>% filter(`Sample ID` %in% c(8599,8592,8594,8593,100499,100505))
