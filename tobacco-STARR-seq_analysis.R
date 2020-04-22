library(readr)
library(tidyr)
library(dplyr)

# load subassembly
subassembly <- read_tsv('path to subassembly .tsv file')

# load experiment data
data.inp <- read_table('path to input barcode count file', col_names = c('count', 'barcode'))
data.out <- read_table('path to cDNA barcode count file', col_names = c('count', 'barcode'))

# read count cutoff
rcc <- 5

data.inp <- filter(data.inp, count >= rcc)
data.out <- filter(data.out, count >= rcc)

# merge data
data.both <- inner_join(data.inp, data.out, by = 'barcode', suffix = c('.inp', '.out'))

# merge with subassembly
data.both <- inner_join(data.both, subassembly, by = 'barcode')
  
# calculate enrichment
data.both <- data.both %>%
  mutate(
    enrichment = log2((count.out / sum(count.out)) / (count.inp / sum(count.inp)))
  )

# normalize
data.both <- data.both %>%
  mutate(
    enrichment = enrichment - median(enrichment[enhancer == 'none'])
  )

# aggregate by variant/fragment (for mutagenesis and pZS*11_4enh fragments)
data.both.ag <- data.both %>%
  group_by(start, stop, strand, length) %>% # adjust grouping variables to fit data
  summarise(
    count = n(),
    sd = sd(enrichment),
    enrichment = median(enrichment)
  )

# calculate "enrichment coverage" (for pZS*11_4enh fragments)
coverage <- tibble(
  pos = seq_len(5283), # length of plasmid
  enrichment = 0,
  raw = 0,
)

for (i in seq(dim(data.both.ag)[1])) {
  coverage$enrichment[data.both.ag$start[i]:data.both.ag$stop[i]] <- coverage$enrichment[data.both.ag$start[i]:data.both.ag$stop[i]] + data.both.ag$enrichment[i]
  coverage$raw[data.both.ag$start[i]:data.both.ag$stop[i]] <- coverage$raw[data.both.ag$start[i]:data.both.ag$stop[i]] + 1
}

coverage <- coverage %>%
  mutate(
    enrichment = enrichment / raw,
  )

# remove undersampled areas
coverage <- coverage %>%
  mutate(
    enrichment = if_else(raw >= 5, enrichment, NA_real_)
  )