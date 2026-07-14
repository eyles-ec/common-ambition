library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(viridis)

#paths 
source("../paths.R")

#set working directory
setwd(wd)

#read in ukhsa surveillance data
#csv made with data pulled from: (https://www.gov.uk/government/statistics/hiv-annual-data)
ukhsa <- read.csv("./ukhsa_data.csv")

#read in ACHC london/sw data
ukhsa_achc <- read.csv("./ukhsa_region_ethn.csv")

#pivot surveillance to long format
hiv_long <- ukhsa %>%
  pivot_longer(
    cols = starts_with("y"),
    names_to = "year",
    values_to = "rate"
  ) %>%
  mutate(
    year = c(
      y15 = 2015, y16 = 2016, y17 = 2017, y18 = 2018, y19 = 2019,
      y20 = 2020, y21 = 2021, y22 = 2022, y23 = 2023, y24 = 2024
    )[year],
    rate = parse_number(rate),
    rate_per_1000 = rate / 100
  )

#pivot achc to long
achc_long <- ukhsa_achc %>%
  pivot_longer(
    cols = starts_with("y"),
    names_to = "year",
    values_to = "rate"
  ) %>%
  mutate(
    year = c(
      y15 = 2015, y16 = 2016, y17 = 2017, y18 = 2018, y19 = 2019,
      y20 = 2020, y21 = 2021, y22 = 2022, y23 = 2023, y24 = 2024
    )[year],
    rate_per_1000 = rate / 100
  )

#generate an 'england, not england' variable for plotting
hiv_long <- hiv_long %>%
  mutate(
    england = ifelse(area == "England", "England", "Region")
  )

#set colour palette with england specifically in red
cols <- viridis(length(unique(hiv_long$area)))
names(cols) <- sort(unique(hiv_long$area))
cols["England"] <- "#D55E00"

#plot the data with the above 
p <- ggplot(
  hiv_long,
  aes(year,
    rate_per_1000,
    colour = area,
    linewidth = england,
    linetype = england)
) +
  geom_line() +
  scale_colour_manual(values = cols) +
  scale_linewidth_manual(
    values = c(
      "Region" = 0.8,
      "England" = 1.2),
    guide = "none"
) +
  scale_linetype_manual(
    values = c(
      "Region" = "solid",
      "England" = "longdash"
    ),
    guide = "none"
) +
  scale_x_continuous(breaks = 2015:2024) +
  labs(
    x = "Year",
    y = "HIV testing rate per 1,000 population",
    colour = "Region"
  ) +
  theme_minimal()

#save plot
ggsave(
  "./Analysis/ukhsa_hiv_testing_rates.png",
  plot = p,
  width = 10,
  height = 6,
  dpi = 300
)


p2<- ggplot(
  achc_long,
  aes(
    x = year,
    y = rate_per_1000,
    colour = area,
    linetype = ethn
  )
) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = cols) +
  scale_x_continuous(breaks = 2015:2024) +
  labs(
    x = "Year",
    y = "HIV testing rate per 1,000 population",
    colour = "Area",
    linetype = "Ethnicity"
  ) +
  theme_minimal()

ggsave(
  "./Analysis/ukhsa_achc_testing_rates.png",
  plot = p2,
  width = 10,
  height = 6,
  dpi = 300
)

