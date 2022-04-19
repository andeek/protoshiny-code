##################################################
## Make hc object for mobility/covid-19 US counties
## data is from nytimes and google
## clustering will be done by mobility
## label will be cases of covid-19 since march
##################################################
start_label <- "aug1"
start_date <- "2020-08-01"
end_date <- "2021-01-15"


## libraries ----
library(tidyverse)  # manipulate data
library(tidycensus) # fips/county/state data map
library(lubridate)  # work with dates
library(protoclust) # devtools::install_github("jacobbien/protoclust")
library(zoo)
library(noncensus)  # devtools::install_github("ramhiser/noncensus")
library(covidcast)  # devtools::install_github("cmu-delphi/covidcast", ref = "main", subdir = "R-packages/covidcast")
# About the package: https://cmu-delphi.github.io/covidcast/covidcastR/

## reproducibility
set.seed(400)

## load mobility data signals
dat <- covidcast_signals(data_source = "safegraph",
                         signal = c("restaurants_visit_prop",
                                    "completely_home_prop"),
                         start_day = start_date,
                         end_day = end_date,
                         geo_type = "county")

dat[[1]] %>%
  rename(restaurant = value) %>%
  select(geo_value, time_value, restaurant) %>%
  inner_join(dat[[2]] %>% rename(home = value) %>% select(geo_value, time_value, home)) -> mobility

## filter out counties with less than 149 days
mobility %>%
  group_by(geo_value) %>%
  summarise(count = n()) %>%
  filter(count == max(count)) %>%
  pull(geo_value) -> counties

## get correlation matrix
num_counties <- length(counties)
corr <- matrix(rep(NA, num_counties^2), nrow = num_counties)

for(i in seq_len(num_counties)) {
  for(j in seq_len(num_counties)) {
    if(i < j) {
      county_a <- mobility %>% 
        filter(geo_value == counties[i]) %>%
        select(time_value, restaurant, home) %>%
        mutate(restaurant = scale(restaurant), home = scale(home))
      
      county_b <- mobility %>% 
        filter(geo_value == counties[j]) %>%
        select(time_value, restaurant, home) %>%
        mutate(restaurant = scale(restaurant), home = scale(home))
      
      ## make sure dates match
      ## don't fill in extra days anymore
      county_a %>% inner_join(county_b, by = "time_value") -> county_ab
      
      county_ab <- county_ab[complete.cases(county_ab),]
      
      county_ab %>%
        select(contains(".x")) %>%
        as.matrix() %>%
        as.numeric() -> county_a
      
      county_ab %>%
        select(contains(".y")) %>%
        as.matrix() %>%
        as.numeric -> county_b
      
      corr[i, j] <- 1 - abs(cor(county_a, county_b))
      cat(paste0("i: ", i, " j: ", j, "\r"))
    }
  }
}
corr[lower.tri(corr)] <- t(corr)[lower.tri(corr)]
diag(corr) <- rep(0, length(diag(corr)))
rownames(corr) <- counties

# many rows are all NA, drop them
if(length(which(colSums(is.na(corr)) > 1000)) > 0) {
  corr_drop_na <- corr[-which(colSums(is.na(corr)) > 1000), -which(colSums(is.na(corr)) > 1000)]
} else {
  corr_drop_na <- corr
}
corr_drop_na <- corr_drop_na[complete.cases(corr_drop_na), complete.cases(corr_drop_na)]

## get the protoclust object
hc <- protoclust(corr_drop_na)

## prepare for plotting labels
mobility_0 <- mobility %>% 
  filter(geo_value %in% rownames(corr_drop_na)) %>%
  ungroup()

mobility_filled_0 <- mobility_0  %>%
  pivot_longer(any_of(c("restaurant", "home")),
               names_to = "var",
               values_to = "prop") %>%
  pivot_wider(values_from = prop, names_from = c(time_value, var), names_sep = "__", values_fill = NA) 

mobility_filled_plot_data <- mobility_filled_0 %>% 
  pivot_longer(!any_of(c("geo_value")), names_to = "date__var", values_to = "prop") %>% 
  separate(date__var, into = c("date", "var"), sep = "__") %>%
  mutate(date = ymd(date)) %>%
  pivot_wider(names_from = var, values_from = prop)

## add regions for colors
data("states")
data("fips_codes")

mobility_filled_plot_data %>%
  left_join(fips_codes %>% mutate(geo_value = paste0(state_code, county_code)) %>% select(geo_value, state, county)) %>%
  left_join(states %>% select(state, region, division), by = c("state" = "state")) %>%
  mutate(county = gsub(" County", "", county)) -> mobility_filled_plot_data

region_colors <- data.frame(regions = unique(states$region), color = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"))
region_colors <- region_colors[!is.na(region_colors$regions),]

## plot mobility data
for(fips in unique(mobility_filled_plot_data$geo_value)) {
  dat <- mobility_filled_plot_data[mobility_filled_plot_data$geo_value == fips,] 
  
  ggplot(dat) +
    geom_point(aes(restaurant, home, colour = date)) +
    theme_void() +
    ggtitle(paste(unique(dat$county), unique(dat$state), sep = ", ")) +
    theme(legend.position = "none",
          plot.title = element_text(size = 8)) +
    scale_color_gradient(low = "grey80", high = region_colors[region_colors$region == unique(dat$region), "color"]) -> p

  
  ## save
  ifelse(!dir.exists(file.path(paste0("examples/covid/daily_plots_mobility_restaurant_home_", start_label, "/"))), dir.create(file.path(paste0("examples/covid/daily_plots_mobility_restaurant_home_", start_label, "/"))), FALSE)
  
  ggsave(p, device = "png", filename = paste0("examples/covid/daily_plots_mobility_restaurant_home_", start_label, "/", fips, ".png"), width = 1, height = 1, units = "in", dpi = 200, bg = "white")
}

## create covid case plot labels
## get county covid plots and save with fips
## data 
## Start by loading data from https://github.com/nytimes/covid-19-data 
covid <- read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")

## diff = # of new cases
covid <- covid %>%
  group_by(county, state, fips) %>%
  mutate(daily_cases = c(0, diff(cases))) %>%
  filter(!is.na(fips))

covid_filled_0 <- covid %>% 
  filter(date <= end_date & date >= start_date) %>%
  ungroup() %>%
  select(-county, -state) %>%
  right_join(mobility_filled_plot_data %>% select(geo_value, county, state) %>% unique(), by = c("fips" = "geo_value")) %>%
  select(-deaths, -cases) %>%
  replace_na(list(date = min(covid$date, na.rm = TRUE), daily_cases = 0)) %>%
  spread(date, daily_cases, fill = 0) 

covid_filled_plot_data <- covid_filled_0 %>% 
  gather(date, daily_cases, -fips, -county, -state) %>% 
  mutate(date = ymd(date))

## make small plots of cases for each county
for(fips in unique(covid_filled_plot_data$fips)) {
  dat <- covid_filled_plot_data[covid_filled_plot_data$fips == fips,] 
  dat %>%
    arrange(date) %>%
    mutate(ma_cases = rollmean(daily_cases, k = 7, fill = NA)) -> dat
  
  dat %>%
    mutate(two_week_diff = ma_cases - dplyr::lag(ma_cases, 14)) %>%
    mutate(two_week_diff_scaled = two_week_diff/max(ma_cases, na.rm = TRUE)) %>% 
    mutate(color = case_when(two_week_diff_scaled < 0 ~ "falling",
                             two_week_diff_scaled < .1 ~ "same",
                             two_week_diff_scaled < .2 ~ "rising_1",
                             two_week_diff_scaled < .3 ~ "rising_2",
                             two_week_diff_scaled >= .3 ~ "rising_3")) -> dat
  
  dat %>%
    filter(date >= start_date) %>%
    ggplot() +
    geom_tile(aes(date, mean(range(daily_cases)), height = diff(range(daily_cases)), fill = color)) +
    geom_point(aes(date, daily_cases), alpha = .5, size = .5) +
    geom_line(aes(date, ma_cases), size = .5) +
    theme_void() +
    ggtitle(paste(unique(dat$county), unique(dat$state), sep = ", ")) +
    theme(legend.position = "none",
          plot.title = element_text(size = 8)) +
    scale_fill_manual(values = c("#e0f3f8", "#ffffbf", "#fee090", "#fc8d59", "#d73027"),
                      breaks = c("falling", "same", "rising_1", "rising_2", "rising_3")) -> p
  
  ## save
  ifelse(!dir.exists(file.path(paste0("examples/covid/covid_daily_plots_mobility_restaurant_home_", start_label, "/"))), dir.create(file.path(paste0("examples/covid/covid_daily_plots_mobility_restaurant_home_", start_label, "/"))), FALSE)
  
  
  ggsave(p, device = "png", filename = paste0("examples/covid/covid_daily_plots_mobility_restaurant_home_", start_label, "/", fips, ".png"), width = 1, height = 1, units = "in", dpi = 200, bg = "white")
}

## add img label to hc object
data.frame(fips = counties) %>% 
  left_join(fips_codes %>% mutate(fips = paste0(state_code, county_code)) %>% select(fips, state, county)) %>%
  mutate(county = gsub(" County", "", county)) %>%
  mutate(label = paste(county, state, sep = ", ")) %>%
  inner_join(data.frame(idx = 1:length(hc$labels), fips = hc$labels)) %>%
  arrange(idx) -> label_df
    
label_df %>%
  pull(fips) %>%
  paste0(".png") -> hc$img 

label_df %>%
  pull(label) -> hc$labels

## save clustering object
save(hc, file = paste0("examples/covid/county_mobility_restaurant_home_", start_label, "-hc.Rdata"))


