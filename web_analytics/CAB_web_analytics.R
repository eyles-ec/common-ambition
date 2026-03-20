#load these libraries. if you don't have them use these lines to install them. select them and press CTRL SHIFT C to remove the # (comment indicators)
# install.packages(c(
#   "googleAnalyticsR",
#   "googleAuthR",
#   "dplyr",
#   "readr",
#   "stringr",
#   "lubridate",
#   "ggplot2",
#   "tidyr"
# ))

library(googleAnalyticsR)
library(googleAuthR)
library(dplyr)
library(readr)
library(stringr)
library(lubridate)
library(ggplot2)
library(tidyr)

#small function to transform na into 0, mathematically necessary
as_num_NA <- function(x) ifelse(is.na(x), 0, as.numeric(x))

#collect paths - save your own paths to a file called paths.R that is ignored by git (.gitignore) 
source("../paths.R")

#set working directory 

setwd(wd)

#account authentication

#Launch the browser to create a log in token
#make sure you have the editor role (though I think admin works too)
#you can answer '1' if you want it to cache those credentials, which means you don't have to log in as often
googleAnalyticsR::ga_auth()

#now our R is authorised to use analytics

#data extraction steps
property_id <- property_id   #GA4 numeric property ID, you find this on the page, this is the one for CAB
#in this case it is saved in the paths.R 

#set the dates you're interested in data for. the dates are in Year Month Day format
start_date <- "2024-01-01"
end_date   <- "2026-03-17"

#for dimensions and metrics see here: https://developers.google.com/analytics/devguides/reporting/data/v1/api-schema for more options
# dimensions - what kind of activity. a row in the table is a unique combination of these descriptors.
dims <- c(
  "date",
  "pagePath", #page viewed
  "sessionSource", #where they came from
  "sessionMedium", #type of traffic
  "deviceCategory", #type of device
  "country", #country
  "screenResolution" #resolution of screen, can proxy bots
)

#counts something that occured
metrics <- c(
  "screenPageViews", #page views
  "sessions", #sessions
  "activeUsers", #number of active users
  "newUsers", #first time users
  "engagedSessions", #sessions with engagement
  "engagementRate", #% sessions engaged
  "bounceRate", # 1- engagement rate
  "userEngagementDuration", #total engagement time
  "eventCount", #total no of events in event name
  "averageSessionDuration" #avg session length
)

#fetch raw data for processing 
#each row represents a 'bucket' of activity, like people characteristed by the 'dimensions' (e.g. they came on  a phone in the UK plus some other things)
ga_raw <- ga_data(
  propertyId = property_id,
  date_range = c(start_date, end_date),
  dimensions = dims,
  metrics    = metrics,
  limit      = -1   #auto‑paging
)

#save the raw data
write_csv(ga_raw, "ga4_raw_full.csv")

#separately pull out yuno clicks. this filters the data based on type (click), then also missing links, then matches on the keyword 'yuno'
#you get a row for each of the dimensions, so you may have more than one row per day (but you can see where the came from before clicking!)
yuno_events <- ga_data(
  propertyId = property_id,
  dimensions = c("date", "linkUrl", "eventName", "sessionSource", "sessionMedium"),
  metrics    = "eventCount",
  date_range = c(start_date, end_date),
  limit      = -1
) %>%
  #keep only click-like events
  filter(eventName %in% c("click", "outbound", "external_link_click", "file_download")) %>%
  #guard against NA linkUrl before pattern matching
  mutate(linkUrl = replace_na(linkUrl, "")) %>%
  #case-insensitive match for "yuno" in the URL
  filter(str_detect(linkUrl, regex("yuno", ignore_case = TRUE)))

#save the yuno data
write.csv(yuno_events, "yuno.csv")

#clean out bot traffic (Hopefully)

ga_clean <- ga_raw %>%
  #this step removes any missing values and converts to numeric type 
  mutate(
    screenPageViews        = as_num_NA(screenPageViews),
    sessions               = as_num_NA(sessions),
    activeUsers            = as_num_NA(activeUsers),
    newUsers               = as_num_NA(newUsers),
    engagedSessions        = as_num_NA(engagedSessions),
    engagementRate         = as_num_NA(engagementRate),
    bounceRate             = as_num_NA(bounceRate),
    userEngagementDuration = as_num_NA(userEngagementDuration),
    eventCount             = as_num_NA(eventCount)
  ) %>%
  #require some engagement (>=10s), you can change this to whatever number you like
  filter(userEngagementDuration >= 10) %>%
  #common bot source filters
  filter(!str_detect(sessionSource,
                     regex("bot|crawl|spider|monitor|uptime|externalhit|amazonaws|semrush|ahrefs",
                           ignore_case = TRUE))) %>%
  #synthetic screen sizes (e.g. no one will have a 0x0 screen)
  filter(screenResolution != "0x0") %>%
  #require more than 1 event OR more than 1 page view
  filter(eventCount > 1 | screenPageViews > 1) %>%
  #exclude auto‑detected crawlers
  filter(is.na(deviceCategory) | deviceCategory != "crawlers")

library(dplyr)
library(stringr)
library(tibble)


#filter to UK only
ga_uk <- ga_clean %>%
  filter(country == "United Kingdom")

#save both all countries and uk only data
write_csv(ga_clean, "ga4_clean_all_countries.csv")
write_csv(ga_uk, "ga4_clean_uk.csv")

#note that for below i use only UK sourced 'events'/'views' but if you want it for all just use the ga_clean data

#daily data table for UK 
daily<- ga_uk %>%
  group_by(date) %>%
  summarise(
    pageviews      = sum(screenPageViews, na.rm = TRUE),
    users          = sum(activeUsers, na.rm = TRUE),
    new_users      = sum(newUsers, na.rm = TRUE),
    returning_users = pmax(users - new_users, 0),
    sessions       = sum(sessions, na.rm = TRUE),
    engaged        = sum(engagedSessions, na.rm = TRUE),
    engagement_rate = ifelse(sessions > 0, engaged / sessions, NA_real_),
    bounce_rate     = 1 - engagement_rate,
    pages_per_sess  = ifelse(sessions > 0, pageviews / sessions, NA_real_),
    avg_sess_dur_s  = ifelse(sessions > 0,
                             sum(userEngagementDuration, na.rm = TRUE) / sessions,
                             NA_real_)
  ) %>%
  ungroup() %>%
  arrange(date)

write_csv(daily, "ga4_daily_uk.csv")

#find the views per page on the site for the whole range of dates
top_pages <- ga_uk %>%
  group_by(pagePath) %>%
  summarise(
    pageviews = sum(screenPageViews, na.rm = TRUE),
    users     = sum(activeUsers, na.rm = TRUE),
    sessions  = sum(sessions, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(pageviews))

#save
write_csv(top_pages, "uk_top_pages.csv")


#traffic sources (e.g. where are people coming from) for the whole range of dates
#if it is (direct) that means someone went directly to the CAB link
traffic_sources <- ga_uk %>%
  group_by(sessionSource, sessionMedium) %>%
  summarise(
    users     = sum(activeUsers, na.rm = TRUE),
    sessions  = sum(sessions, na.rm = TRUE),
    pageviews = sum(screenPageViews, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(sessions))

#save it 
write_csv(traffic_sources, "uk_traffic_sources.csv")

#plot of pageviews
ggplot(daily, aes(as.Date(date), pageviews)) +
  geom_line(color = "darkblue") +
  labs(title = "Daily Pageviews (Cleaned)", x = "Date", y = "Pageviews") +
  theme_minimal()

