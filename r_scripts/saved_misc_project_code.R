
Attempting an animated plot. Doesn't seem to work on my PC for some reasons which I think has to do with the ImageMagick software. I'll figure it out later. 

```{r animated wash smoke map, message=F, warning=F, results="asis", eval = F}

# data wrangle ----
wa_smk_days <- wa_ts %>% 
  select(county, date, smoke5, smoke10, smoke15, smoke_wave) %>% 
  # county needs to be lowercase to join to maps dataframe
  mutate(county = tolower(county)) %>% 
  # need to gather smoke variables by date
  gather(key = smk_cut_method, value = smoke_days, -date, -county) %>% 
  # maybe too many dates, going to filter to just september to october
  filter(date >= "2012-09-15" & date <= "2012-10-01")
  
summary(wa_smk_days)
# join county counts to map df
wa_smk_df <- wa_county_df %>% 
  right_join(wa_smk_days, by = c("subregion" = "county")) %>% 
  filter(smk_cut_method != "smk_day0") %>% 
  # rename smoke day cutoffs 
  mutate(smk_cut_method = 
    ifelse(smk_cut_method == "smoke5","Smoke PM2.5 > 5 ug/m^3",
    ifelse(smk_cut_method == "smoke10", "Smoke PM2.5 > 10 ug/m^3",
    ifelse(smk_cut_method == "smoke15", "Smoke PM2.5 > 15 ug/m^3",
    ifelse(smk_cut_method == "smoke_wave", "Smoke Wave", NA)))),
    # preserve smoke method order for small multiple
    smk_cut_method = factor(smk_cut_method, 
      levels = c("Smoke PM2.5 > 5 ug/m^3", "Smoke PM2.5 > 10 ug/m^3", 
                 "Smoke PM2.5 > 15 ug/m^3", "Smoke Wave")))

# map ----

# smoke days where GWR smoke > 15
smoke_map <- ggplot(wa_smk_df, aes(x=long, y=lat, group=group, frame = date)) +
  # fill with number of smoke days
  geom_polygon(aes(fill = smoke_days), alpha = 0.7) +
  scale_fill_gradient("Smoke-Impacted Days",
                      low = "#0b486b", high = "#f56217") +
  # add county path on top
  geom_path() +
  facet_wrap(~smk_cut_method) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(colour=NA, fill=NA))

gganimate(smoke_map, "output.gif")

ggani
ani.options()

magickPath <-shortPathName("C:\\Program Files\\ImageMagick-7.0.5-Q16\\magick.exe")
ani.options(convert=magickPath)
```