library(scales)

###0. data processing
#load data
raw_global_FoI_data <- read.csv("input.csv")
raw_global_FoI_data$date <- as.Date(raw_global_FoI_data$date)
raw_global_FoI_data$date_import <- as.Date(raw_global_FoI_data$date_import)

global_FoI_data <- raw_global_FoI_data %>%
  filter(censoring==0) #%>%
  #filter(F_i>0)

#select locations to be visualzied
#sampled_location <- sample(levels(factor(global_FoI_data$location)),
#       size=10,
#       replace = F)
#selected_location <- c("United Kingdom","Netherlands","United States of America", "Japan", "Korea, Republic of")
selected_location2 <- c("Austria", "Bahrain", "Canada", "Chile", "Cuba", "Cyprus", "Guam", "Iceland", "Paraguay", "Peru", "United Kingdom","Netherlands","United States of America", "Japan", "Korea, Republic of")

vis_data <- global_FoI_data %>%
  filter(sub_region %in% c("Western Europe","Northern America"))
  #filter(region %in% c("Europe"))
  #filter(location %in% c(selected_location2))#c(selected_location, sampled_location))

###1. FoI over time
ggplot(data = vis_data,
       mapping = aes(x = date, y =location)) +
  geom_raster(aes(fill=F_i)) +
  scale_fill_gradient(low="grey90", high="red", name="Force of infection") +
  geom_point(mapping = aes(x=date_import,
                           y=location))+#, color=region))+
  theme_bw()+
  labs(x="Time [day]", y="Location")
