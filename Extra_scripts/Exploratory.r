# Data project scripts for Limnology.
# This script hosts exploratory code for Lake Pulse data.


# Temperature and chlorophyll

p1 <- ggplot(plankton, aes(x = Temperature_C)) + geom_histogram() # Bimodal

p2 <- ggplot(plankton, aes(x = Chlorophylla_ugL)) + geom_histogram() # Left skew

p1 + p2

ggplot(filter(plankton, Chlorophylla_ugL < 20),
    aes(x = Chlorophylla_ugL)) +
    geom_histogram()

ggplot(plankton, aes(x = log(Chlorophylla_ugL))) + geom_histogram()

ggplot(plankton, aes(x = Temperature_C, y = log(Chlorophylla_ugL))) +
    geom_point() +
    geom_smooth()

# Nutrient concentrations also left skewed

p1 <- ggplot(filter(plankton, Total_phosphorus_ugL < 1800),
    aes(x = Total_phosphorus_ugL)) +
    geom_histogram()
p2 <- ggplot(filter(plankton, Total_nitrogen_mgL < 50),
    aes(x = Total_nitrogen_mgL)) +
    geom_histogram()

p1 + p2

p1 <- ggplot(filter(plankton, Total_phosphorus_ugL < 1800),
    aes(x = log(Total_phosphorus_ugL))) +
    geom_histogram()
p2 <- ggplot(filter(plankton, Total_nitrogen_mgL < 50),
    aes(x = log(Total_nitrogen_mgL))) +
    geom_histogram()

p1 + p2

# Nutrients to chlorophyll

p1 <- ggplot(plankton, aes(x = log(Total_nitrogen_mgL), y = log(Chlorophylla_ugL))) +
    geom_point(alpha = 0.5) +
    geom_smooth()

p2 <- ggplot(plankton, aes(x = log(Total_phosphorus_ugL), y = log(Chlorophylla_ugL))) +
    geom_point(alpha = 0.5) +
    geom_smooth()

p1 + p2

# Latitude

ggplot(plankton, aes(x = latitude)) + geom_histogram()

# Diversity

ggplot(plankton, aes(x = zooplankton_richness)) + geom_histogram(bins = 10)
ggplot(plankton, aes(x = rbr_temperature_mean_tube)) + geom_histogram(bins = 10)

ggplot(plankton, aes(x = rbr_temperature_mean_tube, y = zooplankton_richness)) +
    geom_point(alpha = 0.8) +
    geom_smooth()

ggplot(plankton, aes(x = log(chla), y = zooplankton_richness)) +
    geom_point(alpha = 0.8) +
    geom_smooth()

ggplot(plankton, aes(x = zooplankton_richness, y = log(chla))) +
    geom_point(alpha = 0.8) +
    geom_smooth()


# Scratch ----
# Duplicate measures in the same lake-year combinations.

lakepulse_subset |>
    distinct(
        location_name, latitude, longitude,
        CharacteristicName, year, ResultValue
    ) |>
    count(location_name, latitude, longitude) |>
    arrange(desc(n))

Acton <- lakepulse |> filter(MonitoringLocationName == "Acton Lake")
View(Acton)

# Old lakepulse data

lakepulse <-
    read_csv(here("lakepulse.csv")) |>
    select(! c(DatasetName, MonitoringLocationType))

count(lakepulse, MonitoringLocationID) |> arrange(MonitoringLocationID)
count(lakepulse, MonitoringLocationID, ActivityType)

count(lakepulse, CharacteristicName) |> arrange(desc(n)) |> print(n = 27)
count(lakepulse, CharacteristicName, ResultUnit) |> arrange(desc(n))
count(lakepulse, CharacteristicName, ResultUnit) |>
    filter(!is.na(ResultUnit)) |>
    arrange(desc(n)) |>
    print(n = 27)

# Choose measures to pivot into wider format.

lakepulse_subset <- lakepulse |>
    filter(CharacteristicName %in% c(
        "Chlorophyll a, corrected for pheophytin",
        "Temperature, water",
        "Total Nitrogen, mixed forms",
        "Total Phosphorus, mixed forms"
    )) |>
    select(
        MonitoringLocationID:MonitoringLocationLongitude,
        ActivityStartDate, ActivityDepthHeightMeasure,
        CharacteristicName, ResultValue
    ) |>
    mutate(
        CharacteristicName = case_when(
            CharacteristicName == "Temperature, water" ~ "Temperature_C",
            CharacteristicName == "Chlorophyll a, corrected for pheophytin" ~ "Chlorophylla_ugL",
            CharacteristicName == "Total Nitrogen, mixed forms" ~ "Total_nitrogen_mgL",
            CharacteristicName == "Total Phosphorus, mixed forms" ~ "Total_phosphorus_ugL",
            .default = CharacteristicName),
        year = year(ActivityStartDate), .before = ActivityStartDate
    ) |>
    rename(
        location_id = MonitoringLocationID,
        location_name = MonitoringLocationName,
        longitude = MonitoringLocationLongitude,
        latitude = MonitoringLocationLatitude,
        collection_date = ActivityStartDate,
        collection_depth_m = ActivityDepthHeightMeasure
    )

# Clean to one row per lake.

lakepulse_subset |>
    count(location_id, CharacteristicName, collection_depth_m) |>
    print(n = 50)

# Chlorophyll a and temperature have double values, take mean for now

lakepulse_subset <- lakepulse_subset |>
    mutate(
        mean_result = mean(ResultValue),
        .by = c("location_id", "CharacteristicName")
    ) |>
    distinct(location_id, CharacteristicName, mean_result, .keep_all = TRUE) |>
    select(!c(ResultValue, collection_depth_m))

plankton <-
    pivot_wider(
        lakepulse_subset,
        id_cols = location_id:collection_date,
        names_from = CharacteristicName,
        values_from = mean_result
    )

plankton