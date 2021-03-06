# QA script for Solar Position Algorithm


# Example 1: This is the example given in Ibrahim and Afshin (2008) report
timestamp <- as.POSIXct("2003-10-17 12:30:30")
spa(timestamp, timezone=-7, longitude=-105.1786, latitude=39.742476, elevation=1830.14, pressure=820, temperature=11, surface.slope=30, surface.azimuth.rotation=-10, DeltaT=67)


# Example 2: This example uses the latitude and longitude of Cyprus with Daylight saving time (DST) (i.e. month August)
timestamp <- as.POSIXct("2001-8-8 06:30:30")
spa(timestamp, timezone=3, longitude=33.429859, latitude=35.126413, elevation=135.89, pressure=500, temperature=25, surface.slope=30, surface.azimuth.rotation=-10, DeltaT=64.2)


# Example 3: Same as example 2 but without DST (i.e. the timezone is different - month February)
timestamp <- as.POSIXct("2001-2-8 06:30:30")
spa(timestamp, timezone=2, longitude=33.429859, latitude=35.126413, elevation=135.89, pressure=500, temperature=25, surface.slope=30, surface.azimuth.rotation=-10, DeltaT=64.2)


# Example 4: This example uses a non leap year
timestamp <- as.POSIXct("1995-3-1 17:20:15")
spa(timestamp, timezone=2, longitude=-3.7037902, latitude=40.4167754, elevation=647.67, pressure=1000, temperature=19.4, surface.slope=0, surface.azimuth.rotation=180, DeltaT=60.79)


# Example 5: This example uses a leap year
timestamp <- as.POSIXct("2016-10-18 15:00:00")
spa(timestamp, timezone=-3, longitude=-58.381592, latitude=-34.603722, elevation=30, pressure=1010, temperature=18.17, surface.slope=0, surface.azimuth.rotation=180, DeltaT=68.5)

