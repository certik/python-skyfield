import datetime as dt
import numpy as np
from pytz import timezone
from matplotlib.pylab import plot, savefig, legend, grid, gca
from skyfield.api import load, wgs84

# https://en.wikipedia.org/wiki/Jet_Propulsion_Laboratory_Development_Ephemeris
eph = load('de421.bsp')

earth, sun, moon, venus = eph['earth'], eph['sun'], eph['moon'], eph['venus']
# Los Alamos
observer = wgs84.latlon(35.90369476314685, -106.30851268194805, 2230)

ts = load.timescale()
#t = ts.utc(2024, 3, 30, 16, 30)
zone = timezone('US/Mountain')
now = zone.localize(dt.datetime.now())
t = ts.from_datetime(now)



print("Location: Los Alamos, NM")
print(observer)
print("Time:", t.astimezone(zone))

apparent = (earth + observer).at(t).observe(sun).apparent()
alt, az, distance = apparent.altaz()
print("Sun:")
print("Altitude (0-90):", alt)
print("Azimuth (0-360): ", az)
print()

apparent = (earth + observer).at(t).observe(moon).apparent()
alt, az, distance = apparent.altaz()
print("Moon:")
print("Altitude (0-90):", alt)
print("Azimuth (0-360): ", az)
print()

apparent = (earth + observer).at(t).observe(venus).apparent()
alt, az, distance = apparent.altaz()
print("Venus:")
print("Altitude (0-90):", alt)
print("Azimuth (0-360): ", az)
print()
print()

#---------------

# Fredericksburg, TX
observer = wgs84.latlon(30.274167, -98.871944, 516)

ts = load.timescale()
zone = timezone('US/Central')
now = zone.localize(dt.datetime(2024, 4, 8, 13, 35, 10))
t = ts.from_datetime(now)


print("Location: Fredericksburg, TX")
print(observer)
print("Time:", t.astimezone(zone))

apparent = (earth + observer).at(t).observe(sun).apparent()
alt, az, distance = apparent.altaz()
print("Sun:")
print("Altitude (0-90):", alt)
print("Azimuth (0-360): ", az)
print()

apparent = (earth + observer).at(t).observe(moon).apparent()
alt, az, distance = apparent.altaz()
print("Moon:")
print("Altitude (0-90):", alt)
print("Azimuth (0-360): ", az)
print()

apparent = (earth + observer).at(t).observe(venus).apparent()
alt, az, distance = apparent.altaz()
print("Venus:")
print("Altitude (0-90):", alt)
print("Azimuth (0-360): ", az)
