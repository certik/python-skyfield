import numpy as np
from matplotlib.pylab import plot, savefig, legend, grid
from skyfield.api import load

eph = load('de421.bsp')

#earth, sun, venus = eph['earth'], eph['sun'], eph['venus']
#eph.segments[2].spk_segment.compute_and_differentiate(2414864.4)

# In days, spanning about ~154 years
start_time = 2414864.5
end_time = 2471184.5
N = 100 # steps

planets = [
    # Solar System Barycenter
    "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus",
    "Neptune", "Pluto", "Sun",
    # Earth Barycenter
    "Moon", "Earth",
    # Mercury Barycenter
    "Mercury",
    # Venus Barycenter
    "Venus",
    # Mars Barycenter
    "Mars",
]

for p in range(10):
    xx = []
    yy = []
    for i in range(N):
        time = start_time + (end_time - start_time)*i/(N-1)
        x, y, z = eph.segments[p].spk_segment.compute(time)
        xx.append(x)
        yy.append(y)
    plot(xx, yy, label=planets[p])

legend()
grid()
savefig('solar_system.png')
