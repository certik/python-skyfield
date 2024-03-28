import numpy as np
from matplotlib.pylab import plot, savefig, legend, grid, gca
from skyfield.api import load

eph = load('de421.bsp')

#earth, sun, venus = eph['earth'], eph['sun'], eph['venus']
#eph.segments[2].spk_segment.compute_and_differentiate(2414864.4)

# In days, spanning about ~154 years
start_time = 2414864.5
end_time = 2471184.5
N = 1000 # steps
AU = 149.5978707e6 # km

planets = [
    # Solar System Barycenter
    "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus",
    "Neptune", "Pluto", "Sun",
    # Earth Barycenter
    "Moon", "Earth",
]

Np = 7

data = np.empty((Np, N, 3))
for p in range(Np):
    xx = []
    yy = []
    for i in range(N):
        time = start_time + (end_time - start_time)*i/(N-1)
        # Input: [days]
        # Output: [km]
        x, y, z = eph.segments[p].spk_segment.compute(time)
        data[p, i, :] = [x, y, z]

data = data / AU

for p in range(Np):
    x = data[p,:,1]
    y = data[p,:,2]
    plot(x, y, "-", label=planets[p])

legend()
grid()
gca().set_aspect("equal")
savefig('solar_system.png')
