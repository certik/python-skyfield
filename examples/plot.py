import numpy as np
from matplotlib.pylab import plot, savefig, legend, grid, gca
from scipy.spatial.transform import Rotation
from skyfield.api import load
from skyfield.timelib import Time

# https://en.wikipedia.org/wiki/Jet_Propulsion_Laboratory_Development_Ephemeris
eph = load('de421.bsp')

# In days, spanning about ~154 years (1899 to 2053)
ts = load.timescale()
start_time = 2414864.5
end_time = 2471184.5
# To convert from a date: time = ts.utc(2024, 4, 1)
print("start_time:", Time(ts, start_time).utc)
print("end_time:", Time(ts, end_time).utc)

N = 1000 # steps
AU = 149.5978707e6 # km

planets = [
    # Solar System Barycenter
    "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus",
    "Neptune", "Pluto", "Sun",
    # Earth Barycenter
    "Moon", "Earth",
]

# The system of coordinates seems to be with Earth's north pole pointing up
# along the z-axis. We rotate along the x-axis clock-wise by the Earth's axial
# tilt of 23.44° to get the Earth's orbit into the (x,y) plane.
R = Rotation.from_euler('x', [-23.44], degrees=True).as_matrix()

Np = 10

data = np.empty((Np, N, 3))
for p in range(Np):
    xx = []
    yy = []
    for i in range(N):
        time = start_time + (end_time - start_time)*i/(N-1)
        # Input: [days]
        # Output: [km]
        # To get derivatives: compute_and_differentiate(time)
        x, y, z = eph.segments[p].spk_segment.compute(time)
        X = np.array([x, y, z])
        data[p, i, :] = np.dot(R, X)

data = data / AU

for p in range(Np):
    x = data[p,:,0]
    y = data[p,:,1]
    plot(x, y, "-", label=planets[p])

legend()
grid()
gca().set_aspect("equal")
savefig('solar_system.png')
