import numpy as np
from matplotlib.pylab import plot, savefig, legend, grid, gca
from scipy.spatial.transform import Rotation
from skyfield.api import load

eph = load('de421.bsp')

#earth, sun, venus = eph['earth'], eph['sun'], eph['venus']
#eph.segments[2].spk_segment.compute_and_differentiate(2414864.4)

# In days, spanning about ~154 years (1900 to 2050)
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
