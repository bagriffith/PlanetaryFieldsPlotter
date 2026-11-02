"""Calculate position relative to planetary magnetic dipoles for Voyager 2."""
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
import spiceypy as spice


class Flyby:
    """A time interval for Voyager 2 flyby."""
    def __init__(self, planet: str, interval: tuple[str, str]) -> None:
        self.interval = [spice.str2et(x) for x in interval]
        self.planet = planet

    def position(self, steps: int = 256, inertial: bool = True):
        """Calculate a position time series during the flyby.
        
        Args:
            steps: Number of time steps across the interval.
            inertial: Should the result be in an inertial frame? If not, a
                planet centered corrotating frame.
        """
        t = np.linspace(self.interval[0], self.interval[1], steps)
        frame = f"{self.planet.capitalize()} Barycenter"
        rot_frame = "J2000" if inertial else f"IAU_{self.planet.upper()}"
        pos, _ = spice.spkpos("Voyager 2", t, rot_frame, "NONE", frame)
        return t, pos


PLOT = False

spice.furnsh(["data/kernels/naif0012.tls",
              "data/kernels/Voyager_2.m05016u.merged.bsp",
              "data/kernels/pck00011.tpc"])

FLYBYS = [Flyby("Saturn", ("1981-08-24T15:24:05Z", "1981-08-27T15:24:05Z")),
          Flyby("Uranus", ("1986-01-24T05:58:24Z", "1986-01-25T05:58:24Z"))]

# In km, from Wikipedia
PLANET_R = {'Saturn': 60_268,
            'Uranus': 25_362}

# Geomagnetic Pole, (Latitude, Longitude, Offset(X, Y, Z)) in deg
PLANET_MAG = {
    'Saturn': (89.3, 331.,
               [0., 0., 0.]),  # Ness et al., Science 212, 211 (1981).
    'Uranus': (30., -48. + 150,
               [0., 0., -0.31])  # Ness et al., Science 233, 85 (1986).
                }


def main():
    for fb in FLYBYS:
        # Calculate the position relative to the planet in R_planet units.
        times, positions = fb.position(steps=1000, inertial=False)
        positions -= np.repeat(
            PLANET_R[fb.planet] * \
                np.array(PLANET_MAG[fb.planet][2])[np.newaxis, :],
            positions.shape[0],
            0
        )

        # Rotate the frame to align with the magnetic field.
        m = Rotation.from_euler('ZY',
                                np.deg2rad([[PLANET_MAG[fb.planet][1],
                                            90. - PLANET_MAG[fb.planet][0]]])
                                ).as_matrix()[0].T
        positions = (m @ positions.T).T / PLANET_R[fb.planet]

        # Calculate parameters
        r = np.linalg.norm(positions, axis=1)
        mag_lat = np.arcsin(positions[:, 2] / r)
        l_shell = r / np.cos(mag_lat)**2

        # Save the result
        np.savetxt(f'data/{fb.planet}.csv',
                   np.concatenate([positions,
                                   mag_lat[:, np.newaxis],
                                   l_shell[:, np.newaxis]],
                                  axis=1),
                   fmt='%3.3f',
                   delimiter=',', header='x,y,z,latitude,l-shell')

        if PLOT:
            # Some diagnostic plots used for development/testing.
            # Path side profile
            plt.axis('equal')
            plt.plot(np.linalg.norm(positions[:, :2], axis=1), positions[:, 2])
            plt.show()
            plt.cla()

            # L-shell and magnetic latitude
            plt.plot(times, l_shell)
            plt.twinx().plot(times, np.rad2deg(mag_lat), c='C1')
            plt.show()
            plt.cla()

            fig = plt.figure(figsize=(9, 9))
            ax = fig.add_subplot(111, projection="3d")
            ax.plot(positions[:, 0], positions[:, 1], positions[:, 2])
            ax.scatter([0.], [0.], [0.], c='k')  # Planet Marker
            plt.show()
            plt.cla()


if __name__ == '__main__':
    main()
