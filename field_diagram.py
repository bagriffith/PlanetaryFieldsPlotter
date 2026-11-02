"""Plot Voyager 2 Flyby Magnetic Fields."""
from typing import Optional
import numpy as np
from numpy.typing import ArrayLike
from matplotlib.axes import Axes
from matplotlib.patches import Circle, Annulus
from matplotlib.transforms import Affine2D
import matplotlib.pyplot as plt


class FlybyFieldDiagram(Axes):
    """An interval of time where Voyager 2 had a flyby with a named planet."""
    def plot_flyby(self, path_xy: ArrayLike, path_z: ArrayLike,
                   field_angle: float = 0, rot_angle: float = 0,
                   field_l_shells: Optional[list[float]] = None,
                   show_lines: bool = False, show_vectors: bool = True,
                   flip_field: bool = False, **kwargs):
        """Plot the trajectory, magnetic field vectors, and field lines."""
        # First, rotate everything drawn so it has the correct angle from
        # vertical. Useful to show Uranus is highly tilted.
        trans_rot =  Affine2D().rotate_around(0, 0, -np.deg2rad(rot_angle))

        trans_field = trans_rot + Affine2D()\
            .rotate_around(0, 0, np.deg2rad(field_angle))\
            .scale(-1 if flip_field else 1, 1) # Flip to plot on the left side

        # Add the regular Axes transformations.
        trans_rot = trans_rot + self.transData
        trans_field = trans_field + self.transData

        # Plot the trajectory
        self.plot(path_xy, path_z,
                  lw=kwargs.pop('path_lw', 2),
                  color=kwargs.pop('path_c', 'C0'),
                  transform=trans_field)

        self.add_planet(planet_c=kwargs.pop('planet_c', None),
                        transform=trans_rot)

        if field_l_shells is None:
            # Line density decreases with radius, visualizing field strength.
            field_l_shells = \
                list(-0.1 + np.linspace(0, 3**(-1/3),
                                 kwargs.pop('n_lines', 16) + 1)[1:]**(-3))

        if show_vectors:
            # Vectors are drawn at L shell crossings so they are tangent.
            self.plot_vectors(*self._l_crossings(path_xy, path_z,
                                                 field_l_shells),
                              length=kwargs.pop("vector_length", 2.0),
                              width=kwargs.pop("vector_width", 0.1),
                              head_width=kwargs.pop("vector_head", 0.5),
                              fc=kwargs.pop('vector_c', 'C2'),
                              overhang=kwargs.pop('vector_overhang', 0.25),
                              transform=trans_field)

        if show_lines:
            self.plot_field_lines(field_l_shells,
                                  field_lc=kwargs.pop('field_lc', None),
                                  transform=trans_field)

    def add_planet(self, planet_c='C1',
                   pole_length=4, pole_width=0.25, pole_c='k',
                   **kwargs):
        # Planet
        planet = Circle((0,0), 1., color=planet_c, **kwargs)
        planet.set_zorder(2.1)
        self.add_patch(planet)

        # Rotation Axis
        self.arrow(0, -7*pole_length/16, 0, pole_length,
                   zorder=2.0, overhang=0, width=pole_width, lw=0,
                   length_includes_head=True,
                   head_width=3*pole_width, head_length=3*pole_width,
                   fc=pole_c, **kwargs)

        self.text(0, pole_length/2, '\u21bb',
                  transform_rotates_text=True,
                  fontsize='xx-large',
                  verticalalignment='center',
                  horizontalalignment='center',
                  c=pole_c,
                  **kwargs)

    def plot_vectors(self, xy, z, length: float = 1.0, **kwargs):
        for xy_i, z_i in zip(xy, z):
            b_xy = 3*xy_i*z_i
            r2 = xy_i**2 + z_i**2
            b_z = 3*z_i**2 - r2
            b = np.sqrt(b_xy**2 + b_z**2)

            self.arrow(xy_i, z_i, length*b_xy/b, length*b_z/b,
                       zorder=2.0, length_includes_head=True,
                       shape='full', lw=0, **kwargs)

    def plot_field_lines(self, l_shells: list[float], points: int = 256,
                         field_lc='k', field_lw=0.25, **kwargs):
        field_line_kwargs = {'c': field_lc,
                             'lw': field_lw,
                             'zorder': 1.9,
                             **kwargs}
        
        for l_i in l_shells:
            path = self._l_shell_path(l_i, points)
            self.plot(path[0], path[1], **field_line_kwargs)   # Right side
            self.plot(-path[0], path[1], **field_line_kwargs)  # Left Side

        # Field lines at the poles
        self.plot(np.zeros(points),
                  np.linspace(1, np.max(l_shells), points),
                  **field_line_kwargs)
        self.plot(np.zeros(points),
                  -1*np.linspace(1, np.max(l_shells), points),
                  **field_line_kwargs)

    @staticmethod
    def _l_crossings(path_xy: ArrayLike, path_z: ArrayLike,
                     l_shells: list[float]):
        path_xy = np.array(path_xy)
        path_z = np.array(path_z)

        # Find all crossings
        path_l = (path_xy**2 + path_z**2)**(3/2) / path_xy**2

        crossings = []

        for l_i in l_shells:
            sign = np.sign(path_l - l_i)
            change = np.arange(len(sign))[np.roll(sign, 1) - sign != 0]
            change = change[change != len(sign) - 1] # Last element can't be

            for i in change:
                fraction = (l_i - path_l[i]) / (path_l[i+1] - path_l[i])
                crossings.append(
                    [a[i] + fraction*(a[i+1] - a[i]) for a in [path_xy, path_z]]
                )

        crossings = np.array(crossings)

        return crossings[:, 0], crossings[:, 1]

    @staticmethod
    def _l_shell_path(l_shell: float, points: int = 256):
        limit_theta = np.arccos(np.sqrt(1/l_shell))
        theta_l = np.linspace(-limit_theta, limit_theta, points)
        r = l_shell * np.cos(theta_l)**2
        return np.array([r*np.cos(theta_l), r*np.sin(theta_l)])


FLYBY_COLORS= {'path_c': 'gray',
               'vector_c': 'black',
               'planet_c': 'lightgray',
               'field_lc': 'dimgray'}


def saturn():
    """Plot the Saturn flyby."""
    fig = plt.figure(figsize=(4, 4))
    ax = FlybyFieldDiagram(fig, (0.1, 0.1, 0.8, 0.8))
    fig.add_axes(ax)

    ax.set_aspect('equal')
    ax.set_xlim(-5.5, 20.5)
    ax.set_ylim(-13, 13)

    data = np.loadtxt('data/Saturn.csv', delimiter=',')
    ax.plot_flyby(path_xy=np.linalg.norm(data[:, :2], axis=1),
                  path_z=data[:, 2],
                  show_lines=True,
                  field_angle=0.2,
                  rot_angle=0,
                  **FLYBY_COLORS)

    # Ring
    ring = Annulus((0, 0), (2, 0.65), 0.4,
                   color=FLYBY_COLORS['planet_c'])
    ring.set_zorder(2.2)
    ax.add_patch(ring)

    fig.savefig('Saturn.pdf')


def uranus(answer=False):
    """Plot the Uranus flyby."""
    fig = plt.figure(figsize=(4, 4))
    ax = FlybyFieldDiagram(fig, (0.1, 0.1, 0.8, 0.8))
    fig.add_axes(ax)

    ax.set_aspect('equal')
    ax.set_xlim(-8, 18)
    ax.set_ylim(-15.5, 10.5)

    data = np.loadtxt('data/Uranus.csv', delimiter=',')
    ax.plot_flyby(path_xy=np.linalg.norm(data[:, :2], axis=1),
                  path_z=data[:, 2],
                  field_angle=60.,
                  rot_angle=82.,
                  n_lines=20,
                  show_lines=answer,
                  **FLYBY_COLORS)

    if answer:
        # Mark the equator and magnetic poles
        ax.plot([-20, 20], [8.5, -8.5], ls='--', c='r', zorder=3)
        ax.plot([-0.75, 0.75], [-2, 2], lw=4, c='r', zorder=3)
        ax.text(1, 3, 'N', fontsize='xx-large', c='r', zorder=3)
    else:
        # Show arrow to the sun to illustrate the highly tilted axis.
        ax.arrow(0.1, 0.88, 0.15, 0, transform=ax.transAxes,
                width=.002, head_width=.02, fc='k')
        ax.annotate('To Sun', (0.1, 0.88), xycoords='axes fraction',
                    verticalalignment='bottom')

    fig.savefig(f"Uranus{'-key' if answer else ''}.pdf")


if __name__ == '__main__':
    saturn()
    uranus(answer=False)
    uranus(answer=True)
