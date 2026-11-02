# Voyager 2 Flyby Plotter

## Setup

First, install the required packages. Using a new venv is recommended.

```bash
python install -r requirements.txt
```

## Voyager Position

Next, calculate Voyager 2's trajectory in magnetic coordinates for each flyby.
These are calculated from SPICE kernels, which must be downloaded into a
`data/kernels` directory.

- [naif0012.tls](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls)
- [pck00011.tpc](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc)
- [Voyager_2.m05016u.merged.bsp](https://naif.jpl.nasa.gov/pub/naif/VOYAGER/kernels/spk/Voyager_2.m05016u.merged.bsp)

Then, from the directory containing `data`, run the position script.

```bash
python voyager_position.py
```

## Magnetic Field Diagrams

Finally, create the magnetic field diagrams. It will output three files,
`Saturn.pdf`, `Uranus.pdf`, and `Uranus-key.pdf`.

```bash
python 
```
