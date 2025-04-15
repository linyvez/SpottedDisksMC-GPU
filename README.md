# Modeling of Patchy Disks and Square Plates with GPU-Acceleration using Monte Carlo Method

### Authors

https://github.com/bilyi1andrii \
https://github.com/hhafiya \
https://github.com/linyvez \
https://github.com/shshrg

### Prerequisites

make, cmake, CUDA, gcc

### Idea

This program is meant to model the movement of Patchy Disks and Square Plates in a 2D space. Patchy particles are
particles with surface regions (patches) that have distinct interaction, leading to self-assembly behaviours. Such
models may serve useful in different fields, like biology, physics and chemistry.

### Installation and Usage

To use the program, please follow these steps:

1. Download the contents of this repository
2. Compile the program using the following command:

  ```
  ./compile.sh [OPTIONS]
  ```

Where `[OPTIONS]` specify compilation type:

- `-d` or `--debug-build`
- `-o` or `--optimize-build`
- `-i` or `--relwithdebinfo-build`

  Other:
- `-h` or `--help`
- `c` or `--clean`

3. Run it using the following command:

```
   ./particle_simulation [periodic condition option] [particle type] [patch type]
   ```

Where:

- `periodic condition option`: type of periodic boundary condition that will be used:
  - 0 - no boundary condition,
  - 1 - condition on the X-axis, 
  - 2 - on the Y-axis, 
  - 3 - condition on both X and Y axes);
- `particle type`: 
  - 0 - circular particles,
  - 1 - square particles
- `patch type`: type of the patches

For example:

  ```
  ./particle_simulation 0 0 0
  ```

### Results

As a result of this program, you will get an .xyz file, which can be found in the `data` folder.

### Changing the parameters

If you would like to change some parameters of the simulation, you can change some values in the header files inside the
`include` folder. The most useful of those can be found in `general_config.h` :

- `N`: number of particles inside the simulation cell
- `Lx`: width of the cell
- `Ly`: height of the cell
- `ANIMATION_STEPS`: number of frames
