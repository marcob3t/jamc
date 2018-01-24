# Tasklist

### Preliminary tasks (25/01 - 1.00 PM)

- Update makefile
- Unit tests
  - compute force for a few 2- or 3 particle systems with atoms inside/outside the cutoff;
  - compute one full step time integration for given forces and velocities (no call to force());
  - compute kinetic energy for given velocities and mass;
  - verify that input parameter data is read correctly;
- Update `.travis.yml`

### Individual tasks (26/01 - 6.00 PM)

- Python interface;
- optimize force computation;
- MPI;
- OpenMP;

### Optional tasks

- Morse potential;
- cell list;
- force kernel in CUDA.
