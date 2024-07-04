# Ideal Gas Model

A simple simulation of thermal behaviors of ideal gas molecules.

This model simulates collisions between large number of small spheres (each representing a gas molecule). This model aims to prove the Maxwell distribution of molecule speed.

## Run

Clone the repository and run the following command:

```shell
cargo run --release
```

You should be able to see the simulation started in full screen.

With the current level of optimization, the model can simulate up to around 8,000 molecules at the same time (with FPS around 40). Further optimization is possible.
