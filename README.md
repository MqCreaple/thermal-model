# Ideal Gas Model

A simple simulation of thermal behaviors of ideal gas molecules.

This model simulates collisions between large number of small spheres (each representing a gas molecule). This model aims to prove the Maxwell distribution of molecule speed.

## Run

Clone the repository and run the following command:

```shell
cargo run --release
```

You should be able to see the simulation started in full screen.

You can also run different examples with command:

```shell
cargo run --release --example <example-name>
```

With the current level of optimization, the model can simulate up to around 80,000 molecules with multithreading or 10,000 molecules without multithreading (with FPS around 30 on Intel Core i9). Further optimization is possible.
