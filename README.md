# Capacity Constrained IM



## Compile and Run

#### Compile [using Linux bash]

```bash
cd exp
mkdir build && cd build
cmake ..
make
```

#### Run

```bash
./exp xxx [-r 10000]
```

**xxx**: Name of the graph file in /data/;

**-r 10000**: (optional) Conduct the number of MC simulations (with 10000).



## Experiments

Selecting experiment: Edit the following statement in CMakeLists.txt.

```cmake
add_executable(exp src/exp_batch.cpp)
```

**exp_batch.cpp**: (default) Main experiment.

**exp_query_set.cpp**: Randomly generate one query set.

**exp_gamma.cpp**: Calculate the curvature for a graph and a query set.

