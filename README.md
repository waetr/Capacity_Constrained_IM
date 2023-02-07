# Capacity Constrained IM

Environment: Windows/Linux with $c++14$, $cmake$ and $OpenMP$.

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
./exp xxx [--k aa --eps bb --greedy cc]
```

**xxx**: Name of the graph file in /data/;

**--k aa**: Set the constant $k$ as $aa$. (default as 10)

**--eps bb**: Set the parameter $\epsilon$ of solutions extending *OPIM-C* as $bb$. (default as 0.1)

**--greedy cc**: Specify the type of MC simulation-based solutions as $cc$ to execute. Possible values:

+ **0**: Do not execute these solutions.
+ **1** (default): Execute these solutions with CELF trick.
+ **2**: Execute these solutions with vanilla version.

#### Example

```bash
dnc-corecipient --k 10 --eps 0.05 --greedy 1
```

Run the experiments on *DNC* dataset, with $k=10,\epsilon=0.05$, with the CELF greedy solutions.

## More Experiments

If you want to execute more experiments mentioned in our paper,

edit the following statement in CMakeLists.txt.

```cmake
add_executable(exp src/exp_batch.cpp)
```

**exp_batch.cpp**: (default) Main experiment.

**exp_query_set.cpp**: Randomly generate one query set.

**exp_gamma.cpp**: Calculate the curvature for a graph and a query set.

