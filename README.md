# BIM-Project



## Datasets

| dataset                                                      | type       | N         | M           |
| ------------------------------------------------------------ | ---------- | --------- | ----------- |
| [soc-Pokec](https://snap.stanford.edu/data/soc-Pokec.html)   | Directed   | 1,632,803 | 30,622,564  |
| [soc-LiveJournal1](https://snap.stanford.edu/data/soc-LiveJournal1.html) | Directed   | 4,847,571 | 68,993,773  |
| [com-Orkut](https://snap.stanford.edu/data/com-Orkut.html)   | Undirected | 3,072,441 | 117,185,083 |
| [DNC-corecipient](http://konect.cc/networks/dnc-corecipient/) | Undirected | 906       | 12,085      |

No isolated nodes



## Compile and Run

Compile

```bash
cd exp
mkdir build
cd build
cmake ..
make
```

Run

```
./exp xxx.txt [-r 10000]
```

xxx.txt is the graph file in /data/;

-r is followed by the number of MC simulations.