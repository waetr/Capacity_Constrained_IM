# BIM-project

## 实验用数据集

| Name                                                         | Type       | Nodes     | Edges      |
| ------------------------------------------------------------ | ---------- | --------- | ---------- |
| BlogCatalog                                                  | Undirected | 10,312    | 333,983    |
| [com-DBLP](http://snap.stanford.edu/data/com-DBLP.html)      | Undirected | 317,080   | 1,049,866  |
| [soc-LiveJournal1](http://snap.stanford.edu/data/soc-LiveJournal1.html) | Directed   | 4,847,571 | 68,993,773 |
| [soc-Pokec](http://snap.stanford.edu/data/soc-Pokec.html)    | Directed   | 1,632,803 | 30,622,564 |
| hepPh                                                        | Directed   | 12008     | 237010     |
| facebook                                                     | Directed   | 4039      | 88234      |
| viki-Vote                                                    | Directed   | 7115      | 103689     |

## 编译与运行

需要修改的内容：models.h的MAX_NODE_SIZE

请根据需要修改为数据集可能的最大节点数

```bash
cd exp
mkdir build && cd build
cmake ..
make
./exp com-dblp.csv [-? | -v | -l dblp_mg.txt | -r 10000]
```

Note: 

com-dblp.csv为data文件夹内图文件的文件名；

-l后参数dblp_mg.txt为data文件夹内local single spread文件的文件名；

-r后参数为MC simulation的迭代次数；

-v为打开verbose flag.



## debug日志

2022/8/30 ~ 2022/9/4

问题：用IC model的RR set估计spread时，结果总比MC simulation低5.0~10.0左右

最终推测原因：随机数生成器性质不够好（造成系统误差）。

最终解决方案：在生成random RR set的选取起始点$v$过程，将$v$的随机数生成器从线性同余法（minstd_rand）改为梅森旋转法（mt19937）。对每条边取样时还是选用线性同余法。结果：在保证运行速度不下降的同时，误差显著降低。

## 更新的内容

重构了部分代码。

新增OPIM-Local，OPIM-greedy方法。


## 实验大纲

一个framework,描述option1和option2

输入为：graph, diffusion model, simulation method, set S

输出为：S的邻居节点的子集

graph: 一个图

diffusion model: IC

simulation method: MC, sketch

在option2时需要实现的IM算法：degree, pgrank, CELF

## 更新日志

2022/7/17

1. 写了estimate随机选点的邻居重合度的API。$重合程度=\frac{被两个及以上源点选中的点}{所有被选中的点}$
2. 大致比较了不同的随机数生成器的性能差异：minstd_rand（线性同余法）的速度比mt19937（马特赛特旋转）快约30%。最后选用了minstd_rand生成随机数。*原因：根据线性同余法的性质，在[0,1)范围内近似为平均分布，已经满足simulation的需要。*
3. 改善随机数并且使用-O3加速之后，MC simulation快的飞起，而且r=100和r=1000000的性能差距不大（意味着可以用r=10000正常实验了）
4. 重写了CELF。（原代码在遇到问题(*)时会出bug）

2022/7/18


写了option 2的advance版本。

算法：

1. 把active paticipate set $A$的所有邻居拿出来，记为$\N_A$。*对于每个$v\in \N_A$, $A_v$表示有邻居为$v$的active paticipate的集合。*
2. 按照某种顺序（degree, pageRank, margimal influence(CELF)）对于每个$v\in \N_A$, 尝试将$v$加入到seed里。
3. 每次尝试，任意选择一个$u\in A_v$，标记$u$多选了一个邻居，并且把$v$加入seed里。如果$u$已经选满$k$个邻居了，就把$u$从$A_v$中删除，再从$A_v-\{u\}$中任选一个作为$u$再次尝试，直到$u$还没选满或者$A_v$为空为止。如果$A_v$为空，说明$v$没有对应的active paticipate满足条件了，于是放弃$v$。
4. 把所有$v\in \N_A$依次尝试加入，算法结束。

选用了某一组$size=3$，邻居重合率$=\frac{3}{7}$的paticipate set，$k=2$，在$IC$和$TC-IC$模型上分别测试，输出见test.log。

结论：option 1结果最准确。在邻居重合率高的情况下，option 2漏掉了很多节点，而advance method选的节点数和option 1一样多，且degree,CELF选择的seed和option 1一样。

2022/7/20

修改了一些可能的bug，稍微优化了一些性能有关的细节。

写了顶层experiment方法，目前可以设置好参数直接挂到服务器上跑，输出结果重定向到output/result.out中。

result0.out为上一次本地跑出的部分实验结果。

7/22 2:19

更正了一些逻辑上的问题(CELF advance的while条件，MC simulation的if(deadline=0)部分)，删除了一些冗余变量，增加了一些注释。

修改了Run_simulation方法，使其能传入IM solvers.

多输出了一些内容，便于复现实验。

Problems and Inquiry

**问题：现在的实验方法无法在reasonable amount of time内跑出participant size > 5 且 k > 5的实验数据。**

从早上八点左右截至17：49，实验在服务器上还在size=5部分跑。剩下的size还有{10,20,50,100,200,500}。



观察实验，发现时间基本被CELF以及其advance部分占去。

对CELF的运算时间进行分析：

option 2:运行CELF A.size()次 （Initialization time为计算所有$\sigma(\{u\})$的时间）

```
	CELF done. total time = 15.234, simulation time = 15.234(Initialization time = 14.81200)

​	CELF done. total time = 16.266, simulation time = 16.266(Initialization time = 16.07900)

​	CELF done. total time = 55.531, simulation time = 55.531(Initialization time = 55.35900)

​	CELF done. total time = 49.063, simulation time = 49.063(Initialization time = 48.87500)

​	CELF done. total time = 24.359, simulation time = 24.359(Initialization time = 24.26500)
```

option 2:运行CELF advanced.

```
CELF advanced done. total time = 77.360, simulation time = 77.344(Initialization time = 67.50000)
```

**结论：CELF所有的时间基本在MC simulation上。并且在Initialization部分做的simulation占大多数。**





对MC_simulation，猜想可能是vector存取限制了性能，考虑换成数组。

实验结果：

1. -O3加速下，vector time = 4.516, array time = 4.578
2. 不加速，vector time = 30.437, array time = 22.250

**结论：没有O3的话加速明显，可惜在O3下换不换并没有什么区别**



对MC_simulation，猜想可能是meet次数太多导致性能较慢。

选定一个size=24，邻居总数=8700的集合进行simulation（仅1次）。

实验结果：在deadline=15的TC/IC模型下，平均尝试meet 547000次；在IC模型下，平均尝试activate 123000次.

分析：在IC模型，一个节点最多被其邻居节点各尝试激活1次，时间上限为$O(|V|+|E|)$; 在IC-M模型，每轮迭代一个节点都会被每个active neighbour尝试连接1次，时间上限为$O((|V|+|E|)*deadline)$。当deadline较大时，IC-M模型的simulation比IC模型慢。

残念ながら，即使现在的MC simulation已经做到了理论上的最少操作数，每次完整的simulation(r=10000)还是需要进行很多次meet。

**结论：MC_simulation没法再优化了**



**优化一。**观察到每进行一组数据的实验，都要运行CELF很多次，而每次CELF都要花大量时间计算$\sigma(\{u\})$。

想到可以把$\sigma(\{u\})$用数组存下来。

优化前（先做一次option 2 CELF，再做一次advanced CELF）：

```
CELF method done. total time = 161
CELF advanced done. total time = 75
```

优化后：

```
CELF method done. total time = 68
CELF advanced done. total time = 8
```

**结论：效果拔群！**



**优化二。**之前生成active participant的方式是每当size改变，就完全随机生成。这样存下的$\sigma(\{u\})$没有被充分利用。

假设之前已经算了size = 10的active participant，当size增加到20时，只需要往原来的集合中添加10个点皆可，这样原来节点的邻居$u$的$\sigma(\{u\})$就不需要再算一遍。



优化一，二的缺点：因为CELF的initialization部分被优化掉了，难以比较两个算法在这一部分的时间差异。



*7/21 19:22 已经将优化一，二的version放到服务器上跑了。*

Advanced CELF部分：ap size和k较大时，算法在迭代过程中花费了大量时间。

原因：MC simulation的开销随着seed size增长很快。

```
Seed = {175,232,448,738,963,1886,2091,2240,2442,2520,2785,2799,2802,3560,3857,3907,4373,4495,4838,9918}
size = 20
Simulation result=67.497 time=68.274 meet time=4923270677

Seed ={175,232,445,448,644,666,738,990,1225,1886,2091,2240,2442,2520,2785,2799,2802,3197,3406,3560,3857,3907,4373,4495,4651,4838,4983,4996,7805,8156,8269,9918}  
size = 32
Simulation result=114.294 time=129.433 meet time=9199051402
```

当ap size = 10, k = 10时，advanced CELF需要计算最大size=100的set的MC simulation.

相比option2的celf，只需要计算最大size=10的MC simulation.



7/22

所有点的single spread已经存储到本地文件data/edges_mg.txt.

新增了命令行参数识别工具，目前支持参数如下。

-?，--help，输出命令行帮助并退出程序。

-v，--verbose，在标准流输出verbose message。

-l，--local，采用本地的single spread。

7/23

增加了一些注释，调整了一些代码结构。

新增编译参数：-r，--rounds，MC simulation每次迭代的次数。

2022/7/25 1：00

新增IMM.h，内有所有IMM相关方法。

目前已经实现的方法：

1. RI_Gen: 生成一个RI set或者FI set
2. NodeSelection: 用贪心算法，返回对$\R$的maximum coverage problem的解：点集$S_k$（未debug）
3. FI_simulation: 用FI sketch来近似估计$\sigma(S)$

还剩下的方法：

1. Sampling: 确定需要生成的$\R$的大小，并生成足够多的RI set

2. IMM Algorithm

   ```pseudocode
   Procedure: IMM(G, k, e, l)
   begin procedure
   	l = l * (1 + log(2) / log(n))
   	R = Sampling(G, k, e, l)
   	S* = NodeSelection(R, k)
   	return S*
   end procedure
   ```

性能：

RI_Gen的复杂度为$O(n^*logn^*+m^*)$，其中$n^*$和$m^*$为在集合中被激活的节点数和相关边数。

基于RI_Gen的FI simulation与MC simulation对比，结果准确度相同，并且FI simulation的runtime一般为MC simulation的$\frac{1}{10}$.

**截至目前的所有实验结果在根目录的results文件夹。**

**当ap size=10,k=10时，CELF advanced已经非常的慢了(result_10_10.out, 总时长约24小时)。**

*意见：如果要继续做ap size更大的部分，可以考虑用更快的FI sim代替MC sim?*

7/26

为了防止爆long long情况出现，重新规范了数据类型。

节点index与节点的个数 - node - int64_t

k - int32 - int32_t

取值$\in [0,1]$的整数 - int8_t

其他值 - int64 - int64_t

实现了*IMM Algorithm*，并实现了option 2的IMM版本。调用方式参考main.cpp与test.cpp。

option2思路：在IMM算法的NodeSelection部分， 将"Identify the vertex $v$"的取值范围设为$\N_{a_i}$(active participant的neighbour set)，其余不变。



7/27

实现了IMM-Advance部分。

基本测试场景：数据集blog-catalog，ap size = 500, k = 10

第一次测试：

```pseudocode
seed quality of PAGERANK = 1370.5            time = 0.13      size = 641.333
seed quality of PAGERANK_ADVANCED = 3551.44  time = 0.162     size = 2380.33
seed quality of IMM_NORMAL = 1432.29         time = 495.453   size = 662.667
seed quality of IMM_ADVANCED = 3680.63       time = 0.265667  size = 2489
```

IMM_NORMAL对每个ap node生成约500000个RI set，一共生成500 * 500000个，时间较慢；

IMM_ADVANCED一次性生成约100000个RI set，比IMM_NORMAL中IMM单次生成的数量要少。



注意到类似CELF的优化思路，可以把RI set存下来，让IMM_NORMAL时间大量减少。

第二次测试：

```pseudocode
seed quality of PAGERANK = 1300.6            time = 0.146     size = 603.667
seed quality of PAGERANK_ADVANCED = 3624.72  time = 0.166667  size = 2459.67
seed quality of IMM_NORMAL = 1350.13         time = 4.625     size = 614.667
seed quality of IMM_ADVANCED = 3758.8        time = 0.0416667 size = 2548.67
```

**结论：IMM算法可以拓展到ap size较大的情况，且对比PageRank，IMM生成解集的质量略高。**

进一步的实验有待操作。

