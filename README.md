#### Libcell

一种高维网格表示结构。能表示高维网格(线段，曲面，体结构...)，单形和非单形的网格，流形和非流形的网格...

libcell的遍历接口仿照openmesh,一方面openmesh的遍历接口望文生义，一方面是为了降低学习成本。

#### Dependencies

* Eigen（不再维护相关代码）(Deprecated)
* Algebras（optinal）
* cstructures

基于Eigen库的算法头文件依赖Eigen库，Eigen放到thirdpart文件夹.

我已经不再维护EIgen库的算法，目前我重写了基于Algebras库的算法。

Algebras的张量库可以代替EIgen.

基于Algebras库的算法头文件依赖Algebras库，Algebras放到thirdpart文件夹。

#### install
新建build文件夹，并进入.

```bash
cmake ..
make 
make install
```
#### 算法
利用Algebras库写了任意维的凸包算法，任意维的delauny剖分算法,任意维的最优传输算法，这些算法基于Algebra库的反对称张量，基本做到"无视背景空间维度"下的任意维子流形算法。
#### 可视化

推荐使用Viewer库，代码紧凑，功能齐全,仅仅依赖opengl，跨平台，只需添加必要数据即可显示.

libcell_tools_view.h文件添加了libcell支持Viewer的j小函数，方便你的使用.

#### 遗留的工作

* 各种算法的完善

#### 支持的网格结构

* cell网格结构。市面上没有针对高维网格的数据结构，cell格式文件为仿照off自定义文件格式，非常简单（意味这极低的学习成本）

当网格为单形时，理论上可以推导各个结构层的信息，所以当simplex=1时，cell文件格式只需要点的信息和单形的信息。



当网格为非单形时，理论上你只能获得cell文件格式提供的信息，无法推导上一层或下一层结构信息(除非提供)。所以cell文件格式需要你提供点，"半面"（n-1维流形），由“半面”组合的“胞腔”的信息(n维流形)。

* tools_formats.h提供了off,obj,mesh文件格式的转换

#### 几种网格库的对比

市面上是找不到任意维网格库，所以libcell不但提供这方面的库，也降低的各种概念成本。比如任意维流形的边界算子...,

这是其他库无法比拟的,
* Openmesh 这是开放一种二维流形网格库，支持非单形网格表示(但是概念和使用方法不一，增加学习成本)，相同的团队还开发了三维流形网格库openvolume.缺点是各种概念和工具的转换。优点是速度非常快，遍历接口望文生义，且统一。