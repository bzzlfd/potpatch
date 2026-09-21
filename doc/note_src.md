本篇文档是阅读和开发源码前的 guide

本代码目前没有 style 规范, 您贡献代码就已经让它向规范进步了. 

- [ ] 用法代码参考


## 一些你可能好奇的问题
在阅读之前你可能在好奇一些问题, 这一节汇总这些问题并建议你去阅读下面的哪一部分. 

1. [为什么我安装之后会有一个 `potpatch` 命令, 它是如何工作的](#入口文件)
2. [实现 potentail patch 的关键代码在哪里?](#potential-patching-的关键)
3. [代码结构是什么样, 我开发时应该注意什么](#代码结构)
4. [`OUT.VR`, `atom.config` 中都有 `Lattice` 变量, 如何在 `MaterialSystemInfo` 中检查一致性](#objects)



## 代码结构
### version
一个存放版本的地方, 为了让每个读入版本的地方都得到一致的版本
使用版本的地方有 `__init__`, `parse.py`, `pyproject.toml`

### 入口文件
`__init__.py`, `__main__.py`, `parse.py` 三个是入口相关文件:
`__init__.py` 定义了一个人通过 `import potpatch` 的方式执行时会发生什么;  
`__main__.py` 定义了一个人通过 `python -m potpatch` 的方式执行时会发生什么; `parse` 是 `__main__.py` 执行过程中处理命令行参数和文件参数的模块, 它们的角色是忠实地向后续程序解读它们读到了什么参数 (默认参数通常返回 `None` ). 
在linux下输入 `type potpatch` 可以找到 `potpatch` 可执行文件所在地, 你可以查看 `pip` 安装的 `potpatch` 程序面纱下面是什么样的. 它被设计成和 `python -m potpatch` 有同样的行为.

### objects
`objects.py` 中定义了很多class, 它们中的大多数都和某个 PWmat 文件一一对应, 当使用其中的 `read` method 时, 更是尽力还原原著
除此之外, 
`MaterialSystemInfo` 是这些数据对象的集合。`VR` 和 `Atom` 可以逐步构造，未提供的属性可以保持为 `None`。重复的数据由显式调用的 `MaterialSystemInfo.validate()` 检查；检查不修改对象。晶格来源仍记录在 `Lattice.fromwhere` 中。

**0.2.0 不兼容变更：**移除了 `lattice_check_trigger` 参数和 `__setattr__` 晶格联动。给 `info.lattice`、`info.vr.lattice` 或 `info.atomconfig.lattice` 赋值，只修改指定位置。需要确认一致性时主动调用 `info.validate()`；计算和写出入口也会执行相应检查。详细用法见 [数据检查](./data_validation.md)。


#### `Lattice`
这是一个有单位的二维数组, 
`AL` 存储三个三维矢量 (each row is a vector, 即期望三个向量分别在内存中是连续的); 
`unit` 用于存储单位: 它的值只允许是 `registered_units` 中的值, 相应的 `au_transformer` 记录若干浮点数, 用以进行单位转换, 可以通过 `self.in_unit()` 输出对应单位下的 `AL` 数组.
`fromwhere` 是记录晶格来源的列表。从文件读取的晶格会记录文件路径；调用 `Lattice.copy(appendwhere=...)` 时可以补充来源。
为了方便, 它定义了equal, 除法和乘法. equal用以检查两个Lattice是否相等. 乘法和除法是与supercell相关的概念, bulk Lattice乘以一个包含三个整数的Sequence返回一个新的supercell Lattice, supercell Lattice除以一个bulk Lattice返回一个包含三个整数的Sequence. 
`Lattice` 可以修改；如果它已经放入一个 `MaterialSystemInfo`，修改完成后应重新调用 `validate()`，因为之前的报告只代表检查当时的数据。

#### `VR`
OUT.VR 文件是一个二进制文件, 它的文件格式从 `convert_rho.f90` 推断出来. 其中AL的单位是angstrom. 在没集成进PWmat中的Escan版本里, AL是原子单位.
VR有一个 `vr_fmt` 参数, 它默认是 `PWmat`, 如果不是这个字符串不是严格的`PWmat`时会在读取/写入过程中使用旧版Escan的文件格式(这个代码设计是不是不合理)

它的关键属性是 `lattice` 和 `mesh`, mesh第三个索引是变化最快的索引; `comment` 是对这个对象的特别注释, 可以帮忙区分它是谁; `n123` 可以从 `mesh` 推断出来, 是衍生量, 用 `@property` 可以确保不出错, 但是这样也不能设置这个参数了; 如果是从文件中读入, 它还会记录 `filename` 和 `nnodes`.

VR被设计成多种用法
在 `__init__` 初始化中, 它可以不传入参数, 传入filename, 传入各个属性, 当同时指定filename和关键属性(如`VR(filename=filename, mesh=mesh)`)时, 会先读取filename当中的nnodes, AL, lattice和mesh, 再(在上面例子中)用传入的mesh覆盖从filename中读取的mesh.
之后还可以通过 `self.attr` 的方式设置属性. 

因为在 write 二进制文件, 需要非常小心各种变量的类型. 这里有一个 `revise_keyattribute_type()` 或许可以被用来做写入前的类型检查 (绝赞咕咕中)

VR 定义了乘法, 它是为超胞准备的, 乘以一个包含三个整数的Sequence返回一个新的supercell VR

#### `Atom`
`Atom` 对应 PWmat 的 `IN.ATOM` 输入项所指向的晶体结构文件（通常名为 `atom.config`）。一个 `Atom` 对象保存整个体系的晶格、所有原子的种类、分数坐标和移动标记，并非单个原子。旧类名 `AtomConfig` 仍可导入，作为 `Atom` 的兼容别名；新代码建议使用 `Atom`。`MaterialSystemInfo.atomconfig` 和控制文件中的 `atomconfig` 键继续沿用原名。

`atom.config` 是一个文本文件，
在没集成进PWmat中的Escan版本里, AL是原子单位, 没有section title, 如(LATTICE, POSITION)
有一个 `atoms_fmt` 参数, 它默认是 `PWmat`, 如果不是这个字符串不是严格的`PWmat`时会在读取/写入过程中使用旧版Escan的文件格式(这个代码设计是不是不合理)

它的关键属性有很多
`natoms`: 位于文件头的属性
`lattice`: 位于 Lattice section 的属性
`ityps`, `positions`, `move`: 位于 Position section 的属性
此外还包含可选属性 `filename` 和 `comment`. 

它的用法设计和 `VR` 一样, 它也定义了被用于制作超胞的乘法运算

#### `VATOM`
OUT.VATOM 是一个文本文件

它的任务只有从文件中读取, 此外没有其他工作. 
具体属性可以检查代码注释, 这个文档中的很多内容也会写到代码注释中, 如果那里才是它们应该待的地方

#### `EIGEN`
OUT.EIGEN 是一个二进制文件

它的任务只有从文件中读取, 此外没有其他工作. 
具体属性可以检查代码注释, 这个文档中的很多内容也会写到代码注释中, 如果那里才是它们应该待的地方

#### `MaterialSystemInfo`
它是对上面所有的对象的集合, 此外还有一些 charge, epsilon 等信息

若显式设置 `info.lattice`，它是独立的参考晶格；否则读取 `info.lattice` 时优先返回 `atomconfig.lattice`，其次返回 `vr.lattice`。这个读取便利性不会复制或同步晶格。调用 `validate()` 会比较所有已存在的来源，报告通过、失败或证据不足的项目。



### Potential Patching 的关键  
`patch.py` 和 `correction.py` 是执行 potentail patching 的关键代码, 相较于用户, 这些代码更接近开发者, 它全部使用Hartree原子单位

$$
V_{\mathrm{im}}(\mathbf{r})=\begin{cases}
V_{\mathrm{im}}^{\mathrm{SC}}(\mathbf{r})-V_C(\mathbf{r})+V_{\text {align }}, & \mathbf{r} \in \Omega_{512} \\
V_{\text {bulk }}(\mathbf{r})+1 / \varepsilon r . & \mathbf{r} \notin \Omega_{512} .
\end{cases}.
$$  
where
$$
V_C(\mathbf{r})=\sum_{(i, j, k) \neq(0,0,0)} \frac{1}{\varepsilon \mid \mathbf{r}-\left(i \mathbf{L}_1+j \mathbf{L}_2+k \mathbf{L}_3 \mid\right)}
$$

`patch.py` 的工作是纯纯的 patch , 使得

$$
V_{\mathrm{im}}(\mathbf{r})=\begin{cases}
V_{\mathrm{im}}^{\mathrm{SC}}(\mathbf{r}), & \mathbf{r} \in \Omega_{512} \\
V_{\text {bulk }}(\mathbf{r}) . & \mathbf{r} \notin \Omega_{512} .
\end{cases}.
$$

`correction.py`  中完成了其余的工作:

`gen_charge_correct` 根据 supercell 材料信息和这个函数其中内置的电荷密度分布函数生成两个函数, `minus_V_periodic` 和 `plus_V_single`. 也用 `gen_` 的方式保证了 `minus_V_periodic` 和 `plus_V_single` 用的是同一个电荷密度分布. 
`minus_V_periodic` 的工作是处理 supercell VR mesh , 它对这个势场减去

$$
V_{periodic}(\mathbf{r})=\sum_{(i, j, k) } \frac{1}{\varepsilon \mid \mathbf{r}-\left(i \mathbf{L}_1+j \mathbf{L}_2+k \mathbf{L}_3 \mid\right)}
$$

`plus_V_single` 的工作是处理 suuuupercell VR mesh, 它对势场加上

$$
V_{single}(r) = 1 / \varepsilon r
$$

二者合起来就是

$$
\begin{cases}
-V_C(\mathbf{r}), & \mathbf{r} \in \Omega_{512} \\
+1 / \varepsilon r . & \mathbf{r} \notin \Omega_{512} .
\end{cases}.
$$

具体实现细节上, 因为球对称的电荷密度在球外的势场和点电荷势场无异, 为了方便做 FFT 解 Poisson 方程, 在 $\Omega_{512}$ 内的电荷电荷密度并不是一个点电荷, 而是具有展宽的, 它的形式大概是

$$
\rho(r)=\begin{cases}
\text{sinc}(r/R_0), & r \lt R_0 \\
0 . & r \ge R_0 .
\end{cases}
$$

球对称电荷密度分布, 它的归一化形式以及单个该电荷密度分布产生的势场的解析形式可以在代码中找到,

`edge_match_correct` 函数做的工作是

$$
\begin{cases}
+V_{\text {align }}, & \mathbf{r} \in \Omega_{512} \\
0 . & \mathbf{r} \notin \Omega_{512} .
\end{cases}.
$$

`__main__.potpatch` 中有做patch和修正的顺序




### 一些基础
`constant.py` 存放了一些常数
`utils.py` 存放的是一些杂七杂八的工具, 包含读写Fortran二进制文件的函数实现
`errors.py` 现在很多错误检查在滥用 assert, 但其中一些检查的语义应该用 raise, 这个文件是放置 raise error 的(鸽了)




