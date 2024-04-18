本篇文档是阅读和开发源码前的 guideline 

本代码目前没有 style 规范, 您贡献代码就已经让它向规范进步了. 

- [ ] 用法代码参考


## 一些你可能好奇的问题
在阅读之前你可能在好奇一些问题, 这一节汇总这些问题并建议你去阅读下面的哪一部分. 

1. [为什么我安装之后会有一个 `potpatch` 命令, 它是如何工作的](#入口文件)
2. [实现 potentail patch 的关键代码在哪里?](#potential-patching-的关键)
3. [代码结构是什么样, 我开发时应该注意什么](#代码结构)
4. [`OUT.VR`, `atom.config` 中都有 `Lattice` 变量, 如何在 `MaterialSystemInfo` 中保持一致性](#objects)



## 代码结构
### version
一个存放版本的地方, 为了让每个读入版本的地方都得到一致的版本
使用版本的地方有 __init__, parse.py, pyproject.toml

### 入口文件
`__init__.py`, `__main__.py`, `parse.py` 三个是入口相关文件, 
`__init__.py` 定义了一个人通过 `import potpatch` 的方式执行时会发生什么
`__main__.py` 定义了一个人通过 `python -m potpatch` 的方式执行时会发生什么; `parse` 是 `__main__.py` 执行过程中处理命令行参数和文件参数的模块, 它们的角色是忠实地向后续程序反映它们读到了什么(默认参数通常返回 `None` )
在linux下输入 `type potpatch` 可以找到 `potpatch` 可执行文件所在地, 你可以查看 `pip` 安装的 `potpatch` 程序面纱下面是什么样的, 它被设计成和 `python -m potpatch` 有同样的行为

### objects
`objects.py` 中定义了很多class, 它们中的大多数都和某个 PWmat 文件一一对应, 当使用其中的 `read` method 时, 更是尽力还原原著
除此之外, 
`MaterialSystemInfo` 是这些文件对象的统一集合;  
还有一些很多文件共有的部分被单独抽象出来, 同时在各个地方加入了额外的代码以确保一致性, 比如 Lattice 在 atom.config 和 OUT.VR 文件中都存在, 于是
- 在 `MaterialSystemInfo` 中有一个 check_and_modify 函数, 它被植入这个 `MaterialSystemInfo` 下的 `AtomConfig` 和 `VR` 当中
- `MaterialSystemInfo` 当中的 `__setattr__` 被overwrited, 当set atomconfig/vr attr 时会向这个对象当中植入 check_and_modify; 
- `AtomConfig` 和 `VR` 当中的 `__setattr__` 被重写, 任何尝试设置 AtomConfig 和 VR 当中的 `Lattice` 的动作都会触发 check_and_modify(如果存在) 并进行检查与更新
- `MaterialSystemInfo` 当中的 `__setattr__` 被重写, 尝试 set lattice attr 的行为会更新这个对象中 atomconfig/vr 的 lattice
- `Lattice` 当中有个 `fromwhere` 变量用于存储这个变量从何而来, 比如一个 `Lattice` 是从 `VR` 文件中读入, 在触发 trigger 时被设置在 `MaterialSystemInfo` 上时, `VR` 中的 `Lattice` 会记录它是从某个文件读入的, 而 `MaterialSystemInfo` 中的 `Lattice` 会记录它是 `VR` 从某文件读入, 又触发 trigger 设置在 `MaterialSystemInfo` 之上

`MaterialSystemInfo` 会把它的 lattice trigger 遇到的第一个 lattice 信息记录到自己的 lattice 属性中, 之后所有都会来与这个 lattice 对比. 
如果改变 `MaterialSystemInfo` 中的 lattice, 会自动改变它内部 `VR` 和 `AtomConfig` 的lattice
所以在 `MaterialSystemInfo` 不能同时传入 filename 和 关键属性. 
这时候应该自下而上地创建每个对象. 


#### `Lattice`
这是一个有单位的二维数组, 
`AL` 存储三个三维矢量 (each row is a vector, 即期望三个向量分别在内存中是连续的); 
`unit` 用于存储单位: 它的值只允许是 `registered_units` 中的值, 相应的 `au_transformer` 记录若干浮点数, 用以进行单位转换, 可以通过 `self.in_unit()` 输出对应单位下的 `AL` 数组.
`fromwhere` 是一个列表, 每当它被从一个包含`lattice`属性的对象传入另一个包含`lattice`属性的对象时, 它就会在对Lattice进行 *复制* , 并把新对象的注释append到fromwhere中.
为了方便, 它定义了equal, 除法和乘法. equal用以检查两个Lattice是否相等. 乘法和除法是与supercell相关的概念, bulk Lattice乘以一个包含三个整数的Sequence返回一个新的supercell Lattice, supercell Lattice除以一个bulk Lattice返回一个包含三个整数的Sequence. 
请把它当做最小单元使用, 不要手动修改它实例的属性, 如果你想改变其中的某些属性, 请创建一个新的 `Lattice` 实例

#### `VR`
OUT.VR 文件是一个二进制文件, 它的文件格式从 `convert_rho.f90` 推断出来. 其中AL的单位是angstrom. 在没集成进PWmat中的Escan版本里, AL是原子单位.
VR有一个 `vr_fmt` 参数, 它默认是 `PWmat`, 如果不是这个字符串不是严格的`PWmat`时会在读取/写入过程中使用旧版Escan的文件格式(这个代码设计是不是不合理)

它的关键属性是 `lattice` 和 `mesh`, mesh第三个索引是变化最快的索引; `AL_check_trigger` 存储 `MaterialSystemInfo` 为了检查 `Lattice` 一致性而植入的trigger; `comment` 是对这个对象的特别注释, 可以帮忙区分它是谁; `n123` 可以从 `mesh` 推断出来, 是衍生量, 用 `@property` 可以确保不出错, 但是这样也不能设置这个参数了; 如果是从文件中读入, 它还会记录 `filename` 和 `nndoes`.

VR被设计成多种用法
在 `__init__` 初始化中, 它可以不传入参数, 传入filename, 传入各个属性, 当同时指定filename和关键属性(如`VR(filename=filename, mesh=mesh)`)时, 会先读取filename当中的nnodes, AL, lattice和mesh, 再(在上面例子中)用传入的mesh覆盖从filename中读取的mesh.
之后还可以通过 `self.attr` 的方式设置属性. 

因为在 write 二进制文件, 需要非常小心各种变量的类型. 这里有一个 `revise_keyattribute_type()` 或许可以被用来做写入前的类型检查 (绝赞咕咕中)

VR 定义了乘法, 它是为超胞准备的, 乘以一个包含三个整数的Sequence返回一个新的supercell VR

#### `AtomConfig`
atom.config 是一个文本文件, 
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
OUT.EIGEN 是一个文本文件

它的任务只有从文件中读取, 此外没有其他工作. 
具体属性可以检查代码注释, 这个文档中的很多内容也会写到代码注释中, 如果那里才是它们应该待的地方

#### `MaterialSystemInfo`
它是对上面所有的对象的集合, 此外还有一些 charge, epsilon 等信息

它的用法设计沿袭了 `VR`, `AtomConfig`
但是如果在初始化时同时指定了 ①一个它自身含有 lattice 信息的对象 (如`VR`), ②一个该对象的文件名, 它不一定会按照我们所期望的"覆盖"行为工作, 因为可能会遇到lattice冲突. 



### Potential Patching 的关键  
`patch.py` 和 `correction.py` 是执行 potentail patching 的关键代码, 相较于用户, 这些代码更接近开发者, 它全部使用Hartree原子单位

$$
V_{\mathrm{im}}(\mathbf{r})=\left\{\begin{array}{lr}
V_{\mathrm{im}}^{\mathrm{SC}}(\mathbf{r})-V_C(\mathbf{r})+V_{\text {align }}, & \mathbf{r} \in \Omega_{512} \\
V_{\text {bulk }}(\mathbf{r})+1 / \varepsilon r . & \mathbf{r} \notin \Omega_{512} .
\end{array}\right.
$$  
where
$$
V_C(\mathbf{r})=\sum_{(i, j, k) \neq(0,0,0)} \frac{1}{\varepsilon \mid \mathbf{r}-\left(i \mathbf{L}_1+j \mathbf{L}_2+k \mathbf{L}_3 \mid\right)}
$$

`patch.py` 的工作是纯纯的 patch , 使得

$$
V_{\mathrm{im}}(\mathbf{r})=\left\{\begin{array}{lr}
V_{\mathrm{im}}^{\mathrm{SC}}(\mathbf{r}), & \mathbf{r} \in \Omega_{512} \\
V_{\text {bulk }}(\mathbf{r}) . & \mathbf{r} \notin \Omega_{512} .
\end{array}\right.
$$

`correction.py`  中完成了其余的工作:

`gen_charge_correct` 根据 supercell 材料信息和这个函数其中内置的电荷密度分布函数生成两个函数, `minus_V_periodic` 和 `plus_V_single`. 也用 `gen_` 的方式保证了 `minus_V_periodic` 和 `plus_V_single` 用的是同一个电荷密度分布. 
`minus_V_periodic` 的工作是处理 supercell VR mesh , 它对这个势场减去

$$
V_{periodic}(\mathbf{r})=\sum_{(i, j, k) } \frac{1}{\varepsilon \mid \mathbf{r}-\left(i \mathbf{L}_1+j \mathbf{L}_2+k \mathbf{L}_3 \mid\right)}
$$

`plus_V_single` 的工作是处理 suuuupercell VR mesh, 它对势场加上

$$
V_{single} = 1 / \varepsilon r
$$

二者合起来就是

$$
\left\{\begin{array}{lr}
-V_C(\mathbf{r}), & \mathbf{r} \in \Omega_{512} \\
+1 / \varepsilon r . & \mathbf{r} \notin \Omega_{512} .
\end{array}\right.
$$

具体实现细节上, 因为球对称的电荷密度在球外的势场和点电荷势场无异, 为了方便做 FFT 解 Poisson 方程, 在 $\Omega_{512}$ 内的电荷电荷密度并不是一个点电荷, 而是具有展宽的, 它的形式大概是

$$
\rho(r)=\left\{\begin{cases}{lr}
sinc(r/R_0), & r \lt R_0 \\
0 . & r \ge R_0 .
\end{cases}\right.
$$

球对称电荷密度分布, 它的归一化形式以及单个该电荷密度分布产生的势场的解析形式可以在代码中找到,

`edge_match_correct` 函数做的工作是

$$
\left\{\begin{array}{lr}
+V_{\text {align }}, & \mathbf{r} \in \Omega_{512} \\
0 . & \mathbf{r} \notin \Omega_{512} .
\end{array}\right.
$$

`__main__.potpatch` 中有做patch和修正的顺序


### `mksupcl` 脚本所在
`supercell.py` 放置 `mksupcl` 脚本的代码, 但具体的晶体结构操作 (包括扩胞, 修订原子分数坐标, 原子排序等) 的实现位于 `AtomConfig` class 当中. 

### 一些基础
`constant.py` 存放了一些常数
`utils.py` 存放的是一些杂七杂八的工具, 包含读写Fortran二进制文件的函数实现
`errors.py` 现在很多错误检查在滥用 assert, 但其中一些检查的语义应该用 raise, 这个文件是放置 raise error 的(鸽了)






