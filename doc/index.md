# tutorial
如果你大概了解了 potentail patch 的原理, 而不知道具体应该如何执行计算, 可以去阅读这篇 [tutorial](./tutorial.md)


# 输入文件参数文档
如果你不确定这个程序的某个用法, [note_potpatch](./note_potpatch.md) 处有详细的文档


# 如果想进一步接触代码
按道理在命令行使用 `potpatch` 程序就足够了, 但是我想到有一些人会出于不同的目的会想要进一步接触 `potpatch` 代码: 一种是想要 `import potpatch` 作为自己 python 脚本的一部分的人, 一种是想阅读开发这个代码的人. 
因为目前程序并没有做足够的努力来保证前一种人在使用过程可以和后一种人分离且不遇到困惑, 所以目前这两种人被归入*开发者*. 

在您开发前 [note_src](./note_src.md) 或许会对您有帮助. (额, 文档写成那个鬼样子可能也没有. ) 
(TODO 在未来 note_src 的大部分内容都应该写进代码注释里, note_src 的角色变成介绍一些程序设计思路, 原则, 应该注意到的问题和开发者引路) 

# reference

potential patching

1. [Wang L W. Density functional calculations of shallow acceptor levels in Si[J]. Journal of Applied Physics, 2009, 105(12).](https://doi.org/10.1063/1.3153981)
2. [Zhang G, Canning A, Grønbech-Jensen N, et al. Shallow impurity level calculations in semiconductors using ab initio methods[J]. Physical review letters, 2013, 110(16): 166404.](https://doi.org/10.1103/PhysRevLett.110.166404)
3. [Kang J, Wang L W. First-Principles Calculations of Shallow Acceptor-Carbon Complexes in Si: A Potential-Patching Method with a Hybrid-Functional Correction[J]. Physical Review Applied, 2022, 18(6): 064001.](https://doi.org/10.1103/PhysRevApplied.18.064001)

folded spectrum method (FSM)

1. [Wang L W, Zunger A. Solving Schrödinger’s equation around a desired energy: Application to silicon quantum dots[J]. The Journal of Chemical Physics, 1994, 100(3): 2394-2397.](https://doi.org/10.1063/1.466486)


# notes
1. [atomic unit (a.u.)](./atomic_unit.md)
2. [Ecut2, N123 和 AL 的关系](./Ecut_n123_AL.md)
