
这个 tutorial 将以[汪老师文章][wang] Si bulk 中单个 Al 原子替换一个 Si 原子的 $Si_{Al}$ 作为具体例子演示计算流程并提示注意事项. 

## 1.为 potential patch 准备晶体结构和势场文件
需要提前准备的文件只有 bulk 晶体结果文件. 
如果你在进行一个需要开 SOC 的计算, 可以从本步骤开始使用 SOC 的赝势, 等 patch 之后 Escan 计算再开 SOC. shallow level 对结构影响不大, 这一步开 SOC 与否不会有很大影响, 当然你在这里就打开 SOC 也没有什么问题. 

### 1.1.2 LDA bulk scf
准备 bulk SCF 计算的输入文件, 进行自洽计算.

`etot.input` considerations:
1. `N123`: 计算输出势场文件 `OUT.VR` 实际上存储了一个实空间3维离散网格, 网格数由 晶格常数, `Ecut` 和并行参数 共同影响. 为了后续 patch 过程中网格匹配, 我非常建议在这里显式地设置 `N123`. 你可以看看[这篇笔记](./Ecut_n123_AL.md)进一步了解晶格常数, `Ecut2` 和 `N123` 的关系. 
2. `XCFUNCTIONAL`: 建议采用 `XCFUNCTIONAL = LDA` , PBE泛函会让势场出现很多小锯齿, 这不利于potentail patch. 相应的, 赝势也建议用 LDA 赝势.
3. `CONVERGENCE`: 非常建议设置 `CONVERGENCE=DIFFICULT`, 这样生成的势场在[康老师文章][kang] Fig3 检查中符合得更好. 
4. `OUT.VATOM`: 建议打开 `OUT.VATOM`, 这是一会儿 (调整赝势环节) 要用到的米奇妙妙工具. 你可以把它当作效仿绘制[康老师文章][kang] Fig3 的数据点来源. 
5. `OUT.WG`: 如果你认为波函数文件 `OUT.WG` 没用的话可以设置取消输出它以节省硬盘空间. 

这个例子需要计算VBM, 考虑是否进行非自洽计算得到bulk的VBM, 但是Si的VBM在 $\Gamma$ 点, 自洽计算已经取到了, 可以用 `Gap_Read` 读取


### 1.2.1 LDA supercell relax
从1.1中使用的 atom.config 制作supercell的晶体结果文件, 进行弛豫计算. 

相较于普通的制作超胞, 这里为了 potential patching 效果需要额外固定弛豫时在 patching 边界位置的原子, `potpatch` 程序提供了一个小程序 `potpatch mksupcl` 帮助完成这件事. `example/potpatch.input` 中有控制 mksupcl 行为的参数, 使用 `potpatch mksupcl -i INPUT` 告诉 `potpatch mksupcl` 去哪里找 `INPUT`(在这个例子中, 即`example/potpatch.input`) 文件, [note_potpatch](./note_potpatch.md) 处有关于 `potpatch mksupcl` 的解释. 

关于超胞尺寸, 经验上 (对这个 8 原子正方体 lattice 例子) 使用 4×4×4 超胞就足够满足需求了. 它的原子位置基本不受带电杂质的影响, 它的势场相较于 bulk 只剩下带电杂质产生的库伦势. 请做自己的检查. 

将 [0, 0, 0] (最近) 位置处的原子换成杂质原子. `potpatch` 程序假设杂质原子位于 [0, 0, 0] 进行 patch. 

`etot.input` considerations:
1. `N123`, `Ecut`, `Ecut2`: potpatch 程序对这些参数没什么要求, 可以随便设
2. `IN.PSP`: 相较于bulk计算, 这里引入了新的杂质的赝势
3. `NUM_ELECTRON`: 让体系是 close shell 会让结果更准 (这通常会导致体系带电). 从杂质体系可能的带电状态中选择一个 closed shell 的价电子数作为 `NUM_ELECTRON` 的数值. 在这个例子中, 原本没有缺陷的 4×4×4 supercell 有2048个价电子, 替换一个 Si 原子为 Al 原子后中性体系有 2047 个价电子, acceptor 获得一个电子变成 close shell 后有 2048 个电子, 所以应该设置 `NUM_ELECTRON = 2048`. 
4. `MP_N123`: 超胞可以在倒空间少采样几个点, 这个例子中用单Gamma点已经足够了


### 1.2.2 LDA supercell scf
把 relax 的结果 `final.config` 当作本任务的 `atom.config`

supercell的自洽和bulk的自洽在 `etot.input` 没有很大的区别
 1. `N123`: 因为这个例子使用  4×4×4 超胞, 相应的 `N123 = 128 128 128`
 2. `IN.PSP`, `NUM_ELECTRON`, `MP_N123`: 和上一步弛豫一样



## 2. HSE计算
我们将来会用 HSE 参考对 LDA 的能级进行修正, 不如我们现在就连 HSE 的结果也算出来吧. 
HSE 相对 LDA 的计算, 只需要把 `XCFUNCTIONAL` 改成 `XCFUNCTIONAL = HSE`; 和 `OUT.WG` 一样, 如果你经过考虑认为 `OUT.HSEWR` 不需要, 可以设置 `OUT.HSEWR = F` 以节省硬盘空间. 
虽然我们用的 LDA 赝势, 而 HSE 是 PBE 和 HF 的混合, 按道理应该用 PBE 赝势, 但实践上会看到即使是 LDA 赝势效果也很好. 


## 3.patch 和 escan 计算
接下来会使用 folded spectrum method 计算能级, 它可以很快地 (linear scaling) 计算出已知原子位置和势场的很大体系的少量能量本征值和本征波函数, Escan 中有该方法的一个实现, 它可以通过在 PWmat 当中设置 NONSCF 任务参数启用. 
`potpatch` 程序做的工作是在之前第一步 LDA 计算出的较小体系的原子位置和势场文件 patch 成一个很大体系的原子位置和势场文件以投入Escan计算. 

### 3.0 potpatch
`potpatch`程序的实现中, 在 patch 前后, (假设) 杂质原子一直位于 [0,0,0] 坐标, 这是为什么之前把杂质原子设置在最靠近 [0,0,0] 的位置. 
直接执行 `potpatch -i INPUT` 即开始 potentail patching 过程, `INPUT` 是告诉 `potpatch` 程序它需要的信息的文件. 本例子中是 `example/potpatch.input` 文件. 这个文件参数很简单, 应该看一下就知道应该怎么用了. 
[note_potpatch](./note_potpatch.md) 处有关于 `potpatch` 程序和输入文件的文档. 

### 3.1 LDA Escan bulk
这个计算设置基于普通的 nonscf job, 同样需要准备 `atom.config` 和 `IN.VR`, 这是 bulk 计算, 直接从 LDA bulk scf 当中复制过来就行了. 

`etot.input` 基于普通的 `NONSCF` job :
1. `N123`: 为了确保 PWmat 不出错, 这里的 `N123` 需要使用与scf时同样的设置. 
2. `IN.VR` 和 `IN.KPT`: 需要额外指定. (这个 tutorial 出于通用性的考虑使用 `IN.KPT` 文件, 在 `etot.input` 中正确设置 `MP_N123` 以替换 `IN.KPT` 也可以.)
3. `OUT.WG`: 视情况可以关掉
4. `IN.NONSCF`: 根据 Pwmat 文档要求, 需要设置成 `T`
5. `PRECISION`: 根据 Pwmat 文档要求, 需要设置成 `DOUBLE`
6. `NUM_BAND`: 根据 PWmat 文档, 这个参数决定 Escan 算几个能级, 根据自己的需要设置. 这个例子中我们只关心 1 个能级, 设成 `1` 就好了

`IN.NONSCF`
1. `NONSCF_METH = 2`: 执行 FSM 计算
2. `FSM_EREF`: 可以在第一节 LDA scf 计算之后用 `Gap_Read` 参考 VBM CBM, 设置一个使 FSM 中 $(\hat{H}-\varepsilon_{ref})^2$ 的第一个能级就是 VBM 的值. 
   1. 想让 FSM 算出的第一个能级是 VBM 是因为 $Si_{Al}$ 的 shallow level 位于 VBM 附近. 
   2. 一个可选的建议是让 $\varepsilon_{ref}$ 尽量远离 VBM , 可以让计算收敛得更快. 

`IN.KPT`:
1. Si 的 VBM 在 $\Gamma$ 点, 我们在 `IN.KPT` 文件中设置我们要进行 $\Gamma$ 点的非自洽计算. 

bulk 计算的 `atom.config` 和 `IN.VR` 文件来自 LDA bulk scf 计算, 它并不需要使用 `potpatch` 程序.

按道理说这一步的结果应该和之前 bulk 自洽计算算出来的结果一致. 

### 3.2 LDA Escan suuuupercell
`etot.input`, `IN.NONSCF`, `IN.KPT` 设置与 bulk 计算差别不大, 只需要新增杂质元素的赝势; 根据 patch 后尺寸设置 `N123` 就可以了. 

Escan 的 `atom.config` 和 `IN.VR` 文件通过 `potpatch` 程序生成. 让我们编写 `potpatch` 程序的控制文件吧. 本例子提供了 `potpatch.input` 模板文件, 好了我们编写完了(没有).

控制 `potpatch` 主程序的参数在 `[potpatch]` 这个  header 之下. 
1. `[potpatch.bulk]` 和 `[potpatch.supercell]` 控制输入 `potpatch` 的信息, 包含供 patch 的 bulk 和 supercell 晶体结构 和 势场 文件相对于 `potpatch.input` 文件的位置; 
2. supercell带电情况 `charge` , 它被定义为计算supercell时 PWmat 中 `(setting NUM_ELECTRON) - (default NUM_ELECTRON)` 得到的数值. 在本例子中 `charge = 2048-2047 = 1`
3. 以及该体系的介电常数 `epsilon` , 在本例子中, Si 的相对介电常数是12.34. `[potpatch.supercell]` 下的 `frozen_range` 和 `size` 两个参数与 `mksupcl` 中的概念一样, 不需要设置它们 `potpatch` 也会正常工作, 单独设计这两个参数的原因是, supercell 不一定是 `potpatch mksupcl` 生成的, 这些参数帮助用户确认在 `potpatch` 执行过程中 supercell 与用户预期的一致. 
   1. 你也可以使用 `potpatch --only-inspect` 打印更多 `potpatch` 程序目前了解到的信息, 并在程序开始执行 potentail patch 前终止. 
4. `[potpatch.output]` 控制 `potpatch` 程序输出的信息, 包括 `atom.config` 和 `VR` 文件名, 和想要 patch 生成的 suuuupercell 的尺寸. 

[note_potpatch](./note_potpatch.md) 处有 `potpatch` 程序和输入文件的文档. 

在 `potpatch` 运行过程中, 会有一些输出, 它们可以帮助你判断 potential patching 是否出错了. "standard deviation of diffs at the boundary" 通常只有 几meV, 如果它太大请小心, 考虑使用更大的 supercell 以及检查之前的计算是否含有错误. 

当 suuuupercell 尺寸是 8a 时计算得到的 binding energy 的结果是 81.25meV, 很接近[汪老师文章][wang]中的 80.1 meV, 计算在两张1080ti上花了6分钟. 大成功.


## 4.修正赝势
为了在 LDA 下算准, 接下来做的是调整杂质附近原子的赝势, 让 $E_{\mathrm{im}}^{\mathrm{LDA}+\mathrm{C}}\left(\Omega_{512}\right)$ 变化, 使得下式成立. 

$$
E_{\mathrm{im}}^{\mathrm{HSE}}\left(\Omega_{512}\right)-E_{\mathrm{VBM}}^{\mathrm{HSE}}=E_{\mathrm{im}}^{\mathrm{LDA}+\mathrm{C}}\left(\Omega_{512}\right)-E_{\mathrm{VBM}}^{\mathrm{LDA}}
$$

这一步相当于在原本 LDA 哈密顿量基础上 在杂质原子附近原子的某个轨道 加了一个局域在杂质原子位置附近的微扰哈密顿量. 
这里有一个细节, 实际上这四个 eigenvalue 并不是在同一基准上的, 这是为什么我们一开始计算的时候打开了 `OUT.VATOM`, 我们需要从 `OUT.VATOM` 中获取信息, 把离杂质最远的原子处的势场当做四个体系的基准, 修正上面的公式. 

这个例子没有做有效质量的修正是因为在修正之前已经符合得很好了, 请做自己的检查. 



[wang]: https://doi.org/10.1063/1.3153981
[kang]: https://doi.org/10.1103/PhysRevApplied.18.064001
