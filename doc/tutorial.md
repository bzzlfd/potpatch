
这个 tutorial 将以[汪老师文章][wang] Si bulk 中单个 Al 原子替换一个 Si 原子的 Si:Al 作为具体例子演示计算流程并提示注意事项. 

本文假设读者掌握基本的 PWmat 弛豫/自洽/非自洽计算方法. 

## 1.为 potential patch 准备晶体结构和势场文件
需要提前准备的文件只有 bulk 晶体结构文件. 
如果你在进行一个需要开 SOC 的计算, 可以从本步骤开始使用全相对论赝势, 等 patch 之后 Escan 计算再开 SOC. 浅能级对结构影响不大, 这一步开 SOC 与否不会有很大影响, 当然你在这里就打开 SOC 也没有什么问题. 

### 1.1.2 LDA bulk SCF
准备 bulk SCF 计算的输入文件, 进行自洽计算.

`etot.input` considerations:
1. `N123`: 计算输出势场文件 `OUT.VR` 实际上存储了一个实空间3维离散网格, 网格数由 晶格常数, `Ecut` 和并行参数 共同影响. 为了后续 patch 过程中网格匹配, 我非常建议在这里显式地设置 `N123`. 你可以看看[这篇笔记](./Ecut_n123_AL.md)进一步了解晶格常数, `Ecut2` 和 `N123` 的关系. 
2. `XCFUNCTIONAL`: 建议采用 `XCFUNCTIONAL = LDA` , PBE泛函会让势场出现小锯齿, 这不利于potentail patch. 相应的, 赝势也建议用 LDA 赝势.
3. `CONVERGENCE`: 非常建议设置 `CONVERGENCE=DIFFICULT`, 这样生成的势场在[康老师文章][kang] Fig3 检查中符合得更好. 
4. `OUT.VATOM`: 建议打开 `OUT.VATOM`, 这是一会儿 (调整赝势环节) 要用到的米奇妙妙工具. 你可以把它当作效仿绘制[康老师文章][kang] Fig.3 的数据点来源. 
5. `OUT.WG`: 如果你认为波函数文件 `OUT.WG` 没用的话可以设置取消输出它以节省硬盘空间. 

这个例子需要计算VBM, 考虑是否进行非自洽计算得到bulk的VBM, 但是Si的VBM在 $\Gamma$ 点, 自洽计算已经取到了, 可以用 `Gap_Read` 读取


### 1.2.1 LDA supercell RELAX
从1.1中使用的 atom.config 制作supercell的晶体结构文件, 进行弛豫计算. 

相较于普通的制作超胞, 这里为了 potential patching 效果需要额外固定弛豫时在 patching 边界位置的原子, `potpatch` 程序提供了一个小程序 `potpatch mksupcl` 帮助完成这件事. 

关于超胞尺寸, 经验上 (对这个例子) 使用 4×4×4 超胞就足够满足需求了. 它的原子位置基本不受带电杂质的影响, 它的势场相较于 bulk 只剩下带电杂质产生的库伦势. 请做自己的检查. 

将 `[0, 0, 0]` (最近) 位置处的原子换成杂质原子. `potpatch` 程序假设杂质原子位于 `[0, 0, 0]` 进行 patch. 如果想要掺杂的原子位置不在 `[0, 0, 0]`, 可以使用 `potpatch shift` 子程序对原子位置进行平移. 

`etot.input` considerations:
1. `N123`, `Ecut`, `Ecut2`: 此阶段对这些参数没有要求. 
2. `IN.PSP`: 相较于bulk计算, 这里引入了新的杂质的赝势
3. `MP_N123`: 超胞可以在倒空间少采样几个点, 这个例子中用单Gamma点已经足够了


### 1.2.2 LDA supercell SCF
把 RELAX 的结果 `final.config` 当作本任务的 `atom.config`

supercell的自洽和bulk的自洽在 `etot.input` 没有很大的区别
1. `N123`: 因为这个例子使用  4×4×4 超胞, 相应的 `N123 = 128 128 128`
2. `IN.PSP`, `NUM_ELECTRON`, `MP_N123`: 和上一步弛豫一样
3. `NUM_ELECTRON`: 让体系是 close shell 会让结果更准 (这通常会导致体系带电). 从杂质体系可能的带电状态中选择一个 closed shell 的价电子数作为 `NUM_ELECTRON` 的数值. 在这个例子中, 原本没有缺陷的 4×4×4 supercell 有2048个价电子, 替换一个 Si 原子为 Al 原子后中性体系有 2047 个价电子, acceptor 获得一个电子变成 close shell 后有 2048 个电子, 所以应该设置 `NUM_ELECTRON = 2048`. 



## 2. HSE计算
我们将来会用 HSE 参考对 LDA 的能级进行修正, 不如我们现在就连 HSE 的结果也算出来吧. 
HSE 相对 LDA 的计算, 只需要把 `XCFUNCTIONAL` 改成 `XCFUNCTIONAL = HSE`; 和 `OUT.WG` 一样, 如果你经过考虑认为 `OUT.HSEWR` 不需要, 可以设置 `OUT.HSEWR = F` 以节省硬盘空间. 


## 3.patch 和 escan 计算
接下来会使用 folded spectrum method 计算能级, 它可以很快地 (linear scaling) 计算出已知原子位置和势场的很大体系的少量能量本征值和本征波函数, Escan 中有该方法的一个实现, 它可以通过在 PWmat 当中设置 NONSCF 任务参数启用. 
`potpatch` 程序做的工作是在之前第一步 LDA 计算出的较小体系的原子位置和势场文件 patch 成一个很大体系的原子位置和势场文件以投入Escan计算. 

### 3.1 LDA Escan bulk
这个计算设置基于普通的 nonSCF job, 同样需要准备 `atom.config` 和 `IN.VR`, 这是 bulk 计算, 直接从 LDA bulk SCF 当中复制过来就行了. 

`etot.input` 基于普通的 `NONSCF` job :
1. `N123`: 为了确保 PWmat 不出错, 这里的 `N123` 需要使用与SCF时同样的设置. 
2. `IN.VR` 和 `IN.KPT`: 需要额外指定. (这个 tutorial 出于通用性的考虑使用 `IN.KPT` 文件, 在 `etot.input` 中正确设置 `MP_N123` 以替换 `IN.KPT` 也可以.)
3. `OUT.WG`: 视情况可以关掉
4. `IN.NONSCF`: 根据 Pwmat 文档要求, 需要设置成 `T`
5. `PRECISION`: 根据 Pwmat 文档要求, 需要设置成 `DOUBLE`
6. `NUM_BAND`: 根据 PWmat 文档, 这个参数决定 Escan 算几个能级, 根据自己的需要设置. 这个例子中我们只关心 1 个能级, 设成 `1` 就好了

`IN.NONSCF`
1. `NONSCF_METH = 2`: 执行 FSM 计算
2. `FSM_EREF`: 可以在第一节 LDA SCF 计算之后用 `Gap_Read` 参考 VBM CBM, 设置一个使 FSM 中 $(\hat{H}-\varepsilon_\text{ref})^2$ 的第一个能级就是 VBM 的值. 
   1. 想让 FSM 算出的第一个能级是 VBM 是因为 Si:Al 的浅能级位于 VBM 附近. 
   2. 一个建议是让 $\varepsilon_\text{ref}$ 尽量远离 VBM , 可以让计算收敛得更快. 

`IN.KPT`:
1. Si 的 VBM 在 $\Gamma$ 点, 我们在 `IN.KPT` 文件中设置我们要进行 $\Gamma$ 点的非自洽计算. 

bulk 计算的 `atom.config` 和 `IN.VR` 文件来自 LDA bulk SCF 计算, 它并不需要使用 `potpatch` 程序.

按道理说这一步的结果应该和之前 bulk 自洽计算算出来的结果一致. 

### 3.2 LDA Escan suuuupercell
`etot.input`, `IN.NONSCF`, `IN.KPT` 设置与 Escan bulk 计算差别不大, 只需要新增杂质元素的赝势; 根据 patch 后尺寸设置 `N123` 就可以了. 

Escan 的 `atom.config` 和 `IN.VR` 文件通过 `potpatch` 程序生成. 让我们编写 `potpatch` 程序的控制文件吧. 本例子提供了 `potpatch.input` 模板文件, 好了我们编写完了. 

控制 `potpatch` 主程序的参数在 `[potpatch]` 这个 section 之下. 
1. `[potpatch.bulk]` 和 `[potpatch.supercell]` 控制输入 `potpatch` 的信息, 包含供 patch 的 bulk 和 supercell 晶体结构 和 势场文件位置; 
2. supercell带电情况 `charge` , 它被定义为计算 supercell 时 PWmat 中 `(setting NUM_ELECTRON) - (default NUM_ELECTRON)` 得到的数值. 在本例子中 `charge = 2048-2047 = 1`
3. 以及该体系的介电常数 `epsilon` , 在本例子中, Si 的相对介电常数是12.34. 
4. `[potpatch.supercell]` 下的 `frozen_range` 和 `size` 不是必需的, 不设置它们 `potpatch` 也会正常工作, 设置这些参数帮助用户确认在 `potpatch` 执行过程中 supercell 与是否用户预期的一致. 
   1. `potpatch` 程序打印它目前了解到的信息. 你也可以使用 `potpatch --only-inspect` 让 `potpatch` 仅打印信息, 在开始执行 potentail patch 前退出运行. 
5. `[potpatch.output]` 控制 `potpatch` 程序输出的信息, 包括 `atom.config` 和 `VR` 文件名, 和想要 patch 生成的 suuuupercell 的尺寸. 
6. `[potpatch.check.diff_vatom]` 是可选参数. 当它被设置, 会额外计算 supercell 体系上各原子位势. 
   1. 它同样可以当作效仿绘制[康老师文章][kang] Fig.3 的数据点来源. 
   2. 它和 PWmat OUT.VATOM 文件所使用的计算方法不同. 

[note_potpatch](./note_potpatch.md) 处有该输入文件参数的完整介绍. 

在 `potpatch` 运行过程中, 会有一些输出, 它们可以帮助你判断 potential patching 是否出错了. "standard deviation of diffs at the boundary" 通常只有几 meV, 如果它太大请小心, 考虑使用更大的 supercell 以及检查之前的计算是否含有错误. 

当 suuuupercell 尺寸是 8a 时计算得到的 binding energy 的结果是 81.25 meV, 很接近[汪老师文章][wang]中的 80.1 meV, 我用两张 1080ti 花了 6 分钟完成计算. 大成功.


## 4.修正赝势
为了在 LDA 下算准, 接下来做的是调整杂质附近原子的赝势, 让 $E_{\mathrm{im}}^{\mathrm{LDA}+\mathrm{C}}\left(\Omega_{512}\right)$ 变化, 使得下式成立. 

$$
E_{\mathrm{im}}^{\mathrm{HSE}}\left(\Omega_{512}\right)-E_{\mathrm{VBM}}^{\mathrm{HSE}}=E_{\mathrm{im}}^{\mathrm{LDA}+\mathrm{C}}\left(\Omega_{512}\right)-E_{\mathrm{VBM}}^{\mathrm{LDA}}
$$

这一步或许可以被理解为在原本 LDA 哈密顿量基础上 *在杂质原子附近原子的某个角动量通道* 加了微扰. 

对于 [UPF 格式](https://pseudopotentials.quantum-espresso.org/home/unified-pseudopotential-format) 赝势文件, 比较经济的改赝势方案是更改 `<PP_DIJ>` 这个按角动量分块的分块对角方阵. 
尤其是对于 Kleinman-Bylander 形式的赝势, [这么搞一定不会弄出来 ghost state][gonze]. 
想了解赝势文件中 `<PP_LOCAL>`/`<PP_NONLOCAL>`/`<PP_BETA>`/`<PP_DIJ>` 标签在指代什么变量? 不妨去看看 [Vanderbilt][vanderbilt] 这篇简洁优雅的文章. Quantum Espresso 将超软赝势和模守恒赝势一起处理, 文章中模守恒赝势中的 $B_{ij,l}$ 矩阵兼容地在 UPF 文件中被 $D_{ij,l}$ 矩阵表示. 

需要注意的是, 实际上这四个能量本征值并不是在同一基准上的, 我们需要把离杂质最远的原子位势当做四个体系的基准, 对上面公式的各个能量进行对齐. 

这个例子中有效质量预测结果很好所以没有进行修正, 请做自己的检查. 


[wang]: https://doi.org/10.1063/1.3153981
[kang]: https://doi.org/10.1103/PhysRevApplied.18.064001
[vanderbilt]: https://doi.org/10.1103/PhysRevB.41.7892
[gonze]: https://doi.org/10.1103/PhysRevB.44.8503
