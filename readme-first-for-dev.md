# PWmat Escan 只能算 50000 以下原子
想要复现汪老师的 potential patching 的文章需要扩胞至 20a, 需要计算 64000 个原子. 

催更! 

# 更新赝势部分
1. 我希望这部分也作为 `potpatch` 的一个子命令, 比如 `potpatch modifypsp`, 可以通过一个文件设置赝势(或者通过参数传递), (至于修正赝势的具体形式, 或许可以把它焊死?)
2. 我能想到有两个潜在代码需要完成, UPF格式赝势对象, 寻找近邻原子
3. 具体的设计还需要思考

拜托您了

# 我意识到还没有进一步思考/做的事
1. 这代码被设计可以用来在旧Escan上用, 只要读写文件的时候把默认参数设置成 `fmt="Escan"`, 但是那段检查代码应该严格一点 (不过只要没人知道这个参数, 用默认值, 也没啥问题)
2. 在一些数据类型微妙的地方加一些 type hint 
3. 加上 --verbose 参数, 比如
   1. 我想查看MaterialSystemInfo里面lattice被set的过程
   2. 我想检查时间大户都花了多长时间
   3. 计算过程更多输出信息
4. `potpatch.input` 设置 output 时, 如果不存在这个目录会报错, 算 Escan 的时候要先创建文件夹再运行命令太怪了. 
   1. 或许应该 catch 这个错误, 然后问用户是不是想要创建文件夹
5. `potpatch.input` 文件中, epsilon 参数放在 supercell 内感觉有点不自然, 因为 bulk 当然也是这个参数
   1. 文档, input_file_parse, correction, __main__ 代码需要改
6. tutorial 文档中给出一些参考, 改变什么会使能级有多大变化, potentail patch在算很浅的能级, 很多失误都会对结果产生可观的影响, 很难分辨是自己出错了还是 potentail patch 只能精确到这里. 
7. 文档也写成英文的

绝赞咕咕中, 能用就行.jpg

# 我不确定的事
1. `objects.VR.write_vr()`, `nndoes` 那一段会不会有什么问题, 我用没集成在PWmat上的Escan的时候在这里踩过坑, 或许现在版本如果不小心地设置VR文件中的nnodes还是可能存在问题. 
2. VR 是二进制 Fortran 文件, 每段数据前后有一个 int32 的标注, 记录这段数据有多少字节, 我不确定是不是所有机器上都是int32
3. 大小端也让我担心, 我不懂这具体是怎样的, 但我隐约担心它会导致Fortran二进制文件处理时出错
4. eigen 文件里面 `nref_tot_8`  (从`TDM/TDM.f90`偷的) 是啥. 
5. 我不清楚 $N_i$ 的公式应该是 $\sqrt{2Ecut2}\times AL_i/\pi$ 还是 $\sqrt{2Ecut2} /(ALI_i \times \pi)$
