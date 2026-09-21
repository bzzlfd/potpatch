# ECUT、ECUT2、N123 与晶格的关系

## 两种截断能

PWmat 的 `ECUT` 是波函数的平面波截断能；`ECUT2` 是软电荷密度和电势的截断能，单位都是 Ry。电荷密度由波函数及其复共轭的乘积构成，因此其最高波数可达到波函数最高波数的两倍。平面波动能与波数平方成正比，于是充分保留这部分波数时，理想关系是

$$ECUT2_{\mathrm{ideal}} \simeq 4\,ECUT. \tag{1}$$

这不是所有计算必须满足的等式。PWmat 在 `ACCURACY=NORM` 时默认 `ECUT2=2*ECUT`，在 `HIGH` 或 `VERYHIGH` 时默认 `4*ECUT`；也可以手动设置较小的 `ECUT2`，但应检查所需结果的收敛性。对于需要准确原子受力的弛豫计算，PWmat 推荐 `4*ECUT`。参见 [PWmat 手册的 ECUT 和 ECUT2 小节](https://www.pwmat.com/pwmat-resource/Manual.pdf)。

## 估计 N123

`N123` 是由 `ECUT2` 所需波数范围对应的实空间 FFT 网格。先将第 `i` 个晶格矢量的长度从 Å 换成 Bohr；这里使用的 Bohr 半径约为 `0.529177239 Å`：

$$L_i^{\mathrm{Bohr}}=\frac{L_i^{\mathrm{\AA}}}{0.529177239}. \tag{2}$$

在原子单位下，以 Ry 表示的截断能数值等于截止波数的平方。因此，可用第 `i` 个晶格矢量的模长粗估连续网格数：

$$N_i^{\mathrm{est}}\simeq\frac{L_i^{\mathrm{Bohr}}}{2\pi}2\sqrt{ECUT2_{\mathrm{Ry}}}. \tag{3}$$

**固定 `ECUT2` 时，估计网格数与晶格矢量长度成正比。** 晶胞越长，同样的实空间分辨率需要越多网格点。对于非正交晶胞，我们用下面的例子对网格数量进行粗估:

```text
2.02500000  2.02500000  0.00000000
0.00000000  2.02500000  2.02500000
2.02500000  0.00000000  2.02500000
ECUT2 = 80.0 Ry
```

代入式 (2)、(3) 得到

$$L_i^{\mathrm{Bohr}}=\frac{\sqrt{2}\times2.025}{0.529177239}\simeq5.4118,\qquad
N_i^{\mathrm{est}}\simeq\frac{5.4118\sqrt{80}}{\pi}\simeq15.408. \tag{4}$$

连续估计约为 15.4 格，而 PWmat 的计算参数中 `N123 = 16 16 16` , 与估计相符。

PWmat 还需选取适合 FFT 的整数，并满足 `NODE1` 相关的整除要求。这些式子用于预估网格规模，不能代替程序实际输出；即使 `ECUT2` 相同，不同并行设置也可能得到不同的 `N123`。参见 [PWmat 手册的 N123 小节](https://www.pwmat.com/pwmat-resource/Manual.pdf)。
