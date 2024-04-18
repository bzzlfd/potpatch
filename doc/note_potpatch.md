## `potpatch`
### `potpatch -i INPUT`
主体是一个 `potpatch` 程序, 它是一个假设杂质原子位于 [0, 0, 0] 位置的 potential patching 方法实现. 
它通过 `-i, --input INPUT` 参数指定控制文件, 当没指定时, 默认会在当前工作目录下寻找 `potpatch.input` 文件
当读取完输入文件后, 在先对食材进行读取和检查, 然后开始 patch 和 correction 等计算过程, `-I, --only-inspect` 参数可以让计算在开始之前停止并输出一些它的检查结果(或者在输出结果之前报错), 它或许对你检查错误有帮助.

### `potpatch mksupcl -i INPUT`
`potpatch` 还包含一个子命令 `potpatch mksupcl`.
potentail patching 过程中需要制作超胞, 这个超胞在普通超胞的基础上需要在 supercell 和 bulk 边界处固定住原子位置, 这对 potential patching 有好处. `potpatch mksupcl` 会帮忙完成这些工作. 
它的控制文件同样由 `-i, --input INPUT` 指定. 


## input文件
input 文件是一个 toml 格式文件, 如需改写的 input 文件, 请注意 toml 语法. 但通常你不需要了解 toml 语法, 你只要少量更改模板文件中的少量片段足够了. 

它包含了很多 `potpatch` 主程序/子程序 输入/输出 文件位置, 这里的文件位置被设计成相对于 input 文件所处目录的相对路径, 这样做的好处是即使过了很久你也是可以通过检查输入文件了解自己之前做了什么.
你可以把 `mksupcl`, `modifypsp` 等命令的输入文件放在一个文件, 也可以把它们分开放置. 每个(sub)程序只会关心自己 header 下那部分参数
主体 `potpatch` 程序输入文件有一些参数和 `mksupcl` 的输入参数重叠, 但因为我假设 supercell 并不一定来自 `mksupcl` 子程序, `potpatch` 主体程序也不应该依赖它, 所以重叠参数的角色是检查潜在错误而非左右程序运行结果

### `[potpatch]`
这个 header 下的参数控制 `potpatch` 主程序行为. 当执行 `potpatch` 主程序时, 程序只关心这个 header 下的参数
#### `[potpatch.bulk]`
- `atomconfig` (string)
  - 供 patch 的 bulk 的晶体结构文件位置
- `vr` (string)
  - 供 patch 的 bulk 的势场文件位置

#### `[potpatch.supercell]`
- `atomconfig` (string)
  - 供 patch 的 supercell 的晶体结构文件位置
- `vr` (string)
  - 供 patch 的 supercell 的势场文件位置
- `charge` (integer|float)
  - `(setting NUM_ELECTRON) - (default NUM_ELECTRON)`, `NUM_ELECTRON` 是 PWmat 中的参数
- `epsilon` (float)
  - 该体系的介电常数
- `frozen_range` (float) [optional]
  - 含义为距离 supercell 的 patch 边界内 `frozen_range` *angstrom* 内的原子位置与 bulk 一致
  - 这是会被传入 patch 前检查函数的参数, 目的是帮助用户确认他所 patch 的材料全是如他希望的那样
- `size` (array) [optional]
  - supercell 相较于 bulk 的尺寸
  - 这是会被传入 patch 前检查函数的参数, 目的是帮助用户确认他所 patch 的材料全是如他希望的那样

#### `[potpatch.output]`
- `size` (array)
  - 想要 patch 出的 suuuupercell 尺寸
- `atomconfig` (string)
  - patch 出的 suuuupercell 的晶体结构文件存放位置 (目前代码中, 如果指定了不存在的目录会报错而不是创建该目录)
- `vr` (string)
  - patch 出的 suuuupercell 的势场文件存放位置 (目前代码中, 如果指定了不存在的目录会报错而不是创建该目录)

### `[mksupcl]`
这个 header 下的参数控制 `potpatch mksupcl` 子程序行为. 当执行 `potpatch mksupcl` 子程序时, 程序只关心这个 header 下的参数

- `bulk_atomconfig` (string)
  - 供制作 supercell 的 bulk 的晶体结构文件位置
- `supercell_size` (array)
  - 想要制作的 supercell 相对 bulk 的尺寸
- `frozen_range` (float)
  - 在晶体结构文件中, 将距离 supercell 的 patch 边界内 `frozen_range` *angstrom* 内的原子的 `move` 设置成 `0 0 0`
  - 如果 `frozen_range` < 0, 则是不更新原子位置的普通制作超胞
- `output` (string)
  - 所制作超胞的晶体结构文件存放位置 (目前代码中, 如果指定了不存在的目录会报错而不是创建该目录)



