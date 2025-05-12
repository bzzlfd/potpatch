## `potpatch`
### `potpatch -i INPUT`

`potpatch` 是一个基于杂质原子位于 `[0, 0, 0]` 坐标位置的势场修正 (potential patching) 方法实现工具.

- 通过 `-h` 或 `--help` 查看帮助
- 通过 `-i` 或 `--input` 参数指定控制文件路径
  - 当未指定输入文件时，程序默认在当前工作目录查找 `potpatch.input` 文件
- 调试选项 `-I/--only-inspect`: 
  - 可在计算前中断程序并输出检查结果



## INPUT 文件
输入文件采用 **TOML** 格式, 通常只需修改模板中的关键参数即可, 无需完全掌握 TOML 语法. 

所有相对路径均基于输入文件所在目录解析. 这样做的好处是即使过了很久也可以了解自己之前做了什么. 

所有控制参数在  `[potpatch]` 节下. 

### `[potpatch]`
这个 section 下的参数控制 `potpatch` 主程序行为. 当执行 `potpatch` 主程序时, 程序只关心这个 section 下的参数
#### `[potpatch.bulk]`
- `basedir` (string)
  - 该小节下文件相对目录. 
  - 可以是绝对路径或相对路径, 默认是 `"."`. 当设为相对路径时, 将是基于 INPUT 文件进行解析. 
- `atomconfig` (string)
  - 供 patch 的 bulk 的晶体结构文件位置
- `vr` (string)
  - 供 patch 的 bulk 的势场文件位置

#### `[potpatch.supercell]`
- `basedir` (string)
  - 该小节下文件相对目录. 
  - 可以是绝对路径或相对路径, 默认是 `"."`. 当设为相对路径时, 将是基于 INPUT 文件进行解析. 
- `atomconfig` (string)
  - 供 patch 的 supercell 的晶体结构文件位置
- `vr` (string)
  - 供 patch 的 supercell 的势场文件位置
  
- `frozen_range` (float) [optional]
  - 距离 supercell 的 patch 边界内 `frozen_range` Å 内的原子位置与 bulk 一致
  - 如果它被设置, 程序将在预处理进行额外检查，确保输入体系符合预期
- `size` (array) [optional]
  - supercell 相较于 bulk 的扩胞倍数
  - 如果它被设置, 程序将在预处理进行额外检查，确保输入体系符合预期

#### `[potpatch.correction]`
- `charge` (integer|float)
  - `(setting NUM_ELECTRON) - (default NUM_ELECTRON)`, `NUM_ELECTRON` 是 PWmat 中的控制电子数参数.
- `epsilon` (float|2darray)
  - 该体系的介电常数


#### `[potpatch.output]`
- `basedir` (string)
  - 该小节下文件相对目录. 
  - 可以是绝对路径或相对路径, 默认是 `"."`. 当设为相对路径时, 将是基于 INPUT 文件进行解析. 
- `mkdir` (bool)
  - 是否自动创建目录
  - 默认是 `false`, 此时写入文件的目录不存在将报错.
- `size` (array)
  - 想要 patch 出的 suuuupercell 尺寸
- `atomconfig` (string)
  - patch 出的 suuuupercell 的晶体结构文件存放位置
- `vr` (string)
  - patch 出的 suuuupercell 的势场文件存放位置


#### `[potpatch.check.diff_vatom]`
可选参数, 当它被设置
- `sigma` (float)
  - 高斯加权平均中展宽参数
  - 单位是 Å. 默认是 0.5279 Å. 
- `output` (string)
  - 输出文件名
  - 默认是 `"OUT.diff_vatom"`

