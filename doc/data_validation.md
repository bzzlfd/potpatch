# 数据一致性检查

`VR`、`Atom` 和 `MaterialSystemInfo` 可以从空对象逐步构造。赋值只修改目标属性，不会同步其他对象，也不会自动检查。`MaterialSystemInfo.validate()` 读取当前数据并返回检查报告，同时把这次报告保存在 `last_validation`。

构造函数接收晶格时会复制它，避免多个对象意外共用一块可修改的数组。之后若直接把同一个 `Lattice` 实例赋给多个子对象，仍遵循普通 Python 引用语义。

```python
info.vr.lattice = new_lattice
report = info.validate()
for check in report.checks:
    print(check.name, check.status, check.detail, check.sources)
report.raise_for_errors()
```

每项检查的状态是 `pass`、`fail` 或 `insufficient`。后者表示缺少比较所需的数据。空对象不会因为没有冲突而被判为全部通过。`last_validation` 是上一次检查时的快照；直接修改属性或 NumPy 数组后，应再次调用 `validate()`。

网格和原子数组需要使用相应数值类型的 NumPy 数组；晶格必须是有限且非退化的三维晶格。检查只读取这些数据，不替调用方做类型转换。

计算入口可以通过 `required_paths` 指明必需字段。缺少这些字段会产生 `fail`；仅调用 `validate()` 时，未提供的可选字段允许为空。

```python
report = info.validate(required_paths=(
    "atomconfig.lattice", "atomconfig.natoms", "atomconfig.positions",
    "vr.lattice", "vr.mesh",
))
report.raise_for_errors()
```

主 `potpatch` 流程在检查输入、开始计算和写出结果前重新检查。单独调用 `VR.write()` 或 `Atom.write()` 也会先检查写出所需字段。bulk 与 supercell 的尺寸关系仍由 `inspect_ingredient()` 检查，因为这是两个体系之间、且与具体计算有关的约束。

## 0.2.0 不兼容变更

移除了 `VR` 和当时名为 `AtomConfig`（现为 `Atom`）的类构造函数中的 `lattice_check_trigger` 参数，以及三个对象通过 `__setattr__` 实现的晶格自动同步。设置 `info.lattice` 不再覆盖子对象的晶格；设置 `info.vr.lattice` 也不会改动 `info.atomconfig.lattice`。请在修改完成后主动检查，并在计算或写出前处理报告中的失败项。
