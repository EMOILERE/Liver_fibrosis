# 示例配置文件说明

这个目录包含了不同实验场景的示例配置文件，帮助您快速开始不同类型的模拟实验。

## 📁 文件说明

### 基础配置
- `basic_test.xml` - 基础测试配置，快速验证模拟是否正常工作
- `full_experiment.xml` - 完整48小时实验配置，用于正式研究

### 协同效应研究
- `high_synergy.xml` - 高协同效应配置，用于验证强协同作用
- `low_synergy.xml` - 低协同效应配置，用于对比研究
- `no_synergy.xml` - 无协同效应配置，纯加性效应

### 剂量反应研究
- `low_dose.xml` - 低剂量TGF-β1和外泌体
- `medium_dose.xml` - 中等剂量（标准配置）
- `high_dose.xml` - 高剂量配置

### 时间进程研究
- `short_term.xml` - 短期实验（12小时）
- `medium_term.xml` - 中期实验（24小时）
- `long_term.xml` - 长期实验（72小时）

### 性能优化
- `fast_test.xml` - 快速测试配置，较少细胞和较短时间
- `high_performance.xml` - 高性能配置，优化并行设置

## 🚀 使用方法

### 方法1: 直接复制
```bash
# 复制想要的配置为主配置文件
copy example_configs\high_synergy.xml config\simulation_parameters.xml

# 然后正常编译运行
make clean && make
./liver_fibrosis_igem.exe
```

### 方法2: 使用快速配置脚本
```bash
# 使用Python脚本加载特定配置
python config/quick_config.py --config example_configs/high_synergy.xml --output config/simulation_parameters.xml

# 或者直接指定配置文件
./liver_fibrosis_igem.exe example_configs/high_synergy.xml
```

### 方法3: 参数对比
```bash
# 运行多个配置并比较结果
for config in example_configs/*.xml; do
    echo "Running $config"
    ./liver_fibrosis_igem.exe $config
    # 保存结果到不同目录
done
```

## ⚙️ 配置文件重点参数对比

| 配置文件 | 协同因子 | TGF-β1 (ng/mL) | 模拟时间 | 细胞数 | 用途 |
|----------|----------|----------------|----------|--------|------|
| basic_test | 2.5 | 10.0 | 12h | 50 | 快速测试 |
| high_synergy | 3.5 | 10.0 | 48h | 100 | 强协同验证 |
| low_synergy | 1.5 | 10.0 | 48h | 100 | 弱协同对比 |
| high_dose | 2.5 | 20.0 | 48h | 100 | 高剂量研究 |
| fast_test | 2.5 | 10.0 | 6h | 30 | 极速验证 |

## 🔧 自定义配置

### 创建新配置
1. 复制基础配置文件
2. 修改关键参数
3. 保存为新名称
4. 测试运行

### 常用修改参数
```xml
<!-- 协同效应强度 -->
<max_synergy_factor>2.5</max_synergy_factor>

<!-- TGF-β1浓度 -->
<concentration>10.0</concentration>

<!-- 模拟时间 -->
<simulation_time>2880</simulation_time>

<!-- 实验组选择 -->
<current_experiment>dual_miRNA_1_1</current_experiment>
```

## 📊 结果验证

不同配置应该产生以下预期结果：

### high_synergy.xml
- 双miRNA组应显示明显超过单miRNA组效果之和
- 协同指数 < 0.7
- 纤维化标志物显著降低

### low_synergy.xml  
- 双miRNA组效果接近单miRNA组效果之和
- 协同指数 ≈ 0.9-1.1
- 轻微的协同增效

### no_synergy.xml
- 双miRNA组效果等于单miRNA组效果之和
- 协同指数 ≈ 1.0
- 纯加性效应

## ⚠️ 注意事项

1. **参数合理性**: 确保参数在生物学合理范围内
2. **计算资源**: 高细胞数和长时间模拟需要更多计算资源
3. **结果对比**: 使用相同的随机种子以便结果对比
4. **配置备份**: 修改前备份原配置文件

## 🔍 故障排除

### 常见问题
- **模拟过慢**: 使用fast_test.xml验证
- **协同效应不明显**: 尝试high_synergy.xml
- **内存不足**: 减少细胞数量或域大小
- **编译错误**: 检查配置文件XML格式

### 验证流程
1. 先运行basic_test.xml确保基础功能正常
2. 再运行目标配置进行正式实验
3. 对比不同配置的结果差异
4. 根据需要调整参数

---

有问题请参考主README.md或检查参数校准指南！
