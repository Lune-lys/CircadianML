# ⚡ 快速測試指南

選擇最適合您的測試方式：

## 🚀 方法1: 一鍵自動測試（推薦）

```bash
bash run_all_tests.sh
```

這將運行所有測試並生成完整報告，包括：
- ✅ 環境檢查
- ✅ 9個單元測試
- ✅ Demo演示
- ✅ 2個示例程序
- ✅ 輸出文件驗證

**預期結果**: 看到 "🎉 恭喜！所有測試通過！"

---

## 🎯 方法2: 最簡單測試（30秒）

只想快速驗證項目可運行？

```bash
cd /home/user/CircadianML
export PYTHONPATH=src
python -c "from circadian_analysis import demo_analysis; demo_analysis()"
```

**預期結果**:
```
🧬 Circadian Transcriptomics Analysis Demo
...
✅ Demo completed successfully!
```

---

## 🧪 方法3: 逐步測試

### 步驟1: 環境測試
```bash
python --version  # 應該 ≥ 3.8
pip install -r requirements.txt
```

### 步驟2: 單元測試
```bash
export PYTHONPATH=src
python tests/test_rhythm_detection.py
python tests/test_phase_prediction.py
```

### 步驟3: 功能測試
```bash
python examples/quick_start.py
python examples/full_pipeline.py
```

### 步驟4: 查看結果
```bash
ls -R results/
```

---

## 📊 測試結果查看

### 查看生成的圖片
```bash
# 列出所有生成的圖片
find results/ -name "*.png"

# 在瀏覽器中查看（如果支持）
open results/demo/rhythmic_genes.png  # macOS
xdg-open results/demo/rhythmic_genes.png  # Linux
```

### 查看數據文件
```bash
# 查看JTK結果
head -20 results/quick_start/jtk_results.csv

# 統計節律基因數
awk -F',' '$6 < 0.05 {count++} END {print count}' results/quick_start/jtk_results.csv
```

---

## 🎨 測試自己的數據

創建測試腳本：

```python
# test_my_data.py
import sys
sys.path.insert(0, 'src')

from circadian_analysis.rhythm_detection.jtk_cycle import detect_rhythms
from circadian_analysis.utils.helpers import generate_synthetic_data

# 使用合成數據測試
data = generate_synthetic_data(n_genes=100, n_timepoints=24)

# 或讀取您的數據
# import pandas as pd
# data = pd.read_csv('your_data.csv', index_col=0)

# 運行分析
results = detect_rhythms(data)
rhythmic = results[results['BH.Q'] < 0.05]

print(f"檢測到 {len(rhythmic)} 個節律基因")
```

運行：
```bash
python test_my_data.py
```

---

## ✅ 測試成功標準

如果看到以下結果，說明測試成功：

1. **單元測試**: 9個測試全部顯示 ✓
2. **示例程序**: 運行完成且無錯誤
3. **輸出文件**: `results/` 下有多個PNG和CSV文件
4. **數據驗證**: CSV文件包含正確的列（GeneID, BH.Q, LAG等）

---

## 🐛 常見問題

**Q: ModuleNotFoundError**
```bash
export PYTHONPATH=src  # 設置Python路徑
```

**Q: 缺少依賴包**
```bash
pip install -r requirements.txt
```

**Q: 測試運行太慢**
A: 這是正常的，節律檢測需要計算時間。可以減少基因數量測試。

**Q: 圖片不顯示**
A: 圖片已保存為文件，直接查看 `results/` 目錄。

---

## 📚 更多信息

- 完整測試指南: `TEST_GUIDE.md`
- 項目文檔: `README.md`
- 示例代碼: `examples/`

---

**需要幫助？** 檢查錯誤訊息或查看日誌文件 `/tmp/*.log`
