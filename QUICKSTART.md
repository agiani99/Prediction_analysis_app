# 🚀 Quick Start Guide

## Getting Started in 3 Steps

### Step 1: Install Dependencies (First Time Only)

**Option A: Using pip**
```bash
pip install -r requirements.txt
```

**Option B: Using conda (Recommended for RDKit)**
```bash
conda create -n compound_dashboard python=3.9
conda activate compound_dashboard
conda install -c conda-forge rdkit
pip install streamlit pandas numpy plotly seaborn scikit-learn scipy umap-learn
```

### Step 2: Launch the Dashboard

**Linux/Mac:**
```bash
./launch_dashboard.sh
```

**Windows:**
```cmd
launch_dashboard.bat
```

**Manual Launch:**
```bash
streamlit run compound_prediction_dashboard.py
```

### Step 3: Upload Your Data

1. Click "Browse files" in the sidebar
2. Upload your predictions CSV
3. Upload your protein mapping CSV
4. Start exploring!

---

## 🎯 Common Use Cases

### Use Case 1: Explore Kinase Predictions

1. **Model Selection:**
   - Mode: inhibition
   - Model: xgb
   - Type: cls

2. **Family Filter:**
   - Select: Kinase

3. **Visualize:**
   - Go to "Dimensionality Reduction" tab
   - Choose UMAP
   - Color by a specific kinase
   - Look for clusters of active compounds

### Use Case 2: Find Similar Compounds

1. **Navigate to "Compound Search" tab**

2. **Select a hit compound:**
   - Choose CompoundID from dropdown

3. **Set similarity threshold:**
   - Number of similar compounds: 10-20

4. **Find Similar:**
   - Click "Find Similar Compounds"
   - Review structures and predictions

5. **Export results**

### Use Case 3: Correlation Analysis

1. **Select protein family:**
   - Example: GPCR

2. **Go to "Correlation Analysis" tab**

3. **Examine heatmap:**
   - Look for highly correlated proteins (red/blue)
   - Identify protein clusters

4. **Review statistics:**
   - Check mean correlation
   - Find top correlated pairs

5. **Download correlation matrix**

### Use Case 4: Compare Models

1. **Run analysis with:**
   - RF regression

2. **Export results**

3. **Change to:**
   - XGBoost regression

4. **Compare:**
   - Look at correlation differences
   - Check prediction distributions
   - Identify model agreement/disagreement

---

## 📊 Understanding the Outputs

### Dimensionality Reduction Plots

- **PCA**: Linear projections, use for quick overview
- **t-SNE**: Cluster visualization, use for pattern discovery
- **UMAP**: Balance of local and global structure

**What to look for:**
- ✅ Clear clusters = distinct activity profiles
- ✅ Separation by color = good predictions
- ⚠️ Scattered points = diverse compounds

### Correlation Heatmaps

**Color Scale:**
- 🔴 Red (positive): Proteins with similar activity profiles
- 🔵 Blue (negative): Proteins with opposite activity
- ⚪ White (zero): No correlation

**Interpretation:**
- High correlation (|ρ| > 0.8): Consider redundancy
- Moderate correlation (0.5 < |ρ| < 0.8): Related targets
- Low correlation (|ρ| < 0.3): Independent targets

### Similarity Scores

**Cosine Similarity:**
- 1.0 = Identical prediction profile
- 0.5 = Moderately similar
- 0.0 = Completely different

**Use for:**
- Finding backup compounds
- Identifying structure-activity relationships
- Scaffold hopping

---

## 🎨 Customization Tips

### Adjust Performance

**For large datasets:**
```python
# In dashboard settings
sample_size = 500  # Reduce for faster plots
top_n = 5          # Reduce similar compounds
```

### Filter Data

**Focus on active compounds (classification):**
- Min probability: 0.7
- Filter by protein with known activity

**Focus on potent compounds (regression):**
- Set appropriate value ranges
- Filter by specific proteins

### Export Options

**Available exports:**
- ✅ Correlation matrices (CSV)
- ✅ Similar compounds (CSV)
- ✅ Filtered predictions (CSV)
- ✅ Chemical structures (PNG via screenshots)

---

## 🔧 Troubleshooting

### Dashboard won't start
```bash
# Check Streamlit installation
streamlit --version

# Reinstall if needed
pip install --upgrade streamlit
```

### Chemical structures not showing
```bash
# Install/reinstall RDKit
conda install -c conda-forge rdkit
# OR
pip install rdkit-pypi
```

### Out of memory
```bash
# Reduce sample size
# Filter by specific family
# Process fewer proteins at once
```

### Plots not rendering
```bash
# Clear Streamlit cache
# Restart dashboard
# Check browser console for errors
```

---

## 💡 Pro Tips

1. **Start small**: Begin with one family, then expand
2. **Use filters**: Don't load all 100K compounds at once
3. **Save views**: Download interesting subsets for reanalysis
4. **Compare families**: Look for cross-family patterns
5. **Document findings**: Export correlations and similar compounds

---

## 📚 Further Reading

- **Streamlit Documentation**: https://docs.streamlit.io
- **RDKit Documentation**: https://www.rdkit.org/docs
- **Plotly Documentation**: https://plotly.com/python
- **UMAP Guide**: https://umap-learn.readthedocs.io

---

## 🆘 Getting Help

**Common Questions:**

**Q: Can I analyze multiple families at once?**
A: Yes! Select "All" in the family filter.

**Q: How do I interpret t-SNE plots?**
A: Focus on cluster separation, not absolute positions.

**Q: What's the difference between RF and XGB?**
A: Different algorithms; compare both for robust predictions.

**Q: Can I use this with my own data?**
A: Yes! Just follow the CSV format specification.

**Q: How many compounds can I analyze?**
A: Tested with 100K+, but performance depends on your hardware.

---

**Happy Exploring! 🧬**
