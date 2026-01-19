# 🧬 Compound Prediction Dashboard

A comprehensive Streamlit dashboard for exploring and analyzing compound prediction data across multiple proteins and families.

## ✨ Features

### 1. **Interactive Model Selection**
- Choose between inhibition/activation predictions
- Select Random Forest (RF) or XGBoost (XGB) models
- Toggle between regression and classification predictions
- Filter by protein families (Kinase, GPCR, Protease, etc.)

### 2. **Dimensionality Reduction Visualizations**
- **PCA (Principal Component Analysis)**: 2D and 3D projections with explained variance
- **t-SNE**: Non-linear dimensionality reduction for pattern discovery
- **UMAP**: Uniform Manifold Approximation and Projection for cluster visualization
- Color by prediction values or explore patterns
- Interactive plots with molecular property tooltips
- Chemical structure viewer for selected compounds

### 3. **Correlation Analysis**
- **Spearman correlation heatmaps** between proteins within families
- Interactive correlation matrices
- Statistical summaries (mean, median, max correlations)
- Top correlated protein pairs identification
- Export correlation tables for further analysis

### 4. **Chemical Structure Integration**
- Hover tooltips showing:
  - Molecular Weight (MW)
  - LogP (lipophilicity)
  - Hydrogen Bond Acceptors/Donors (HBA/HBD)
  - Topological Polar Surface Area (TPSA)
- 2D structure rendering from SMILES
- Structure viewer for compound exploration

### 5. **Similar Compound Finder**
- Find compounds with similar prediction profiles
- Cosine similarity-based search
- Adjustable number of results (5-50)
- Side-by-side structure comparison
- Prediction value comparison tables
- Export similar compound lists

### 6. **Data Explorer**
- Interactive data table with filtering
- Dynamic column selection
- Value range filtering
- Summary statistics
- Export filtered datasets

## 📋 Requirements

- Python 3.8+
- Dependencies listed in `requirements.txt`

## 🚀 Installation

### Step 1: Install Dependencies

```bash
pip install -r requirements.txt
```

### Step 2: Prepare Your Data

Ensure you have two CSV files:

1. **Predictions CSV** (`predictions_EUOS_HTS.csv`):
   - Must contain: `CompoundID`, `SMILES` columns
   - Prediction columns following the pattern:
     ```
     {protein_name}_{inhibition|activation}_{rf|xgb}_{reg|cls}.pkl[__class|__proba]
     ```

2. **Protein Mapping CSV** (`Column_names_pkl_annotated_UID.csv`):
   - Columns: `HGNC_name`, `UniProtID`, `Family`
   - Maps protein names to families

## 🎯 Usage

### Launch the Dashboard

```bash
streamlit run compound_prediction_dashboard.py
```

The dashboard will open in your default web browser at `http://localhost:8501`

### Workflow

1. **Upload Files**
   - Use the sidebar to upload both CSV files
   - Wait for successful loading confirmation

2. **Configure Model**
   - Select Mode: Inhibition or Activation
   - Select Model: RF (Random Forest) or XGB (XGBoost)
   - Select Type: Regression or Classification

3. **Choose Protein Family**
   - Select from available families (Kinase, GPCR, Protease, etc.)
   - Or select "All" to analyze all proteins

4. **Explore Tabs**

   **Tab 1: Dimensionality Reduction**
   - Choose reduction method (PCA/t-SNE/UMAP)
   - Optionally color by specific protein predictions
   - Adjust sample size for performance
   - View 2D or 3D projections
   - Click on compounds to view structures

   **Tab 2: Correlation Analysis**
   - View Spearman correlation heatmap
   - Explore correlation statistics
   - Identify highly correlated protein pairs
   - Download correlation matrices

   **Tab 3: Compound Search**
   - Select a query compound
   - Set number of similar compounds to find
   - View ranked similar compounds with structures
   - Compare prediction profiles
   - Download results

   **Tab 4: Data Explorer**
   - Browse prediction data table
   - Apply filters by value ranges
   - View summary statistics
   - Export filtered datasets

## 📊 Data Format Examples

### Predictions CSV
```csv
CompoundID,SMILES,Protein1_inhibition_rf_reg.pkl,Protein1_inhibition_xgb_reg.pkl,...
COMP001,CCO,0.456,0.478,...
COMP002,CC(C)O,0.234,0.256,...
```

### Classification Predictions
```csv
CompoundID,SMILES,Protein1_activation_rf_cls.pkl__class,Protein1_activation_rf_cls.pkl__proba,...
COMP001,CCO,1,0.875,...
COMP002,CC(C)O,0,0.234,...
```

### Protein Mapping CSV
```csv
HGNC_name,UniProtID,Family
Protein kinase 1,P12345,Kinase
Cytochrome P450 3A4,P08684,Cytochrome_P450
5-hydroxytryptamine receptor 1A,P08908,GPCR
```

## 🎨 Features in Detail

### Dimensionality Reduction Methods

- **PCA**: Best for linear relationships, shows explained variance
- **t-SNE**: Excellent for cluster visualization, non-linear
- **UMAP**: Fast, preserves global structure, good for large datasets

### Correlation Analysis

- **Spearman's ρ**: Non-parametric correlation (robust to outliers)
- **Significance testing**: P-values for each correlation
- **Strong correlation threshold**: |ρ| > 0.7

### Similarity Search

- **Cosine similarity**: Measures angle between prediction vectors
- **Range**: 0 (dissimilar) to 1 (identical)
- **Application**: Find compounds with similar activity profiles

## 🔧 Performance Tips

1. **For large datasets (>10K compounds)**:
   - Use sampling in dimensionality reduction
   - Start with 1000-2000 samples
   - Increase gradually if needed

2. **For many proteins (>50)**:
   - Filter by specific families
   - Focus on proteins of interest
   - Use correlation analysis to identify redundant proteins

3. **For faster loading**:
   - Consider pre-filtering your CSV files
   - Remove unused prediction columns
   - Save filtered versions for repeated use

## 📁 File Structure

```
.
├── compound_prediction_dashboard.py  # Main dashboard application
├── requirements.txt                  # Python dependencies
├── README.md                        # This file
└── data/                            # Your data files (not included)
    ├── predictions_EUOS_HTS.csv
    └── Column_names_pkl_annotated_UID.csv
```

## 🐛 Troubleshooting

### Issue: RDKit installation fails
```bash
# Try conda installation instead
conda install -c conda-forge rdkit
```

### Issue: Out of memory errors
- Reduce sample size in visualizations
- Filter by specific protein families
- Process data in batches

### Issue: Slow performance
- Use smaller sample sizes for plots
- Reduce number of selected proteins
- Consider upgrading to a machine with more RAM

## 📝 Column Naming Convention

The dashboard automatically parses column names following this pattern:

```
{protein_name}_{mode}_{model}_{type}.pkl[__suffix]
```

Where:
- `protein_name`: Protein identifier (can contain underscores)
- `mode`: `inhibition` or `activation`
- `model`: `rf` (Random Forest) or `xgb` (XGBoost)
- `type`: `reg` (Regression) or `cls` (Classification)
- `suffix` (for classification only):
  - `__class`: Predicted class (0 or 1)
  - `__proba`: Probability of class 1

Example:
```
Cytochrome_P450_3A4_activation_rf_cls.pkl__class
Cytochrome_P450_3A4_activation_rf_cls.pkl__proba
Cytochrome_P450_3A4_activation_rf_reg.pkl
```

## 🤝 Contributing

Feel free to extend this dashboard with additional features:
- More dimensionality reduction methods
- Additional molecular descriptors
- Custom clustering algorithms
- Export to different formats
- Integration with external databases

## 📄 License

This dashboard is provided as-is for research and analysis purposes.

## 🙏 Acknowledgments

Built with:
- Streamlit for the web interface
- RDKit for chemistry functionality
- Plotly for interactive visualizations
- scikit-learn for dimensionality reduction
- UMAP for manifold learning

---

**Version**: 1.0.0  
**Last Updated**: January 2025  
**Author**: Claude (Anthropic)
