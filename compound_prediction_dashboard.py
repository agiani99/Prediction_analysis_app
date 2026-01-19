import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap
from scipy.stats import spearmanr
from scipy.spatial.distance import cdist
import re
from io import BytesIO
import base64
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, Descriptors
import warnings
warnings.filterwarnings('ignore')

# --- NEW: make page config resilient (multipage / reruns) ---
try:
    from streamlit.errors import StreamlitAPIException
except Exception:  # pragma: no cover
    StreamlitAPIException = Exception

def _safe_set_page_config():
    try:
        # Streamlit expects `layout`, not `page_layout` (varies by version).
        try:
            st.set_page_config(
                page_title="Compound Prediction Dashboard",
                layout="wide",
                initial_sidebar_state="expanded",
            )
        except TypeError:
            # Fallback for any Streamlit build that uses a different kwarg name.
            st.set_page_config(
                page_title="Compound Prediction Dashboard",
                page_layout="wide",
                initial_sidebar_state="expanded",
            )
    except StreamlitAPIException:
        # Likely already set by another page/script; don't hard-fail.
        pass

_safe_set_page_config()

# Custom CSS for better styling
st.markdown("""
<style>
    .main-header {
        font-size: 3rem;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 2rem;
    }
    .metric-card {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        margin: 0.5rem 0;
    }
    .stTabs [data-baseweb="tab-list"] {
        gap: 2rem;
    }
    .stTabs [data-baseweb="tab"] {
        height: 3rem;
        padding: 0 2rem;
    }
</style>
""", unsafe_allow_html=True)

@st.cache_data
def load_data(predictions_file, protein_mapping_file):
    """Load prediction data and protein family mapping"""
    # Load predictions
    df_pred = pd.read_csv(predictions_file)
    
    # Load protein mapping
    df_proteins = pd.read_csv(protein_mapping_file)
    
    return df_pred, df_proteins

@st.cache_data
def parse_column_names(_columns):
    """Parse column names to extract protein, mode, model, and type"""
    parsed = []

    for col in _columns:
        if col in ['CompoundID', 'SMILES']:
            continue

        # Pattern: {protein}_{inh|act}_{rf|xgb}_{reg|cls}[__proba|__class]
        match = re.match(
            r'(.+?)_(inhibition|activation)_(rf|xgb)_(reg|cls)\.pkl(?:__(class|proba))?',
            col
        )

        if match:
            protein_name = match.group(1)
            mode = match.group(2)
            model = match.group(3)
            pred_type = match.group(4)
            suffix = match.group(5) if match.group(5) else None

            parsed.append({
                'column': col,
                'protein': protein_name,
                'mode': mode,
                'model': model,
                'type': pred_type,
                'suffix': suffix
            })

    return pd.DataFrame(parsed)

@st.cache_data
def create_protein_family_mapping(df_proteins, parsed_cols):
    """Create mapping between protein names and families"""
    mapping = {}
    
    for _, row in df_proteins.iterrows():
        protein_name = row['HGNC_name']
        families = str(row['Family']).split(';')
        
        for _, col_info in parsed_cols.iterrows():
            col_protein = col_info['protein'].replace('_', ' ')
            
            if protein_name.lower() in col_protein.lower() or \
               col_protein.lower() in protein_name.lower():
                for family in families:
                    if protein_name not in mapping:
                        mapping[protein_name] = {
                            'family': family.strip(),
                            'uniprot': row['UniProtID']
                        }
                break
    
    return mapping

def smiles_to_img(smiles, size=(200, 200)):
    """Convert SMILES to base64 encoded image"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        img = Draw.MolToImage(mol, size=size)
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode()
        return f"data:image/png;base64,{img_str}"
    except:
        return None

def get_molecular_properties(smiles):
    """Calculate basic molecular properties"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {}
        
        return {
            'MW': round(Descriptors.MolWt(mol), 2),
            'LogP': round(Descriptors.MolLogP(mol), 2),
            'HBA': Descriptors.NumHAcceptors(mol),
            'HBD': Descriptors.NumHDonors(mol),
            'TPSA': round(Descriptors.TPSA(mol), 2)
        }
    except:
        return {}

def filter_columns_by_selection(parsed_cols, mode, model, pred_type):
    """Filter columns based on user selection"""
    filtered = parsed_cols[
        (parsed_cols['mode'] == mode) &
        (parsed_cols['model'] == model) &
        (parsed_cols['type'] == pred_type)
    ]
    
    # For classification, only use probability columns
    if pred_type == 'cls':
        filtered = filtered[filtered['suffix'] == 'proba']
    else:
        filtered = filtered[filtered['suffix'].isna()]
    
    return filtered['column'].tolist()

def perform_dimensionality_reduction(data, method='PCA', n_components=2, random_state=42):
    """Perform PCA, t-SNE, or UMAP"""
    if method == 'PCA':
        reducer = PCA(n_components=n_components, random_state=random_state)
        transformed = reducer.fit_transform(data)
        explained_var = reducer.explained_variance_ratio_
        return transformed, explained_var
    
    elif method == 't-SNE':
        reducer = TSNE(n_components=n_components, random_state=random_state, 
                       perplexity=min(30, len(data)-1))
        transformed = reducer.fit_transform(data)
        return transformed, None
    
    elif method == 'UMAP':
        reducer = umap.UMAP(n_components=n_components, random_state=random_state)
        transformed = reducer.fit_transform(data)
        return transformed, None
    
    return None, None

def calculate_spearman_correlation(data):
    """Calculate Spearman correlation matrix"""
    corr_matrix = np.zeros((data.shape[1], data.shape[1]))
    p_values = np.zeros((data.shape[1], data.shape[1]))
    
    for i in range(data.shape[1]):
        for j in range(data.shape[1]):
            if i == j:
                corr_matrix[i, j] = 1.0
                p_values[i, j] = 0.0
            else:
                corr, pval = spearmanr(data[:, i], data[:, j], nan_policy='omit')
                corr_matrix[i, j] = corr
                p_values[i, j] = pval
    
    return corr_matrix, p_values

def find_similar_compounds(compound_idx, prediction_vectors, top_n=10):
    """Find compounds with similar prediction vectors"""
    target_vector = prediction_vectors[compound_idx].reshape(1, -1)
    
    # Calculate cosine similarity (1 - cosine distance)
    distances = cdist(target_vector, prediction_vectors, metric='cosine')[0]
    
    # Get indices of most similar compounds (excluding self)
    similar_indices = np.argsort(distances)[1:top_n+1]
    similarities = 1 - distances[similar_indices]
    
    return similar_indices, similarities

def create_interactive_scatter(df_plot, x_col, y_col, color_col, hover_data, title):
    """Create interactive scatter plot with structure tooltips"""
    fig = px.scatter(
        df_plot,
        x=x_col,
        y=y_col,
        color=color_col,
        hover_data=hover_data,
        title=title,
        template='plotly_white',
        height=600
    )
    
    fig.update_traces(
        marker=dict(size=8, line=dict(width=0.5, color='white')),
        selector=dict(mode='markers')
    )
    
    fig.update_layout(
        hovermode='closest',
        font=dict(size=12),
        title_font_size=16
    )
    
    return fig

def create_correlation_heatmap(corr_matrix, labels, title):
    """Create correlation heatmap using plotly"""
    fig = go.Figure(data=go.Heatmap(
        z=corr_matrix,
        x=labels,
        y=labels,
        colorscale='RdBu_r',
        zmid=0,
        zmin=-1,
        zmax=1,
        text=np.round(corr_matrix, 2),
        texttemplate='%{text}',
        textfont={"size": 8},
        colorbar=dict(title="Spearman ρ")
    ))
    
    fig.update_layout(
        title=title,
        xaxis_title="Proteins",
        yaxis_title="Proteins",
        height=700,
        width=800,
        template='plotly_white'
    )
    
    fig.update_xaxes(tickangle=45)
    
    return fig

@st.cache_data
def _read_csv(source, usecols=None, nrows=None):
    """
    source: either a filesystem path (str) or bytes (from UploadedFile.getvalue()).
    """
    if isinstance(source, (bytes, bytearray)):
        buf = BytesIO(source)
        # pandas may read multiple times; always rewind
        buf.seek(0)
        src = buf
    else:
        src = source  # path-like

    # Try faster engine when available; fall back gracefully.
    try:
        return pd.read_csv(src, usecols=usecols, nrows=nrows, engine="pyarrow")
    except Exception:
        if isinstance(src, BytesIO):
            src.seek(0)
        return pd.read_csv(src, usecols=usecols, nrows=nrows)

def main():
    # Header
    st.markdown('<h1 class="main-header">🧬 Compound Prediction Dashboard</h1>',
                unsafe_allow_html=True)

    # Sidebar for file upload and settings
    with st.sidebar:
        st.header("⚙️ Configuration")

        st.subheader("📁 Data Files")

        st.info(
            "Streamlit’s file uploader size limit is controlled by server config "
            "(`server.maxUploadSize` in MB). For multi‑GB CSVs, use 'Local path'."
        )

        pred_source_mode = st.radio(
            "Predictions source",
            options=["Upload", "Local path"],
            horizontal=True,
        )
        map_source_mode = st.radio(
            "Protein mapping source",
            options=["Upload", "Local path"],
            horizontal=True,
        )

        predictions_bytes_or_path = None
        proteinmap_bytes_or_path = None

        if pred_source_mode == "Upload":
            predictions_file = st.file_uploader(
                "Upload Predictions CSV",
                type=["csv"],
                help="Subject to Streamlit server upload limit."
            )
            if predictions_file is not None:
                predictions_bytes_or_path = predictions_file.getvalue()
        else:
            predictions_path = st.text_input(
                "Predictions CSV path (on the machine running Streamlit)",
                placeholder=r"C:\path\to\predictions.csv",
            )
            if predictions_path:
                predictions_bytes_or_path = predictions_path

        if map_source_mode == "Upload":
            protein_mapping_file = st.file_uploader(
                "Upload Protein Mapping CSV",
                type=["csv"],
                help="Subject to Streamlit server upload limit."
            )
            if protein_mapping_file is not None:
                proteinmap_bytes_or_path = protein_mapping_file.getvalue()
        else:
            proteinmap_path = st.text_input(
                "Protein Mapping CSV path (on the machine running Streamlit)",
                placeholder=r"C:\path\to\protein_mapping.csv",
            )
            if proteinmap_path:
                proteinmap_bytes_or_path = proteinmap_path

        if predictions_bytes_or_path and proteinmap_bytes_or_path:
            st.success("✅ Sources configured")
        else:
            st.warning("⚠️ Please provide both data sources to begin")
            st.stop()

    # --- CHANGED: read only headers first (fast, low memory) ---
    with st.spinner("Loading headers..."):
        df_pred_header = _read_csv(predictions_bytes_or_path, nrows=0)
        df_proteins = _read_csv(proteinmap_bytes_or_path)

        parsed_cols = parse_column_names(df_pred_header.columns)
        protein_family_map = create_protein_family_mapping(df_proteins, parsed_cols)

    # Main controls
    st.sidebar.header("🎛️ Model Selection")
    
    col1, col2, col3 = st.sidebar.columns(3)
    
    with col1:
        mode = st.selectbox(
            "Mode",
            options=['inhibition', 'activation'],
            help="Select inhibition or activation predictions"
        )
    
    with col2:
        model = st.selectbox(
            "Model",
            options=['rf', 'xgb'],
            help="Random Forest or XGBoost"
        )
    
    with col3:
        pred_type = st.selectbox(
            "Type",
            options=['reg', 'cls'],
            help="Regression or Classification"
        )
    
    # Filter columns based on selection
    selected_columns = filter_columns_by_selection(parsed_cols, mode, model, pred_type)
    
    if len(selected_columns) == 0:
        st.error("No columns found for the selected combination!")
        st.stop()
    
    # Extract protein names and families for selected columns
    protein_info = []
    for col in selected_columns:
        col_data = parsed_cols[parsed_cols['column'] == col].iloc[0]
        protein_name = col_data['protein'].replace('_', ' ')
        
        # Find family
        family = 'Unknown'
        for prot_key, prot_data in protein_family_map.items():
            if prot_key.lower() in protein_name.lower():
                family = prot_data['family']
                break
        
        protein_info.append({
            'column': col,
            'protein': protein_name,
            'family': family
        })
    
    df_protein_info = pd.DataFrame(protein_info)
    
    # Family selection
    st.sidebar.header("🔬 Protein Family Filter")
    available_families = sorted(df_protein_info['family'].unique())
    
    selected_family = st.sidebar.selectbox(
        "Select Protein Family",
        options=['All'] + available_families,
        help="Filter by protein family for analysis"
    )
    
    # Filter by family
    if selected_family != 'All':
        family_columns = df_protein_info[df_protein_info['family'] == selected_family]['column'].tolist()
    else:
        family_columns = selected_columns
    
    if len(family_columns) == 0:
        st.warning(f"No proteins found in family: {selected_family}")
        st.stop()
    
    # --- CHANGED: now read only the columns we actually need ---
    required_cols = ["CompoundID", "SMILES"] + family_columns
    with st.spinner("Loading selected prediction columns..."):
        df_pred = _read_csv(predictions_bytes_or_path, usecols=required_cols)

    # Prepare prediction matrix
    prediction_matrix = df_pred[family_columns].values

    # Remove rows with NaN values
    valid_mask = ~np.isnan(prediction_matrix).any(axis=1)
    prediction_matrix_clean = prediction_matrix[valid_mask]
    df_pred_clean = df_pred[valid_mask].reset_index(drop=True)
    
    # Display metrics
    st.header("📊 Dataset Overview")
    
    metric_col1, metric_col2, metric_col3, metric_col4 = st.columns(4)
    
    with metric_col1:
        st.metric("Total Compounds", f"{len(df_pred_clean):,}")
    
    with metric_col2:
        st.metric("Selected Proteins", len(family_columns))
    
    with metric_col3:
        st.metric("Protein Family", selected_family)
    
    with metric_col4:
        st.metric("Model", f"{model.upper()} {pred_type.upper()}")
    
    # Main tabs
    tab1, tab2, tab3, tab4 = st.tabs([
        "🎯 Dimensionality Reduction",
        "🔥 Correlation Analysis",
        "🔍 Compound Search",
        "📈 Data Explorer"
    ])
    
    # TAB 1: Dimensionality Reduction
    with tab1:
        st.header("Dimensionality Reduction Visualization")
        
        reduction_col1, reduction_col2 = st.columns([1, 3])
        
        with reduction_col1:
            st.subheader("Settings")
            
            reduction_method = st.selectbox(
                "Method",
                options=['PCA', 't-SNE', 'UMAP'],
                help="Select dimensionality reduction method"
            )
            
            if reduction_method == 'PCA':
                n_components = st.slider("Components", 2, 3, 2)
            else:
                n_components = 2
            
            color_by = st.selectbox(
                "Color by",
                options=['None'] + family_columns,
                help="Select a protein to color points by prediction"
            )
            
            # --- CHANGED: make slider safe for small datasets ---
            n_rows = len(df_pred_clean)
            if n_rows < 2:
                st.error("Not enough compounds after filtering (need at least 2).")
                st.stop()

            min_sample = 2 if n_rows < 100 else 100
            max_sample = max(min_sample, min(10000, n_rows))
            default_sample = min(1000, n_rows)

            sample_size = st.slider(
                "Sample Size",
                min_value=min_sample,
                max_value=max_sample,
                value=default_sample,
                step=100 if max_sample >= 100 else 1,
                help="Number of compounds to plot (for performance)"
            )
            
            random_state = st.number_input("Random Seed", value=42, step=1)
            
            if st.button("🚀 Generate Plot", type="primary"):
                st.session_state['generate_reduction'] = True
        
        with reduction_col2:
            if 'generate_reduction' in st.session_state and st.session_state['generate_reduction']:
                with st.spinner(f"Performing {reduction_method}..."):
                    # Sample data
                    if sample_size < len(prediction_matrix_clean):
                        sample_indices = np.random.choice(
                            len(prediction_matrix_clean),
                            size=sample_size,
                            replace=False
                        )
                        sample_matrix = prediction_matrix_clean[sample_indices]
                        sample_df = df_pred_clean.iloc[sample_indices].reset_index(drop=True)
                    else:
                        sample_matrix = prediction_matrix_clean
                        sample_df = df_pred_clean
                    
                    # Perform reduction
                    transformed, explained_var = perform_dimensionality_reduction(
                        sample_matrix,
                        method=reduction_method,
                        n_components=n_components,
                        random_state=random_state
                    )
                    
                    # Prepare plot data
                    plot_df = sample_df[['CompoundID', 'SMILES']].copy()
                    plot_df['Component 1'] = transformed[:, 0]
                    plot_df['Component 2'] = transformed[:, 1]
                    
                    if n_components == 3:
                        plot_df['Component 3'] = transformed[:, 2]
                    
                    # Add molecular properties for hover
                    mol_props = []
                    for smiles in plot_df['SMILES']:
                        props = get_molecular_properties(smiles)
                        mol_props.append(props)
                    
                    for prop in ['MW', 'LogP', 'HBA', 'HBD', 'TPSA']:
                        plot_df[prop] = [mp.get(prop, np.nan) for mp in mol_props]
                    
                    # Color by prediction if selected
                    if color_by != 'None':
                        plot_df['Prediction'] = sample_df[color_by].values
                        color_col = 'Prediction'
                    else:
                        plot_df['Cluster'] = 1
                        color_col = 'Cluster'
                    
                    # Create plot
                    hover_data = {
                        'CompoundID': True,
                        'MW': ':.2f',
                        'LogP': ':.2f',
                        'TPSA': ':.2f'
                    }
                    
                    if color_by != 'None':
                        hover_data['Prediction'] = ':.3f'
                    
                    if n_components == 2:
                        fig = create_interactive_scatter(
                            plot_df,
                            'Component 1',
                            'Component 2',
                            color_col,
                            hover_data,
                            f"{reduction_method} Projection - {selected_family}"
                        )
                    else:
                        # 3D plot
                        fig = px.scatter_3d(
                            plot_df,
                            x='Component 1',
                            y='Component 2',
                            z='Component 3',
                            color=color_col,
                            hover_data=hover_data,
                            title=f"{reduction_method} 3D Projection - {selected_family}",
                            template='plotly_white',
                            height=700
                        )
                    
                    st.plotly_chart(fig, use_container_width=True)
                    
                    # Show explained variance for PCA
                    if reduction_method == 'PCA' and explained_var is not None:
                        st.info(
                            f"📊 Explained Variance: "
                            f"PC1={explained_var[0]*100:.1f}%, "
                            f"PC2={explained_var[1]*100:.1f}%"
                            + (f", PC3={explained_var[2]*100:.1f}%" if n_components == 3 else "")
                        )
                    
                    # Add structure viewer
                    st.subheader("🔬 Structure Viewer")
                    selected_point = st.selectbox(
                        "Select Compound",
                        options=plot_df['CompoundID'].tolist(),
                        key='structure_viewer'
                    )
                    
                    if selected_point:
                        smiles = plot_df[plot_df['CompoundID'] == selected_point]['SMILES'].values[0]
                        img_data = smiles_to_img(smiles, size=(400, 400))
                        
                        if img_data:
                            st.image(img_data, caption=f"Compound: {selected_point}")
                            st.code(smiles, language=None)
    
    # TAB 2: Correlation Analysis
    with tab2:
        st.header("Spearman Correlation Analysis")
        
        st.info(
            f"📊 Analyzing correlations between {len(family_columns)} proteins "
            f"in the {selected_family} family"
        )
        
        with st.spinner("Calculating Spearman correlations..."):
            corr_matrix, p_values = calculate_spearman_correlation(prediction_matrix_clean)
            
            # Prepare labels
            protein_labels = [
                df_protein_info[df_protein_info['column'] == col]['protein'].values[0]
                for col in family_columns
            ]
            
            # Truncate labels if too long
            protein_labels_short = [
                label[:30] + '...' if len(label) > 30 else label
                for label in protein_labels
            ]
        
        # Display heatmap
        fig_corr = create_correlation_heatmap(
            corr_matrix,
            protein_labels_short,
            f"Protein-Protein Correlation ({selected_family})"
        )
        
        st.plotly_chart(fig_corr, use_container_width=True)
        
        # Correlation statistics
        st.subheader("📈 Correlation Statistics")
        
        # Get upper triangle values (excluding diagonal)
        upper_tri_indices = np.triu_indices_from(corr_matrix, k=1)
        upper_tri_values = corr_matrix[upper_tri_indices]
        
        stat_col1, stat_col2, stat_col3, stat_col4 = st.columns(4)
        
        with stat_col1:
            st.metric("Mean Correlation", f"{np.mean(upper_tri_values):.3f}")
        
        with stat_col2:
            st.metric("Median Correlation", f"{np.median(upper_tri_values):.3f}")
        
        with stat_col3:
            st.metric("Max Correlation", f"{np.max(upper_tri_values):.3f}")
        
        with stat_col4:
            strong_corr = np.sum(np.abs(upper_tri_values) > 0.7)
            st.metric("Strong Correlations (|ρ|>0.7)", strong_corr)
        
        # Top correlated pairs
        st.subheader("🔝 Top Correlated Protein Pairs")
        
        # Create dataframe of correlations
        corr_pairs = []
        for i in range(len(protein_labels)):
            for j in range(i+1, len(protein_labels)):
                corr_pairs.append({
                    'Protein 1': protein_labels[i],
                    'Protein 2': protein_labels[j],
                    'Correlation': corr_matrix[i, j],
                    'Abs Correlation': abs(corr_matrix[i, j]),
                    'P-value': p_values[i, j]
                })
        
        df_corr_pairs = pd.DataFrame(corr_pairs)
        df_corr_pairs_sorted = df_corr_pairs.sort_values('Abs Correlation', ascending=False)
        
        st.dataframe(
            df_corr_pairs_sorted.head(20),
            use_container_width=True,
            hide_index=True
        )
        
        # Download correlation matrix
        csv = df_corr_pairs_sorted.to_csv(index=False)
        st.download_button(
            label="📥 Download Full Correlation Table",
            data=csv,
            file_name=f"correlation_matrix_{selected_family}_{mode}_{model}_{pred_type}.csv",
            mime="text/csv"
        )
    
    # TAB 3: Compound Search
    with tab3:
        st.header("🔍 Similar Compound Finder")
        
        st.write(
            "Find compounds with similar prediction profiles based on "
            "cosine similarity of prediction vectors"
        )
        
        search_col1, search_col2 = st.columns([1, 2])
        
        with search_col1:
            st.subheader("Search Settings")
            
            # Search by compound ID
            compound_search = st.selectbox(
                "Select Compound",
                options=df_pred_clean['CompoundID'].tolist(),
                help="Choose a compound to find similar ones"
            )
            
            top_n = st.slider(
                "Number of Similar Compounds",
                min_value=5,
                max_value=50,
                value=10,
                step=5
            )
            
            if st.button("🔎 Find Similar Compounds", type="primary"):
                st.session_state['search_compound'] = compound_search
                st.session_state['top_n'] = top_n
        
        with search_col2:
            if 'search_compound' in st.session_state:
                compound_id = st.session_state['search_compound']
                top_n = st.session_state['top_n']
                
                # Find compound index
                compound_idx = df_pred_clean[
                    df_pred_clean['CompoundID'] == compound_id
                ].index[0]
                
                # Find similar compounds
                with st.spinner("Finding similar compounds..."):
                    similar_indices, similarities = find_similar_compounds(
                        compound_idx,
                        prediction_matrix_clean,
                        top_n=top_n
                    )
                
                # Display query compound
                st.subheader("Query Compound")
                query_smiles = df_pred_clean.iloc[compound_idx]['SMILES']
                query_img = smiles_to_img(query_smiles, size=(300, 300))
                
                if query_img:
                    st.image(query_img, caption=f"Query: {compound_id}")
                    st.code(query_smiles, language=None)
                
                # Display similar compounds
                st.subheader(f"Top {top_n} Similar Compounds")
                
                similar_df = df_pred_clean.iloc[similar_indices][
                    ['CompoundID', 'SMILES']
                ].copy()
                similar_df['Similarity'] = similarities

                # --- CHANGED: preserve source row index; avoid iloc/iterrows mismatch ---
                similar_df = similar_df.reset_index().rename(columns={'index': 'SourceIndex'})

                # Display as grid
                for _, row in similar_df.iterrows():
                    with st.expander(
                        f"🧪 {row['CompoundID']} - Similarity: {row['Similarity']:.3f}",
                        expanded=False
                    ):
                        col_a, col_b = st.columns([1, 2])

                        with col_a:
                            img = smiles_to_img(row['SMILES'], size=(200, 200))
                            if img:
                                st.image(img)

                        with col_b:
                            st.code(row['SMILES'], language=None)

                            # Show prediction values
                            st.write("**Predictions:**")
                            pred_values = df_pred_clean.loc[row['SourceIndex'], family_columns]
                            pred_df = pd.DataFrame({
                                'Protein': [
                                    df_protein_info[df_protein_info['column'] == col]['protein'].values[0]
                                    for col in family_columns
                                ],
                                'Value': pred_values.values
                            })
                            st.dataframe(pred_df, hide_index=True, height=200)
                
                # Download similar compounds
                csv_similar = similar_df.to_csv(index=False)
                st.download_button(
                    label="📥 Download Similar Compounds",
                    data=csv_similar,
                    file_name=f"similar_to_{compound_id}.csv",
                    mime="text/csv"
                )
    
    # TAB 4: Data Explorer
    with tab4:
        st.header("📊 Data Explorer")
        
        st.subheader("Prediction Data Table")
        
        # Display options
        show_cols = st.multiselect(
            "Select columns to display",
            options=['CompoundID', 'SMILES'] + family_columns,
            default=['CompoundID'] + family_columns[:5]
        )
        
        # Filter options
        st.subheader("Filters")
        
        filter_col1, filter_col2 = st.columns(2)
        
        with filter_col1:
            if pred_type == 'reg':
                min_pred = st.number_input(
                    "Minimum Prediction Value",
                    value=float(prediction_matrix_clean.min()),
                    step=0.1
                )
                max_pred = st.number_input(
                    "Maximum Prediction Value",
                    value=float(prediction_matrix_clean.max()),
                    step=0.1
                )
            else:
                min_pred = st.slider(
                    "Minimum Probability",
                    min_value=0.0,
                    max_value=1.0,
                    value=0.0,
                    step=0.05
                )
                max_pred = st.slider(
                    "Maximum Probability",
                    min_value=0.0,
                    max_value=1.0,
                    value=1.0,
                    step=0.05
                )
        
        with filter_col2:
            filter_protein = st.selectbox(
                "Filter by protein",
                options=['All'] + family_columns
            )
        
        # Apply filters
        display_df = df_pred_clean[show_cols].copy()
        
        if filter_protein != 'All':
            mask = (
                (df_pred_clean[filter_protein] >= min_pred) &
                (df_pred_clean[filter_protein] <= max_pred)
            )
            display_df = display_df[mask]
        
        st.write(f"Showing {len(display_df)} compounds")
        st.dataframe(display_df, use_container_width=True, height=400)
        
        # Summary statistics
        st.subheader("📈 Summary Statistics")
        
        stats_df = df_pred_clean[family_columns].describe()
        st.dataframe(stats_df, use_container_width=True)
        
        # Download data
        csv_data = display_df.to_csv(index=False)
        st.download_button(
            label="📥 Download Filtered Data",
            data=csv_data,
            file_name=f"predictions_{selected_family}_{mode}_{model}_{pred_type}.csv",
            mime="text/csv"
        )

if __name__ == "__main__":
    main()
