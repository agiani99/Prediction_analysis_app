#!/bin/bash

# Compound Prediction Dashboard Launcher
# This script checks dependencies and launches the Streamlit dashboard

echo "🧬 Compound Prediction Dashboard Launcher"
echo "=========================================="
echo ""

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "❌ Python 3 is not installed. Please install Python 3.8 or higher."
    exit 1
fi

echo "✅ Python found: $(python3 --version)"

# Check if pip is installed
if ! command -v pip &> /dev/null && ! command -v pip3 &> /dev/null; then
    echo "❌ pip is not installed. Please install pip."
    exit 1
fi

echo "✅ pip found"

# Check if streamlit is installed
if ! python3 -c "import streamlit" 2>/dev/null; then
    echo "⚠️  Streamlit not found. Installing dependencies..."
    pip install -r requirements.txt
    
    if [ $? -ne 0 ]; then
        echo "❌ Failed to install dependencies. Please run manually:"
        echo "   pip install -r requirements.txt"
        exit 1
    fi
    echo "✅ Dependencies installed successfully"
else
    echo "✅ Streamlit found"
fi

# Check if RDKit is installed
if ! python3 -c "import rdkit" 2>/dev/null; then
    echo "⚠️  RDKit not found. Attempting to install..."
    pip install rdkit 2>/dev/null || pip install rdkit-pypi 2>/dev/null
    
    if [ $? -ne 0 ]; then
        echo "⚠️  Could not install RDKit via pip."
        echo "    For best results, install via conda:"
        echo "    conda install -c conda-forge rdkit"
        echo ""
        echo "    Continuing anyway - some features may not work..."
    else
        echo "✅ RDKit installed successfully"
    fi
else
    echo "✅ RDKit found"
fi

echo ""
echo "🚀 Launching dashboard..."
echo "   The dashboard will open in your default browser"
echo "   Press Ctrl+C to stop the server"
echo ""

# Launch streamlit
streamlit run compound_prediction_dashboard.py

# Cleanup on exit
echo ""
echo "👋 Dashboard stopped. Goodbye!"
