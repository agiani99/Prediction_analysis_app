@echo off
REM Compound Prediction Dashboard Launcher for Windows
REM This script checks dependencies and launches the Streamlit dashboard

echo.
echo 🧬 Compound Prediction Dashboard Launcher
echo ==========================================
echo.

REM Check if Python is installed
python --version >nul 2>&1
if errorlevel 1 (
    echo ❌ Python is not installed. Please install Python 3.8 or higher.
    echo    Download from: https://www.python.org/downloads/
    pause
    exit /b 1
)

echo ✅ Python found
python --version

REM Check if pip is installed
pip --version >nul 2>&1
if errorlevel 1 (
    echo ❌ pip is not installed. Please install pip.
    pause
    exit /b 1
)

echo ✅ pip found

REM Check if streamlit is installed
python -c "import streamlit" >nul 2>&1
if errorlevel 1 (
    echo ⚠️  Streamlit not found. Installing dependencies...
    pip install -r requirements.txt
    
    if errorlevel 1 (
        echo ❌ Failed to install dependencies. Please run manually:
        echo    pip install -r requirements.txt
        pause
        exit /b 1
    )
    echo ✅ Dependencies installed successfully
) else (
    echo ✅ Streamlit found
)

REM Check if RDKit is installed
python -c "import rdkit" >nul 2>&1
if errorlevel 1 (
    echo ⚠️  RDKit not found. Attempting to install...
    pip install rdkit-pypi
    
    if errorlevel 1 (
        echo ⚠️  Could not install RDKit via pip.
        echo     For best results, install via conda:
        echo     conda install -c conda-forge rdkit
        echo.
        echo     Continuing anyway - some features may not work...
    ) else (
        echo ✅ RDKit installed successfully
    )
) else (
    echo ✅ RDKit found
)

echo.
echo 🚀 Launching dashboard...
echo    The dashboard will open in your default browser
echo    Press Ctrl+C to stop the server
echo.

REM Launch streamlit
streamlit run compound_prediction_dashboard.py

REM Cleanup on exit
echo.
echo 👋 Dashboard stopped. Goodbye!
pause
