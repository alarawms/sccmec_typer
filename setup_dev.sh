#!/bin/bash
set -e

echo "Setting up development environment..."

# 1. Create Python virtual environment
if [ ! -d ".venv" ]; then
    echo "Creating virtual environment (.venv)..."
    python3 -m venv .venv
fi

# 2. Install dependencies
echo "Installing Python dependencies..."
source .venv/bin/activate
pip install pytest

# 3. Download minimap2
MINIMAP2_VER="2.28"
MINIMAP2_DIR="minimap2-${MINIMAP2_VER}_x64-linux"
MINIMAP2_TAR="${MINIMAP2_DIR}.tar.bz2"
MINIMAP2_URL="https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VER}/${MINIMAP2_TAR}"

if [ ! -f "bin/minimap2" ]; then
    echo "Downloading minimap2 v${MINIMAP2_VER}..."
    wget -q "$MINIMAP2_URL" -O "$MINIMAP2_TAR"
    tar -xjf "$MINIMAP2_TAR"
    
    echo "Installing minimap2 to bin/..."
    mv "$MINIMAP2_DIR/minimap2" bin/
    
    # Cleanup
    rm "$MINIMAP2_TAR"
    rm -rf "$MINIMAP2_DIR"
else
    echo "minimap2 already exists in bin/."
fi

echo "------------------------------------------------"
echo "Setup complete!"
echo "Activate the environment with: source .venv/bin/activate"
echo "Run tests with: ./run_tests.sh (create this script if needed) or just run the python script."
