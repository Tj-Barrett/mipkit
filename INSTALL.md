# MIPkit Installation Guide

## Quick Install

### Standard Installation
```bash
pip install -e .
```

### Development Installation (includes testing tools)
```bash
pip install -e .[dev]
```

### Testing Only
```bash
pip install -e .[test]
```

### Documentation Only
```bash
pip install -e .[docs]
```

## Installation Options Explained

### Optional Dependencies

The package now supports optional dependency groups:

- **`[dev]`** - Full development environment
  - Testing: pytest, pytest-cov, pytest-mock, pytest-xdist, pytest-timeout
  - Code quality: black, flake8, mypy, pylint, isort
  - Coverage: coverage, coverage-badge
  - Documentation: sphinx, sphinx-rtd-theme
  - Development tools: ipython, ipdb

- **`[test]`** - Testing tools only
  - pytest, pytest-cov, pytest-mock

- **`[docs]`** - Documentation generation
  - sphinx, sphinx-rtd-theme, sphinx-autodoc-typehints

## Step-by-Step Installation

### 1. Clone the Repository
```bash
git clone https://github.com/tj-barrett/bmbt-mipkit.git
cd MIPkit
```

### 2. Create Virtual Environment (Recommended)
```bash
# Using venv
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Or using conda
conda create -n mipkit python=3.11
conda activate mipkit
```

### 3. Install MIPkit

For users:
```bash
pip install -e .
```

For developers:
```bash
pip install -e .[dev]
```

### 4. Verify Installation
```bash
# Check MIPkit is installed
MIPkit --help

# Run tests (if dev dependencies installed)
pytest

# Check version
python -c "import MIPkit; print(MIPkit.__version__)"
```

## Dependencies

### Core Dependencies (automatically installed)
- numpy >= 1.20.0
- scipy >= 1.7.0
- scikit-learn >= 1.0.0
- joblib >= 1.0.0
- pandas >= 1.3.0
- matplotlib >= 3.4.0
- rich >= 10.0.0
- rdkit >= 2022.03.1
- pyyaml >= 5.4.0
- seaborn >= 0.11.0
- networkx >= 2.6.0
- pyqt5 >= 5.15.0

### External Tools (must be installed separately)
- **GROMACS** - Molecular dynamics simulations
- **Open Babel** - Chemical file format conversion
- **AutoDock Vina** - Molecular docking
- **ACPYPE** - Topology generation

## Installing External Tools

### GROMACS
```bash
# Ubuntu/Debian
sudo apt-get install gromacs

# macOS
brew install gromacs

# Or from source: http://www.gromacs.org/
```

### Open Babel
```bash
# Ubuntu/Debian
sudo apt-get install openbabel

# macOS
brew install open-babel

# Or via conda
conda install -c conda-forge openbabel
```

### AutoDock Vina
```bash
# Download from: https://github.com/ccsb-scripps/AutoDock-Vina/releases
# Or via conda
conda install -c conda-forge autodock-vina
```

### ACPYPE
```bash
pip install acpype
# Or via conda
conda install -c conda-forge acpype
```

## RDKit Installation

RDKit is a core dependency. If pip installation fails, use conda:

```bash
conda install -c conda-forge rdkit
```

## Development Setup

For contributing to MIPkit:

```bash
# Clone and install in editable mode with dev dependencies
git clone https://github.com/tj-barrett/bmbt-mipkit.git
cd MIPkit
pip install -e .[dev]

# Set up pre-commit hooks (optional)
pip install pre-commit
pre-commit install

# Run tests
pytest

# Check code quality
black --check MIPkit/
flake8 MIPkit/
mypy MIPkit/

# Format code
black MIPkit/
isort MIPkit/
```

## Troubleshooting

### RDKit Installation Issues
If RDKit fails to install via pip, use conda:
```bash
conda create -n mipkit python=3.11
conda activate mipkit
conda install -c conda-forge rdkit
pip install -e .[dev]
```

### Qt/PyQt5 Issues
If PyQt5 causes issues:
```bash
# Try PyQt6 instead
pip uninstall pyqt5
pip install pyqt6
```

### Permission Errors
Use `--user` flag or virtual environment:
```bash
pip install --user -e .
```

### Import Errors
Ensure you're in the virtual environment:
```bash
which python  # Should point to venv/bin/python
python -c "import MIPkit"  # Should not raise ImportError
```

## Updating

### Update MIPkit
```bash
cd MIPkit
git pull origin main
pip install -e .[dev] --upgrade
```

### Update Dependencies
```bash
pip install --upgrade -e .[dev]
```

## Uninstallation

```bash
pip uninstall MIPkit
```

## Testing Your Installation

Run the test suite to verify everything is working:

```bash
# Quick test
pytest tests/unit/test_utils.py -v

# Full test suite
pytest

# With coverage
pytest --cov=MIPkit --cov-report=html
```

## Getting Help

If you encounter issues:

1. Check the [GitHub Issues](https://github.com/tj-barrett/bmbt-mipkit/issues)
2. Read the [documentation](https://github.com/tj-barrett/bmbt-mipkit)
3. Ensure external tools (GROMACS, etc.) are properly installed
4. Verify all dependencies are installed: `pip list`

## Platform-Specific Notes

### Windows
- Use `venv\Scripts\activate` instead of `source venv/bin/activate`
- Some external tools may need manual installation
- Consider using WSL (Windows Subsystem for Linux) for better compatibility

### macOS
- Install Xcode Command Line Tools: `xcode-select --install`
- Use Homebrew for external tools when possible
- M1/M2 Macs may need rosetta or conda for some packages

### Linux
- Ubuntu/Debian: Use `apt-get` for system packages
- Ensure development headers are installed: `sudo apt-get install python3-dev`
- For HPC environments, use conda or module load
