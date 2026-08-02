"""
Pytest configuration and shared fixtures for MIPkit tests.
"""
import os
import sys
import tempfile
import pytest
from pathlib import Path

# Add MIPkit to path
sys.path.insert(0, str(Path(__file__).parent.parent))


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_pdb_content():
    """Sample PDB file content for testing."""
    return """ATOM      1  C   MAA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  C   MAA A   1       1.540   0.000   0.000  1.00  0.00           C
ATOM      3  O   MAA A   1       2.240   1.000   0.000  1.00  0.00           O
ATOM      4  O   MAA A   1       2.240  -1.220   0.000  1.00  0.00           O
ATOM      5  H   MAA A   1      -0.360   1.026   0.000  1.00  0.00           H
ATOM      6  H   MAA A   1      -0.360  -0.513  -0.889  1.00  0.00           H
ATOM      7  H   MAA A   1      -0.360  -0.513   0.889  1.00  0.00           H
ATOM      8  H   MAA A   1       3.200  -1.220   0.000  1.00  0.00           H
END
"""


@pytest.fixture
def sample_pdb_file(temp_dir, sample_pdb_content):
    """Create a sample PDB file for testing."""
    pdb_file = temp_dir / "test.pdb"
    pdb_file.write_text(sample_pdb_content)
    return pdb_file


@pytest.fixture
def mock_pargs():
    """Mock pargs object for testing."""
    class MockPargs:
        def __init__(self):
            self.gmx = "gmx"
            self.obabel = "obabel"
            self.vina = "vina"
            self.temp = 300
            self.cutoff_CC = 2.5
            self.cutoff_CO = 3.0
            self.dt = {"min": None, "nvt": 0.001, "npt": 0.001, "rxn": 0.001}
            self.nsteps = {"min": 10000, "nvt": 50000, "npt": 50000, "rxn": 50000}
            self.clean = False
            self.explicit = False
            self.implicit = True

    return MockPargs()
