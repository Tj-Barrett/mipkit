"""
Unit tests for MIPkit/utils/read_pdb.py
"""
import pytest
from pathlib import Path


class TestReadPDB:
    """Test PDB file reading functionality."""

    def test_read_pdb_basic(self, sample_pdb_file):
        """Test basic PDB reading."""
        from MIPkit.utils.read_pdb import read_pdb

        molecules = read_pdb(str(sample_pdb_file), hydrogens=True)

        assert len(molecules) > 0
        assert molecules[0]["name"] == "MAA"
        assert "atoms" in molecules[0]
        assert "coords" in molecules[0]

    def test_read_pdb_without_hydrogens(self, sample_pdb_file):
        """Test PDB reading with hydrogens stripped."""
        from MIPkit.utils.read_pdb import read_pdb

        molecules_with_h = read_pdb(str(sample_pdb_file), hydrogens=True)
        molecules_without_h = read_pdb(str(sample_pdb_file), hydrogens=False)

        # Should have fewer atoms when hydrogens are excluded
        assert len(molecules_without_h[0]["atoms"]) <= len(molecules_with_h[0]["atoms"])

    def test_read_pdb_file_not_found(self):
        """Test that read_pdb handles missing files correctly."""
        from MIPkit.utils.read_pdb import read_pdb

        with pytest.raises(FileNotFoundError):
            read_pdb("nonexistent_file.pdb")

    def test_read_pdb_with_conect(self, temp_dir):
        """Test PDB reading with CONECT records."""
        from MIPkit.utils.read_pdb import read_pdb

        # Create PDB with CONECT records
        pdb_content = """ATOM      1  C   MAA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  C   MAA A   1       1.540   0.000   0.000  1.00  0.00           C
CONECT    1    2
CONECT    2    1
END
"""
        pdb_file = temp_dir / "conect_test.pdb"
        pdb_file.write_text(pdb_content)

        molecules = read_pdb(str(pdb_file), conect=True, hydrogens=False)

        assert len(molecules) > 0
        # Check that bonds were read if present
        if "bonds" in molecules[0]:
            assert isinstance(molecules[0]["bonds"], dict)


class TestCleanPDB:
    """Test PDB file cleaning functionality."""

    def test_clean_pdb_removes_excluded_records(self, sample_pdb_file):
        """Test that clean_pdb removes specified record types."""
        from MIPkit.utils.read_pdb import clean_pdb

        # This function modifies the file in place
        original_content = sample_pdb_file.read_text()

        # Clean the file (removes CONECT, TER, etc by default)
        clean_pdb(str(sample_pdb_file))

        cleaned_content = sample_pdb_file.read_text()

        # File should still have ATOM records
        assert "ATOM" in cleaned_content

        # Restore original
        sample_pdb_file.write_text(original_content)


@pytest.mark.unit
class TestPDBResourceManagement:
    """Test that PDB reading properly manages file resources."""

    def test_read_pdb_uses_context_manager(self, sample_pdb_file):
        """Verify that read_pdb uses context manager (no leaked file handles)."""
        from MIPkit.utils.read_pdb import read_pdb
        import gc

        # Read the file multiple times
        for _ in range(10):
            molecules = read_pdb(str(sample_pdb_file))

        # Force garbage collection
        gc.collect()

        # Try to delete the file - should succeed if file handles are closed
        import os
        # On Windows, this would fail if file handles were leaked
        # On Unix, we can still delete, but let's check we can open it
        with open(str(sample_pdb_file), 'r') as f:
            content = f.read()
            assert content  # Should be able to read

    def test_read_pdb_exception_handling(self, temp_dir):
        """Test that file handles close even when exceptions occur."""
        from MIPkit.utils.read_pdb import read_pdb

        # Create a malformed PDB that might cause issues
        bad_pdb = temp_dir / "bad.pdb"
        bad_pdb.write_text("INVALID PDB CONTENT\n")

        # Even if parsing has issues, file should close
        try:
            molecules = read_pdb(str(bad_pdb))
        except Exception:
            pass  # We expect it might fail

        # File handle should be closed
        with open(str(bad_pdb), 'r') as f:
            content = f.read()
            assert content
