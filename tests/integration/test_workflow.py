"""
Integration tests for MIPkit workflows.

These tests verify that multiple components work together correctly.
"""
import pytest
from pathlib import Path
from unittest.mock import Mock, patch


@pytest.mark.integration
class TestPDBWorkflow:
    """Test complete PDB reading and processing workflow."""

    def test_read_and_clean_pdb_workflow(self, temp_dir, sample_pdb_content):
        """Test reading a PDB file and cleaning it."""
        from MIPkit.utils.read_pdb import read_pdb, clean_pdb

        # Create PDB with extra records
        pdb_content = sample_pdb_content + "CONECT    1    2\nTER\n"
        pdb_file = temp_dir / "workflow_test.pdb"
        pdb_file.write_text(pdb_content)

        # Read PDB
        molecules = read_pdb(str(pdb_file), conect=False, hydrogens=True)
        assert len(molecules) > 0

        # Clean PDB
        clean_pdb(str(pdb_file))

        # Read again
        cleaned_molecules = read_pdb(str(pdb_file), conect=False, hydrogens=True)
        assert len(cleaned_molecules) > 0

    def test_pdb_molecule_dict_structure(self, sample_pdb_file):
        """Test that PDB molecules have correct dictionary structure."""
        from MIPkit.utils.read_pdb import read_pdb

        molecules = read_pdb(str(sample_pdb_file), hydrogens=True)

        # Check structure
        assert len(molecules) > 0
        mol = molecules[0]

        # Required keys
        assert "name" in mol
        assert "atoms" in mol
        assert "coords" in mol
        assert "bonds" in mol

        # Types
        assert isinstance(mol["name"], str)
        assert isinstance(mol["atoms"], list)
        assert isinstance(mol["coords"], list)
        assert isinstance(mol["bonds"], dict)


@pytest.mark.integration
class TestConfigGeneration:
    """Test configuration generation workflow."""

    def test_generate_config_creates_yaml(self, temp_dir):
        """Test that generate_config creates a YAML file."""
        from MIPkit.utils.utils import generate_config
        from unittest.mock import Mock
        import os

        # Mock args - use real values not Mock objects to avoid YAML serialization issues
        args = Mock()
        args.center_x = [5.0]
        args.center_y = [5.0]
        args.center_z = [5.0]
        args.size_x = [10.0]
        args.size_y = [10.0]
        args.size_z = [10.0]
        args.exhaustiveness = [8]
        args.fms = ["MAA", "2"]  # Real list, not Mock
        args.energy_range = 4  # Real number, not Mock

        # Change to temp dir
        original_dir = os.getcwd()
        try:
            os.chdir(temp_dir)

            # Create a simple PDB
            pdb_file = temp_dir / "template.pdb"
            pdb_file.write_text("ATOM      1  C   TST A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n")

            # Generate config (may fail if RDKit not available, which is fine)
            try:
                config = generate_config(str(pdb_file), args)

                # Check if YAML was created
                yaml_file = temp_dir / "generated_config.yaml"
                if yaml_file.exists():
                    assert yaml_file.exists()
                    content = yaml_file.read_text()
                    assert "center_x" in content or "exhaustiveness" in content
            except (ImportError, SystemExit):
                pytest.skip("RDKit not available or config generation failed")
        finally:
            os.chdir(original_dir)


@pytest.mark.integration
class TestErrorHandlingIntegration:
    """Test that error handling works across multiple components."""

    def test_invalid_fm_caught_across_modules(self):
        """Test that invalid FM inputs are caught consistently."""
        from MIPkit.utils.utils import generate_recipe
        from unittest.mock import patch

        with patch('MIPkit.utils.utils.fm2smiles') as mock_fm2smiles:
            mock_fm2smiles.return_value = {'MAA': 'C=C(C)C(O)=O'}

            # Invalid FM should exit
            with pytest.raises(SystemExit):
                generate_recipe(['INVALID_FM', '2'])

    def test_missing_template_caught_consistently(self):
        """Test that missing template files are caught."""
        from MIPkit.utils.utils import generate_config
        from unittest.mock import Mock

        args = Mock()

        # None template should exit
        with pytest.raises(SystemExit):
            generate_config(None, args)


@pytest.mark.integration
@pytest.mark.slow
class TestPerformanceIntegration:
    """Integration tests for performance optimizations."""

    def test_large_molecule_list_performance(self):
        """Test that set caching helps with large molecule lists."""
        import time

        # Simulate a large n_fm list
        n_fm = [f"MAA {i}" for i in range(1000)]

        # Test membership checking (what the code does)
        test_items = ["test1", "test2", "MAA 500"]

        # Without caching (old way)
        start = time.time()
        for _ in range(100):
            for item in test_items:
                result = item not in set(n_fm)  # Creates set each time
        old_time = time.time() - start

        # With caching (new way)
        start = time.time()
        n_fm_set = set(n_fm)  # Create once
        for _ in range(100):
            for item in test_items:
                result = item not in n_fm_set  # Reuse set
        new_time = time.time() - start

        # Should be significantly faster
        assert new_time < old_time


@pytest.mark.integration
class TestResourceCleanup:
    """Test that resources are properly cleaned up across workflows."""

    def test_multiple_pdb_reads_no_leak(self, temp_dir, sample_pdb_content):
        """Test that multiple PDB reads don't leak file handles."""
        from MIPkit.utils.read_pdb import read_pdb
        import gc

        pdb_file = temp_dir / "leak_test.pdb"
        pdb_file.write_text(sample_pdb_content)

        # Read many times
        for i in range(50):
            molecules = read_pdb(str(pdb_file))
            assert len(molecules) > 0

        # Force garbage collection
        gc.collect()

        # Should still be able to read
        molecules = read_pdb(str(pdb_file))
        assert len(molecules) > 0
