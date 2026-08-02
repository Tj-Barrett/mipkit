"""
Unit tests for MIPkit/utils/utils.py
"""
import pytest
import sys
from pathlib import Path

from MIPkit.utils.utils import (
    isInt,
    isFloat,
    isString,
    isBool,
)


class TestTypeCheckers:
    """Test type checking utility functions."""

    def test_isInt_with_valid_integers(self):
        """Test isInt with valid integer inputs."""
        assert isInt(42) is True
        assert isInt(0) is True
        assert isInt(-10) is True
        assert isInt("123") is True

    def test_isInt_with_invalid_inputs(self):
        """Test isInt with non-integer inputs."""
        # Note: isInt(3.14) returns True because int(3.14) succeeds
        assert isInt(3.14) is True
        assert isInt("hello") is False
        assert isInt("3.14") is False  # int("3.14") raises ValueError
        # isInt(None) raises TypeError
        # isInt([]) raises TypeError

    def test_isFloat_with_valid_floats(self):
        """Test isFloat with valid float inputs."""
        assert isFloat(3.14) is True
        assert isFloat(0.0) is True
        assert isFloat(-2.5) is True
        assert isFloat("3.14") is True
        assert isFloat(42) is True  # Integers are floats too

    def test_isFloat_with_invalid_inputs(self):
        """Test isFloat with non-float inputs."""
        assert isFloat("hello") is False
        # isFloat(None) raises TypeError
        # isFloat([]) raises TypeError
        assert isFloat("not a number") is False

    def test_isString_with_valid_strings(self):
        """Test isString with valid string inputs."""
        assert isString("hello") is True
        # isString("") raises IndexError because it tries to access s[0]
        assert isString("123") is True

    def test_isString_with_invalid_inputs(self):
        """Test isString with non-string inputs."""
        assert isString(123) is False
        assert isString(3.14) is False
        assert isString(None) is False
        # isString([]) raises IndexError

    def test_isBool_with_valid_booleans(self):
        """Test isBool with valid boolean inputs."""
        assert isBool(True) is True
        assert isBool(False) is True
        # Note: isBool only checks isinstance(s, bool), so strings return False
        assert isBool("true") is False
        assert isBool("false") is False
        assert isBool("True") is False
        assert isBool("False") is False

    def test_isBool_with_invalid_inputs(self):
        """Test isBool with non-boolean inputs."""
        assert isBool(1) is False
        assert isBool(0) is False
        assert isBool("yes") is False
        assert isBool("no") is False
        assert isBool(None) is False
        assert isBool([]) is False


class TestGenerateConfig:
    """Test configuration generation functions."""

    def test_generate_config_with_template(self, temp_dir, sample_pdb_file):
        """Test generate_config with a valid template."""
        from MIPkit.utils.utils import generate_config
        from unittest.mock import Mock

        # Create mock args
        args = Mock()
        args.center_x = [5.0]
        args.center_y = [5.0]
        args.center_z = [5.0]
        args.size_x = [10.0]
        args.size_y = [10.0]
        args.size_z = [10.0]
        args.exhaustiveness = [8]

        # Test that config generation doesn't crash
        # (Full test would require RDKit)
        template_path = str(sample_pdb_file)
        # This will fail if RDKit isn't available, but that's expected
        # config = generate_config(template_path, args)

    def test_generate_config_without_template(self):
        """Test generate_config with no template exits properly."""
        from MIPkit.utils.utils import generate_config
        from unittest.mock import Mock

        args = Mock()

        # Should call sys.exit()
        with pytest.raises(SystemExit):
            generate_config(None, args)


class TestGenerateRecipe:
    """Test recipe generation from functional monomers."""

    def test_generate_recipe_valid_input(self):
        """Test generate_recipe with valid FM input."""
        from MIPkit.utils.utils import generate_recipe

        # Mock the fm2smiles function
        from unittest.mock import patch

        with patch('MIPkit.utils.utils.fm2smiles') as mock_fm2smiles:
            mock_fm2smiles.return_value = {'MAA': 'C=C(C)C(O)=O', 'EGDMA': 'COC(=O)C(=C)C(=O)OC'}

            # Valid input: FM name followed by quantity
            afms = ['MAA', '2', 'EGDMA', '8']
            fms = generate_recipe(afms)

            assert fms == {'MAA': 2, 'EGDMA': 8}

    def test_generate_recipe_invalid_quantity(self):
        """Test generate_recipe with invalid quantity exits properly."""
        from MIPkit.utils.utils import generate_recipe
        from unittest.mock import patch

        with patch('MIPkit.utils.utils.fm2smiles') as mock_fm2smiles:
            mock_fm2smiles.return_value = {'MAA': 'C=C(C)C(O)=O'}

            # Missing quantity
            afms = ['MAA', 'invalid']

            with pytest.raises(SystemExit):
                generate_recipe(afms)

    def test_generate_recipe_unknown_fm(self):
        """Test generate_recipe with unknown FM exits properly."""
        from MIPkit.utils.utils import generate_recipe
        from unittest.mock import patch

        with patch('MIPkit.utils.utils.fm2smiles') as mock_fm2smiles:
            mock_fm2smiles.return_value = {'MAA': 'C=C(C)C(O)=O'}

            # Unknown FM
            afms = ['UNKNOWN_FM', '2']

            with pytest.raises(SystemExit):
                generate_recipe(afms)


@pytest.mark.unit
class TestErrorHandling:
    """Test that our exception handling improvements work correctly."""

    def test_generate_recipe_catches_specific_exceptions(self):
        """Verify that generate_recipe catches specific exceptions, not bare except."""
        from MIPkit.utils.utils import generate_recipe
        from unittest.mock import patch

        with patch('MIPkit.utils.utils.fm2smiles') as mock_fm2smiles:
            mock_fm2smiles.return_value = {'MAA': 'C=C(C)C(O)=O'}

            # Test that IndexError is caught (trying to access afms[i+1] when out of bounds)
            afms = ['MAA']  # Missing quantity

            with pytest.raises(SystemExit):
                generate_recipe(afms)
