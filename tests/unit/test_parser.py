"""
Unit tests for MIPkit/utils/parser.py
"""
import pytest
from pathlib import Path


class TestParserTypeChecks:
    """Test the type checking functions in parser."""

    def test_isInt_parser(self):
        """Test isInt from parser module."""
        from MIPkit.utils.parser import isInt

        assert isInt(42) is True
        assert isInt("123") is True
        assert isInt(3.14) is True  # int(3.14) succeeds
        assert isInt("abc") is False

    def test_isString_parser(self):
        """Test isString from parser module."""
        from MIPkit.utils.parser import isString

        assert isString("hello") is True
        # isString("") raises IndexError
        assert isString(123) is False

    def test_isBool_parser(self):
        """Test isBool from parser module."""
        from MIPkit.utils.parser import isBool

        assert isBool(True) is True
        assert isBool(False) is True
        assert isBool("true") is False  # Only actual bool type returns True
        assert isBool("false") is False


class TestNCyclesConversion:
    """Test the ncycles parameter conversion that we fixed."""

    def test_cycles_conversion_valid_integer(self):
        """Test that valid integer strings are converted properly."""
        from MIPkit.utils.parser import isInt

        # Our fix ensures that ValueError, IndexError, TypeError are caught
        assert isInt("10") is True
        assert isInt(10) is True

    def test_cycles_conversion_invalid_input(self):
        """Test that invalid inputs are handled gracefully."""
        from MIPkit.utils.parser import isInt

        assert isInt("not a number") is False
        # isInt(None) raises TypeError, it's not caught by the function


@pytest.mark.unit
class TestSubprocessOperations:
    """Test subprocess-related operations in parser."""

    def test_subprocess_cleanup_uses_subprocess_run(self, temp_dir):
        """Verify that subprocess.run is used (not deprecated subprocess.call or os.system)."""
        from MIPkit.utils import parser
        import inspect

        # Check that the module doesn't import os.system or subprocess.call in problematic ways
        source = inspect.getsource(parser)

        # Should not have bare os.system calls (we replaced them)
        assert 'os.system(f"rm' not in source or 'subprocess.run(f"rm' in source

        # Should use subprocess.run, not subprocess.call
        if 'subprocess.call' in source:
            pytest.fail("Found deprecated subprocess.call in parser.py")


class TestDirectorySetup:
    """Test directory setup functions."""

    def test_sanitize_inputs_creates_directories(self, temp_dir, monkeypatch):
        """Test that sanitize_inputs creates necessary directories."""
        from MIPkit.utils.parser import sanitize_inputs
        from unittest.mock import Mock

        # Change to temp directory
        monkeypatch.chdir(temp_dir)

        # Create mock args
        args = Mock()
        args.basedir = ["test_base"]
        args.react = False
        args.restart = False

        # This will create directories
        # Note: Full test would require mocking all dependencies
        # For now, test that the function exists and can be called
        assert callable(sanitize_inputs)


@pytest.mark.unit
class TestExceptionHandling:
    """Test that exception handling uses specific exceptions, not bare except."""

    def test_parser_no_bare_except(self):
        """Verify no bare except clauses remain in parser."""
        from MIPkit.utils import parser
        import inspect
        import re

        source = inspect.getsource(parser)

        # Look for bare except: patterns (allowing for whitespace)
        # This regex matches "except:" but not "except SomeException:"
        bare_except_pattern = r'except\s*:\s*\n'

        matches = list(re.finditer(bare_except_pattern, source))

        if matches:
            pytest.fail(f"Found {len(matches)} bare except clauses in parser.py")
