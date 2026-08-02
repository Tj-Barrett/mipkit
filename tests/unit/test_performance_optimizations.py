"""
Unit tests for performance optimizations we implemented.
"""
import pytest
from unittest.mock import Mock, patch
import time


class TestSetCachingOptimization:
    """Test the set caching optimization in complex.py."""

    def test_set_caching_reduces_set_creations(self):
        """Verify that set is created once, not in every iteration."""
        from MIPkit.run import complex as complex_module
        import inspect

        source = inspect.getsource(complex_module)

        # Look for the pattern: n_fm_set = set(n_fm)
        # This should appear before the if statements
        assert 'n_fm_set = set(n_fm)' in source

        # Count how many times we create the set vs how many times we use it
        set_creation_count = source.count('n_fm_set = set(n_fm)')
        set_usage_count = source.count('not in n_fm_set')

        # We should use the cached set multiple times
        assert set_usage_count > set_creation_count

    def test_no_inline_set_conversion_in_loops(self):
        """Verify that we don't have set(n_fm) inside conditional checks."""
        from MIPkit.run import complex as complex_module
        import inspect

        source = inspect.getsource(complex_module)

        # Should not have "not in set(n_fm)" pattern
        # (we replaced it with "not in n_fm_set")
        problematic_pattern = 'not in set(n_fm)'

        if problematic_pattern in source:
            pytest.fail("Found inline set conversion in loop - not using cached set")

    @pytest.mark.slow
    def test_set_caching_performance_benefit(self):
        """Demonstrate that set caching improves performance."""
        # Simulate the old way vs new way

        # Old way: create set in loop
        n_fm = ['MAA 1', 'EGDMA 2', 'MAA 3', 'EGDMA 4'] * 100  # 400 items

        start = time.time()
        for i in range(1000):
            ta = 'test'
            # Old way - creates set every iteration
            result = ta not in set(n_fm)
        old_time = time.time() - start

        # New way: cache set
        start = time.time()
        n_fm_set = set(n_fm)  # Create once
        for i in range(1000):
            ta = 'test'
            # New way - use cached set
            result = ta not in n_fm_set
        new_time = time.time() - start

        # New way should be significantly faster
        # Allow some variance, but should be at least 2x faster
        assert new_time < old_time * 0.8, f"Cached version ({new_time:.4f}s) not faster than uncached ({old_time:.4f}s)"


class TestVariableNamingImprovements:
    """Test that variable naming improvements were applied."""

    def test_line_variable_renamed_from_underscore_l(self):
        """Verify that _l was renamed to line in read_pdb."""
        from MIPkit.utils import read_pdb
        import inspect

        source = inspect.getsource(read_pdb)

        # Should have "for line_idx, line in enumerate(lines):"
        assert 'for line_idx, line in enumerate(lines):' in source

        # Should not have confusing _l variable
        # Allow _l in comments or strings, but not as a variable name
        lines = [line.strip() for line in source.split('\n')
                if not line.strip().startswith('#') and '_l' in line]

        # Filter out acceptable uses
        problematic_lines = [line for line in lines
                           if 'for' in line and '_l in enumerate' in line]

        if problematic_lines:
            pytest.fail(f"Found confusing _l variable in: {problematic_lines}")


@pytest.mark.unit
class TestCodeQualityMetrics:
    """Test overall code quality improvements."""

    def test_no_bare_except_in_utils(self):
        """Verify no bare except clauses in utils.py."""
        from MIPkit.utils import utils
        import inspect
        import re

        source = inspect.getsource(utils)

        # Look for bare except patterns
        bare_except_pattern = r'except\s*:\s*\n'
        matches = list(re.finditer(bare_except_pattern, source))

        if matches:
            pytest.fail(f"Found {len(matches)} bare except clauses in utils.py")

    def test_no_bare_except_in_complex(self):
        """Verify no bare except clauses in complex.py."""
        from MIPkit.run import complex as complex_module
        import inspect
        import re

        source = inspect.getsource(complex_module)

        bare_except_pattern = r'except\s*:\s*\n'
        matches = list(re.finditer(bare_except_pattern, source))

        if matches:
            pytest.fail(f"Found {len(matches)} bare except clauses in complex.py")

    def test_specific_exception_types_used(self):
        """Verify that specific exception types are used."""
        from MIPkit.utils import utils
        import inspect

        source = inspect.getsource(utils)

        # Should have specific exception types
        assert any(exc in source for exc in ['IndexError', 'ValueError', 'TypeError', 'KeyError'])

    def test_context_managers_for_files(self):
        """Verify that file operations use context managers."""
        from MIPkit.utils import utils
        import inspect

        source = inspect.getsource(utils)

        # Count with open() statements
        with_open_count = source.count('with open(')

        # Count bare open() statements (not in with)
        # This is approximate but gives us an idea
        bare_open_lines = [line for line in source.split('\n')
                          if '= open(' in line and 'with open(' not in line]

        # Most file operations should use context managers
        if len(bare_open_lines) > with_open_count:
            pytest.fail(f"More bare open() calls ({len(bare_open_lines)}) than with open() ({with_open_count})")
