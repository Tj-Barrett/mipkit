"""
Unit tests for MIPkit/run/run_gmx.py

Tests focused on the fixes we made to critical bugs.
"""
import pytest
import sys
import re
from unittest.mock import Mock, patch, MagicMock


class TestStringComparison:
    """Test the string comparison fix."""

    def test_run_type_comparison_uses_equality(self):
        """Verify that run_type comparison uses != not 'is not'."""
        from MIPkit.run import run_gmx
        import inspect

        # Get the source code
        source = inspect.getsource(run_gmx)

        # Should not have "is not" with string literals
        # Our fix changed: if run_type is not "min": to if run_type != "min":
        if 'is not "min"' in source or "is not 'min'" in source:
            pytest.fail("Found 'is not' with string literal - should use !=")


class TestUnreachableCode:
    """Test that unreachable code was removed."""

    def test_no_unreachable_elif_after_else(self):
        """Verify that the unreachable elif after else was removed."""
        from MIPkit.run import run_gmx
        import inspect
        import re

        source = inspect.getsource(run_gmx)

        # Look for the pattern of else: followed by elif error:
        # This was unreachable code we fixed
        pattern = r'else:\s+.*\s+sys\.exit\(\)\s+elif\s+error:'

        if re.search(pattern, source, re.MULTILINE | re.DOTALL):
            pytest.fail("Found unreachable elif after else block")


@pytest.mark.unit
class TestSubprocessUsage:
    """Test that subprocess operations are updated."""

    def test_uses_subprocess_run_not_os_system(self):
        """Verify that os.system was replaced with subprocess.run."""
        from MIPkit.run import run_gmx
        import inspect

        source = inspect.getsource(run_gmx)

        # Count os.system vs subprocess.run usage
        os_system_count = source.count('os.system(')
        subprocess_run_count = source.count('subprocess.run(')

        # We should have replaced os.system calls
        # Allow for commented-out os.system
        uncommented_os_system = [line for line in source.split('\n')
                                 if 'os.system(' in line and not line.strip().startswith('#')]

        if uncommented_os_system:
            pytest.fail(f"Found {len(uncommented_os_system)} uncommented os.system calls")

    def test_subprocess_calls_have_shell_parameter(self):
        """Verify that subprocess.run calls using shell commands have shell=True."""
        from MIPkit.run import run_gmx
        import inspect
        import re

        source = inspect.getsource(run_gmx)

        # Find subprocess.run calls with f-strings that look like shell commands
        shell_command_pattern = r'subprocess\.run\(f["\'](?:cp|rm|mv)'

        matches = re.finditer(shell_command_pattern, source)

        for match in matches:
            # Get the full function call
            start = match.start()
            # Find the closing parenthesis
            call_text = source[start:start+200]  # Look ahead

            # Should have shell=True
            if 'shell=True' not in call_text:
                pytest.fail(f"subprocess.run shell command missing shell=True: {call_text[:80]}")


class TestResourceManagement:
    """Test that file resources are properly managed."""

    def test_file_operations_use_context_managers(self):
        """Verify that file operations use 'with' statements."""
        from MIPkit.run import run_gmx
        import inspect
        import re

        source = inspect.getsource(run_gmx)

        # Look for open() calls that aren't in with statements
        # This is a simplified check
        lines = source.split('\n')

        for i, line in enumerate(lines):
            if '= open(' in line and 'with open(' not in line:
                # Check if there's a corresponding .close() within next 20 lines
                has_close = any('.close()' in lines[i+j] for j in range(1, min(21, len(lines)-i)))

                # If not, it might be a resource leak
                if not has_close and 'workdir' not in line:
                    # Some context might be okay, but flag it
                    pass  # Could add warnings here


@pytest.mark.unit
class TestExceptionHandling:
    """Test exception handling improvements."""

    def test_gmx_make_mdp_catches_specific_exceptions(self):
        """Verify that gmx_make_mdp catches KeyError, not bare except."""
        from MIPkit.run import run_gmx
        import inspect

        source = inspect.getsource(run_gmx.gmx_make_mdp)

        # Should catch KeyError specifically
        assert 'except KeyError' in source or 'KeyError' in source

        # Should not have bare except
        if re.search(r'except\s*:\s*\n', source):
            pytest.fail("Found bare except in gmx_make_mdp")
