# MIPkit Test Suite

Comprehensive test suite for MIPkit, covering unit tests, integration tests, and performance validations.

## Test Structure

```
tests/
├── unit/                          # Unit tests for individual modules
│   ├── test_utils.py             # Tests for utility functions
│   ├── test_read_pdb.py          # Tests for PDB reading/writing
│   ├── test_parser.py            # Tests for argument parsing
│   ├── test_run_gmx.py           # Tests for GROMACS operations
│   └── test_performance_optimizations.py  # Tests for optimizations
├── integration/                   # Integration tests
│   └── test_workflow.py          # End-to-end workflow tests
├── fixtures/                      # Test data and fixtures
├── conftest.py                   # Shared pytest fixtures
└── README.md                     # This file
```

## Running Tests

### Install Test Dependencies

```bash
pip install pytest pytest-cov pytest-mock
```

### Run All Tests

```bash
# From the project root
pytest

# With coverage report
pytest --cov=MIPkit --cov-report=html

# Verbose output
pytest -v
```

### Run Specific Test Categories

```bash
# Unit tests only
pytest tests/unit/

# Integration tests only
pytest tests/integration/

# Tests with specific markers
pytest -m unit
pytest -m integration
pytest -m "not slow"
```

### Run Specific Test Files

```bash
# Test utilities
pytest tests/unit/test_utils.py

# Test PDB reading
pytest tests/unit/test_read_pdb.py

# Test a specific test function
pytest tests/unit/test_utils.py::TestTypeCheckers::test_isInt_with_valid_integers
```

## Test Markers

Tests are marked with the following pytest markers:

- `@pytest.mark.unit` - Unit tests for individual functions
- `@pytest.mark.integration` - Tests that combine multiple components
- `@pytest.mark.slow` - Tests that take longer to run
- `@pytest.mark.requires_rdkit` - Tests that need RDKit installed
- `@pytest.mark.requires_external` - Tests needing external tools (gmx, obabel, etc.)

## What's Tested

### Critical Bug Fixes
- ✅ String comparison bug (`is not` → `!=`)
- ✅ Bare exit statement fix (`exit` → `sys.exit()`)
- ✅ Bare exception handlers (specific exceptions)
- ✅ Unreachable code removal
- ✅ Resource leak fixes (context managers)

### API Modernization
- ✅ `subprocess.call()` → `subprocess.run()`
- ✅ `os.system()` → `subprocess.run()`
- ✅ Proper `shell=True` parameter usage

### Performance Optimizations
- ✅ Set caching in loops
- ✅ Variable naming improvements
- ✅ Performance benchmarks

### Functionality
- ✅ Type checking functions (isInt, isFloat, isString, isBool)
- ✅ PDB file reading and writing
- ✅ Configuration generation
- ✅ Recipe generation from functional monomers
- ✅ Error handling and validation

## Coverage Goals

Target coverage: **>80%** for core utility modules

Current focus areas:
- `MIPkit/utils/utils.py` - Core utilities
- `MIPkit/utils/read_pdb.py` - PDB I/O
- `MIPkit/utils/parser.py` - Argument parsing
- `MIPkit/run/run_gmx.py` - GROMACS operations

## Writing New Tests

### Test Naming Convention

- Test files: `test_*.py` or `*_test.py`
- Test classes: `Test*`
- Test functions: `test_*`

### Example Test

```python
import pytest

class TestMyFeature:
    """Test suite for my feature."""

    def test_basic_functionality(self):
        """Test that basic functionality works."""
        result = my_function(input_data)
        assert result == expected_output

    def test_error_handling(self):
        """Test that errors are handled properly."""
        with pytest.raises(ValueError):
            my_function(invalid_input)

    @pytest.mark.slow
    def test_performance(self):
        """Test performance characteristics."""
        import time
        start = time.time()
        result = my_function(large_dataset)
        elapsed = time.time() - start
        assert elapsed < 1.0  # Should complete in under 1 second
```

### Using Fixtures

```python
def test_with_temp_file(temp_dir):
    """Test using the temp_dir fixture."""
    test_file = temp_dir / "test.txt"
    test_file.write_text("test content")
    assert test_file.exists()

def test_with_sample_pdb(sample_pdb_file):
    """Test using the sample PDB file fixture."""
    from MIPkit.utils.read_pdb import read_pdb
    molecules = read_pdb(str(sample_pdb_file))
    assert len(molecules) > 0
```

## Continuous Integration

To run tests in CI/CD:

```bash
# Install dependencies
pip install -e .[dev]
pip install pytest pytest-cov

# Run tests with XML output for CI
pytest --cov=MIPkit --cov-report=xml --junitxml=junit.xml

# Generate coverage badge
coverage-badge -o coverage.svg
```

## Troubleshooting

### Tests fail with ImportError

```bash
# Make sure MIPkit is installed
pip install -e .

# Or add to PYTHONPATH
export PYTHONPATH="${PYTHONPATH}:$(pwd)"
```

### RDKit-dependent tests fail

Some tests require RDKit. Either:
1. Install RDKit: `conda install -c conda-forge rdkit`
2. Skip RDKit tests: `pytest -m "not requires_rdkit"`

### External tool tests fail

Tests requiring gmx, obabel, vina can be skipped:
```bash
pytest -m "not requires_external"
```

## Contributing

When adding new features:
1. Write tests for new functionality
2. Ensure existing tests still pass: `pytest`
3. Check coverage: `pytest --cov=MIPkit --cov-report=term-missing`
4. Aim for >80% coverage on new code

## Test Maintenance

- Review and update tests when fixing bugs
- Add regression tests for reported issues
- Keep test data minimal but representative
- Document complex test scenarios
- Remove obsolete tests when refactoring
