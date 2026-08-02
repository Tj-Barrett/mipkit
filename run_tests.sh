#!/bin/bash
# Test runner script for MIPkit

set -e  # Exit on error

echo "======================================"
echo "MIPkit Test Suite"
echo "======================================"
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if pytest is installed
if ! command -v pytest &> /dev/null; then
    echo -e "${RED}Error: pytest is not installed${NC}"
    echo "Install with: pip install pytest pytest-cov pytest-mock"
    exit 1
fi

# Default: run all tests
TEST_SCOPE="${1:-all}"

case "$TEST_SCOPE" in
    "unit")
        echo -e "${YELLOW}Running unit tests only...${NC}"
        pytest tests/unit/ -v
        ;;
    "integration")
        echo -e "${YELLOW}Running integration tests only...${NC}"
        pytest tests/integration/ -v
        ;;
    "fast")
        echo -e "${YELLOW}Running fast tests (excluding slow tests)...${NC}"
        pytest -m "not slow" -v
        ;;
    "coverage")
        echo -e "${YELLOW}Running tests with coverage report...${NC}"
        pytest --cov=MIPkit --cov-report=term-missing --cov-report=html
        echo ""
        echo -e "${GREEN}Coverage report generated in htmlcov/index.html${NC}"
        ;;
    "all")
        echo -e "${YELLOW}Running all tests...${NC}"
        pytest -v
        ;;
    *)
        echo -e "${RED}Unknown test scope: $TEST_SCOPE${NC}"
        echo "Usage: $0 [unit|integration|fast|coverage|all]"
        exit 1
        ;;
esac

# Check exit code
if [ $? -eq 0 ]; then
    echo ""
    echo -e "${GREEN}✓ All tests passed!${NC}"
else
    echo ""
    echo -e "${RED}✗ Some tests failed${NC}"
    exit 1
fi
