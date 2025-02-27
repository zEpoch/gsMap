# conftest.py
import subprocess
from pathlib import Path

import pytest


def pytest_addoption(parser):
    """Add custom command line options"""
    parser.addoption(
        "--run-real-data",
        action="store_true",
        default=False,
        help="run tests that require real data",
    )
    parser.addoption(
        "--example-data-dir",
        action="store",
        default=None,
        help="Path to existing gsMap example data directory",
    )
    parser.addoption(
        "--resource-dir",
        action="store",
        default=None,
        help="Path to existing gsMap resource directory",
    )


@pytest.fixture(scope="session")
def example_data_dir(request, tmp_path_factory):
    """Get or create example data directory"""
    example_dir = request.config.getoption("--example-data-dir")
    if example_dir:
        return Path(example_dir)

    # If no directory provided, download and set up in temporary directory
    base_dir = tmp_path_factory.mktemp("gsmap_test")
    example_data_url = "https://yanglab.westlake.edu.cn/data/gsMap/gsMap_example_data.tar.gz"

    # Download and extract example data
    subprocess.run(["wget", example_data_url], cwd=base_dir, check=True)
    subprocess.run(["tar", "-xvzf", "gsMap_example_data.tar.gz"], cwd=base_dir, check=True)

    return base_dir / "gsMap_example_data"


@pytest.fixture(scope="session")
def resource_dir(request, tmp_path_factory):
    """Get or create resource directory"""
    resource_dir = request.config.getoption("--resource-dir")
    if resource_dir:
        return Path(resource_dir)

    # If no directory provided, download and set up in temporary directory
    base_dir = tmp_path_factory.mktemp("gsmap_resource")
    resource_url = "https://yanglab.westlake.edu.cn/data/gsMap/gsMap_resource.tar.gz"

    # Download and extract resource files
    subprocess.run(["wget", resource_url], cwd=base_dir, check=True)
    subprocess.run(["tar", "-xvzf", "gsMap_resource.tar.gz"], cwd=base_dir, check=True)

    return base_dir / "gsMap_resource"


@pytest.fixture(scope="session")
def work_dir(tmp_path_factory):
    """Create working directory for test outputs"""
    work_dir = tmp_path_factory.mktemp("Mouse_Embryo")
    return work_dir




def pytest_configure(config):
    """Configure pytest with custom markers and default options"""
    # Add marker descriptions
    config.addinivalue_line(
        "markers", "real_data: mark test that requires real gsMap data (disabled by default)"
    )

    # Set default options (equivalent to addopts in pytest.ini)
    if not config.option.verbose:
        config.option.verbose = True


def pytest_collection_modifyitems(config, items):
    """Skip real_data tests by default unless --run-real-data is specified"""
    if not config.getoption("--run-real-data"):
        skip_real_data = pytest.mark.skip(reason="need --run-real-data option to run")
        for item in items:
            if "real_data" in item.keywords:
                item.add_marker(skip_real_data)
