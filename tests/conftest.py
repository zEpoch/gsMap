import logging
from pathlib import Path
import copy

import pytest

from gsMap.config import RunAllModeConfig


def pytest_addoption(parser):
    """Add custom command line options"""
    parser.addoption(
        "--run-real-data",
        action="store_true",
        default=False,
        help="run tests that require real data",
    )
    parser.addoption(
        "--test-data",
        action="store",
        default=None,
        help="Path to test data directory for gsMap tests",
    )
    parser.addoption(
        "--work-dir",
        action="store",
        default=None,
        help="Path to working directory for test outputs (defaults to a temporary directory)",
    )


@pytest.fixture(scope="session")
def test_data_dir(request):
    """Get test data directory"""
    test_dir = request.config.getoption("--test-data")
    if not test_dir:
        pytest.skip("--test-data not provided")

    test_dir = Path(test_dir)
    if not test_dir.exists():
        pytest.skip(f"Test data directory does not exist: {test_dir}")

    return test_dir


@pytest.fixture(scope="session")
def example_data_dir(test_data_dir):
    """Get path to example data directory"""
    example_dir = test_data_dir / "gsMap_example_data"
    if not example_dir.exists():
        pytest.skip(f"Example data directory not found: {example_dir}")

    return example_dir


@pytest.fixture(scope="session")
def resource_dir(test_data_dir):
    """Get path to resource directory"""
    resource_dir = test_data_dir / "gsMap_resource"
    if not resource_dir.exists():
        pytest.skip(f"Resource directory not found: {resource_dir}")

    return resource_dir


@pytest.fixture(scope="session")
def work_dir(request, tmp_path_factory):
    """
    If --work-dir is provided, use that directory instead of a temporary one.
    """
    custom_dir = request.config.getoption("--work-dir")

    if custom_dir:
        work_dir = Path(custom_dir)
        work_dir.mkdir(parents=True, exist_ok=True)
        return work_dir

    # Otherwise, use a temporary directory
    return tmp_path_factory.mktemp("Mouse_Embryo")


@pytest.fixture(scope="session")
def homolog_file(resource_dir):
    """Get path to mouse to human homolog file"""
    homolog_path = resource_dir / "homologs/mouse_human_homologs.txt"
    if not homolog_path.exists():
        pytest.skip(f"Homolog file not found: {homolog_path}")
    return homolog_path


@pytest.fixture(scope="session")
def subsampled_h5ad_file1(example_data_dir):
    """Get path to first subsampled h5ad file"""
    h5ad_path = example_data_dir / "ST/E16.5_E1S1.MOSTA_subsampled.h5ad"
    if not h5ad_path.exists():
        pytest.skip(f"Subsampled h5ad file not found: {h5ad_path}")
    return h5ad_path


@pytest.fixture(scope="session")
def subsampled_h5ad_file2(example_data_dir):
    """Get path to second subsampled h5ad file"""
    h5ad_path = example_data_dir / "ST/E16.5_E2S11.MOSTA_subsampled.h5ad"
    if not h5ad_path.exists():
        pytest.skip(f"Subsampled h5ad file not found: {h5ad_path}")
    return h5ad_path


@pytest.fixture(scope="session")
def iq_sumstats_file(example_data_dir):
    """Get path to IQ GWAS summary statistics file"""
    sumstats_path = example_data_dir / "GWAS/IQ_NG_2018.sumstats.gz"
    if not sumstats_path.exists():
        pytest.skip(f"Summary statistics file not found: {sumstats_path}")
    return sumstats_path


@pytest.fixture(scope="session")
def reference_panel(resource_dir):
    """Get path to reference panel file"""
    ref_panel_path_prefix = (resource_dir / "LD_Reference_Panel_subset/1000G.EUR.QC.subset").as_posix()
    for chromosome in range(1, 23):
        bim_file = Path(f"{ref_panel_path_prefix}.{chromosome}.bim")
        assert bim_file.exists(), f"Reference panel file not found: {bim_file}"
    return ref_panel_path_prefix


@pytest.fixture(scope="session")
def base_config(work_dir, resource_dir, subsampled_h5ad_file1, homolog_file, iq_sumstats_file, reference_panel):
    """Create a base RunAllModeConfig fixture"""
    basic_config = RunAllModeConfig(
        workdir=work_dir,
        sample_name="test_sample",  # This will be overridden in specific test fixtures
        annotation="annotation",
        data_layer="count",
        homolog_file=str(homolog_file),
        hdf5_path=str(subsampled_h5ad_file1),
        trait_name="IQ",
        sumstats_file=str(iq_sumstats_file),
        max_processes=4,
        gsMap_resource_dir=str(resource_dir),
    )
    basic_config.bfile_root = reference_panel
    return basic_config


@pytest.fixture
def stepbystep_config(base_config):
    """Create a config for step-by-step pipeline tests"""
    config = copy.deepcopy(base_config)
    config.sample_name = "step_by_step_test"

    return config


@pytest.fixture
def quickmode_config(base_config):
    """Create a config for quick mode tests"""
    config = copy.deepcopy(base_config)
    config.sample_name = "quick_mode_test"
    return config


@pytest.fixture
def conditional_config(base_config):
    """Create a config for conditional analysis tests"""
    config = copy.deepcopy(base_config)
    config.sample_name = "conditional_analysis_test"
    return config


@pytest.fixture
def customlatent_config(base_config):
    """Create a config for customized latent tests"""
    config = copy.deepcopy(base_config)
    config.sample_name = "custom_latent_test"
    return config


@pytest.fixture
def biorep_config1(base_config, subsampled_h5ad_file1):
    """Create a config for first biological replicate"""
    config = copy.deepcopy(base_config)
    config.sample_name = "bio_rep_test_1"
    config.hdf5_path = str(subsampled_h5ad_file1)
    return config


@pytest.fixture
def biorep_config2(base_config, subsampled_h5ad_file2):
    """Create a config for second biological replicate"""
    config = copy.deepcopy(base_config)
    config.sample_name = "bio_rep_test_2"
    config.hdf5_path = str(subsampled_h5ad_file2)
    return config


@pytest.fixture
def additional_baseline_dir(test_data_dir):
    """Get path to additional baseline annotations directory"""
    baseline_dir = test_data_dir / "gsMap_additional_annotation"
    if not baseline_dir.exists():
        pytest.skip(f"Additional baseline annotation directory not found: {baseline_dir}")
    return baseline_dir


def pytest_configure(config):
    """Configure pytest with custom markers and default options"""
    # Add marker descriptions
    config.addinivalue_line(
        "markers", "real_data: mark test that requires real data (disabled by default)"
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