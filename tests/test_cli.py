# test_cli.py
import logging
import shlex
import sys
from unittest.mock import patch

import pytest

from gsMap.main import main


def parse_bash_command(command: str) -> list[str]:
    """Convert multi-line bash command to argument list for sys.argv"""
    cleaned_command = command.replace("\\\n", " ")
    cleaned_command = " ".join(cleaned_command.splitlines())
    cleaned_command = " ".join(cleaned_command.split())
    return shlex.split(cleaned_command)


@pytest.mark.real_data
def test_gsmap_step_by_step_pipeline(stepbystep_config):
    """Test complete gsMap pipeline with real data"""
    logger = logging.getLogger("test_gsmap_pipeline")
    config = stepbystep_config
    logger.info(f"Running step-by-step pipeline test with config: {config.sample_name}")

    # Step 1: Find latent representations
    logger.info("Step 1: Finding latent representations")
    command = f"""
    gsmap run_find_latent_representations \
        --workdir '{config.workdir}' \
        --sample_name {config.sample_name} \
        --input_hdf5_path '{config.hdf5_path}' \
        --annotation '{config.annotation}' \
        --data_layer '{config.data_layer}' \
        --n_comps '{config.n_comps}'
    """
    with patch.object(sys, "argv", parse_bash_command(command)):
        main()

    # Verify Step 1
    assert config.hdf5_with_latent_path.exists(), "Latent representation h5ad file not created"
    assert config.hdf5_with_latent_path.stat().st_size > 0, (
        "Latent representation h5ad file is empty"
    )

    # Step 2: Latent to gene
    logger.info("Step 2: Mapping latent representations to genes")
    command = f"""
    gsmap run_latent_to_gene \
        --workdir '{config.workdir}' \
        --sample_name {config.sample_name} \
        --annotation '{config.annotation}' \
        --latent_representation 'latent_GVAE' \
        --num_neighbour {config.num_neighbour} \
        --num_neighbour_spatial {config.num_neighbour_spatial} \
        --homolog_file '{config.homolog_file}'
    """
    with patch.object(sys, "argv", parse_bash_command(command)):
        main()

    # Verify Step 2
    assert config.mkscore_feather_path.exists(), "Mkscore feather file not created"
    assert config.mkscore_feather_path.stat().st_size > 0, "Mkscore feather file is empty"

    # Step 3: Generate LDScores
    logger.info("Step 3: Generating LDScores")
    command = f"""
    gsmap run_generate_ldscore \
        --workdir '{config.workdir}' \
        --sample_name {config.sample_name} \
        --chrom all \
        --bfile_root '{config.bfile_root}' \
        --keep_snp_root '{config.keep_snp_root}' \
        --gtf_annotation_file '{config.gtffile}' \
        --gene_window_size 50000
    """
    with patch.object(sys, "argv", parse_bash_command(command)):
        main()

    # Verify Step 3
    # verify ldscore files
    ldscore_chunk1_file = (
        config.ldscore_save_dir
        / f"{config.sample_name}_chunk1"
        / f"{config.sample_name}.22.l2.ldscore.feather"
    )
    assert ldscore_chunk1_file.exists(), "LDScore chunk1 file not created"
    assert ldscore_chunk1_file.stat().st_size > 0, "LDScore chunk1 file is empty"
    assert config.ldscore_save_dir.is_dir(), "LDScore directory not created"
    assert any(config.ldscore_save_dir.iterdir()), "LDScore directory is empty"

    # Step 4: Spatial LDSC
    logger.info("Step 4: Running spatial LDSC")
    command = f"""
    gsmap run_spatial_ldsc \
        --workdir '{config.workdir}' \
        --sample_name {config.sample_name} \
        --trait_name '{config.trait_name}' \
        --sumstats_file '{config.sumstats_file}' \
        --w_file '{config.w_file}' \
        --num_processes {config.max_processes}
    """
    with patch.object(sys, "argv", parse_bash_command(command)):
        main()

    # Verify Step 4
    spatial_ldsc_result = config.get_ldsc_result_file(config.trait_name)
    assert spatial_ldsc_result.exists(), "Spatial LDSC results not created"
    assert spatial_ldsc_result.stat().st_size > 0, "Spatial LDSC results file is empty"

    # Step 5: Cauchy combination
    logger.info("Step 5: Running Cauchy combination test")
    command = f"""
    gsmap run_cauchy_combination \
        --workdir '{config.workdir}' \
        --sample_name {config.sample_name} \
        --trait_name '{config.trait_name}' \
        --annotation '{config.annotation}'
    """
    with patch.object(sys, "argv", parse_bash_command(command)):
        main()

    # Verify Step 5
    cauchy_result = config.get_cauchy_result_file(config.trait_name)
    assert cauchy_result.exists(), "Cauchy combination results not created"
    assert cauchy_result.stat().st_size > 0, "Cauchy combination results file is empty"

    # Step 6: Generate report
    logger.info("Step 6: Generating final report")
    command = f"""
    gsmap run_report \
        --workdir '{config.workdir}' \
        --sample_name {config.sample_name} \
        --trait_name '{config.trait_name}' \
        --annotation '{config.annotation}' \
        --sumstats_file '{config.sumstats_file}' \
        --top_corr_genes 50
    """
    with patch.object(sys, "argv", parse_bash_command(command)):
        main()

    # Verify Step 6
    report_file = config.get_gsMap_report_file(config.trait_name)
    assert report_file.exists(), "Final report not created"
    assert report_file.stat().st_size > 0, "Final report file is empty"

    # Verify report directory structure
    report_dir = config.get_report_dir(config.trait_name)
    assert report_dir.is_dir(), "Report directory not created"
    assert any(report_dir.iterdir()), "Report directory is empty"

    logger.info("Pipeline test completed successfully")


@pytest.mark.real_data
def test_gsmap_quick_mode(quickmode_config):
    """Test the gsMap quick_mode pipeline with real data"""
    logger = logging.getLogger("test_gsmap_quick_mode")
    logger.info("Starting quick_mode pipeline test")
    config = quickmode_config

    # Test the quick_mode command
    command = f"""
    gsmap quick_mode \
        --workdir '{config.workdir}' \
        --homolog_file '{config.homolog_file}' \
        --sample_name {config.sample_name} \
        --gsMap_resource_dir '{config.gsMap_resource_dir}' \
        --hdf5_path '{config.hdf5_path}' \
        --annotation '{config.annotation}' \
        --data_layer '{config.data_layer}' \
        --sumstats_file '{config.sumstats_file}' \
        --trait_name '{config.trait_name}' \
        --max_processes {config.max_processes}
    """

    with patch.object(sys, "argv", parse_bash_command(command)):
        main()

    # Verify output files and directories

    # Verify find_latent_representations step
    latent_rep_file = config.hdf5_with_latent_path
    assert latent_rep_file.exists(), "Latent representation h5ad file not created"
    assert latent_rep_file.stat().st_size > 0, "Latent representation h5ad file is empty"

    # Verify latent_to_gene step
    mkscore_file = config.mkscore_feather_path
    assert mkscore_file.exists(), "Mkscore feather file not created"
    assert mkscore_file.stat().st_size > 0, "Mkscore feather file is empty"

    # Verify generate_ldscore step (in quick mode, it uses symbolic links)
    ldscore_dir = config.ldscore_save_dir
    assert ldscore_dir.exists(), "LDScore directory not created"
    assert (ldscore_dir / "baseline").exists(), "Baseline annotation directory not created"
    assert (ldscore_dir / "SNP_gene_pair").exists(), "SNP_gene_pair directory not created"

    # Verify spatial_ldsc step
    spatial_ldsc_result = config.get_ldsc_result_file(config.trait_name)
    assert spatial_ldsc_result.exists(), "Spatial LDSC results not created"
    assert spatial_ldsc_result.stat().st_size > 0, "Spatial LDSC results file is empty"

    # Verify cauchy_combination step
    cauchy_result = config.get_cauchy_result_file(config.trait_name)
    assert cauchy_result.exists(), "Cauchy combination results not created"
    assert cauchy_result.stat().st_size > 0, "Cauchy combination results file is empty"

    # Verify report generation
    report_file = config.get_gsMap_report_file(config.trait_name)
    assert report_file.exists(), "Final report not created"
    assert report_file.stat().st_size > 0, "Final report file is empty"

    # Verify report directory structure
    report_dir = config.get_report_dir(config.trait_name)
    assert report_dir.is_dir(), "Report directory not created"
    assert any(report_dir.iterdir()), "Report directory is empty"

    # Verify all key visualizations are present
    gsmap_plot_dir = config.get_gsMap_plot_save_dir(config.trait_name)
    assert gsmap_plot_dir.exists(), "gsMap plot directory not created"
    assert any(gsmap_plot_dir.iterdir()), "gsMap plot directory is empty"

    manhattan_plot_path = config.get_manhattan_html_plot_path(config.trait_name)
    assert manhattan_plot_path.exists(), "Manhattan plot not created"

    gss_plot_dir = config.get_GSS_plot_dir(config.trait_name)
    assert gss_plot_dir.exists(), "GSS plot directory not created"
    assert any(gss_plot_dir.iterdir()), "GSS plot directory is empty"

    logger.info("Quick mode pipeline test completed successfully")
