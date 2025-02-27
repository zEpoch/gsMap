# test_gsmap.py
import logging
import shlex
import sys
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from gsMap.main import main


def parse_bash_command(command: str) -> list[str]:
    cleaned_command = command.replace("\\\n", " ")
    cleaned_command = " ".join(cleaned_command.splitlines())
    cleaned_command = " ".join(cleaned_command.split())
    return shlex.split(cleaned_command)


@pytest.mark.real_data
def test_gsmap_step_by_step_pipeline(example_data_dir, resource_dir, work_dir):
    """Test complete gsMap pipeline with real data"""
    logger = logging.getLogger("test_gsmap_pipeline")
    sample_name = "step_by_step_E16.5_E1S1"
    annotation = "annotation"
    data_layer = "count"
    homolog_file = resource_dir / "homologs/mouse_human_homologs.txt"
    gM_slices = None
    hdf5_path = example_data_dir / "ST/E16.5_E1S1.MOSTA.h5ad"
    trait_name = "IQ"
    sumstats_file = f'{example_data_dir}/GWAS/IQ_NG_2018.sumstats.gz'

    from gsMap.config import RunAllModeConfig

    run_all_config = RunAllModeConfig(
        workdir=work_dir,
        sample_name=sample_name,
        annotation=annotation,
        data_layer=data_layer,
        homolog_file=homolog_file,
        hdf5_path=str(hdf5_path),
        gM_slices=gM_slices,
        trait_name=trait_name,
        sumstats_file=str(sumstats_file),
        max_processes=10,
        gsMap_resource_dir=resource_dir,
    )

    # set bfile to subset of 50 individuals of 1000GEUR
    run_all_config.bfile_root = resource_dir / "LD_Reference_Panel_subset/1000G.EUR.QC.subset"

    # Step 1: Find latent representations
    logger.info("Step 1: Finding latent representations")
    command = f"""
    gsmap run_find_latent_representations \
        --workdir '{run_all_config.workdir}' \
        --sample_name {run_all_config.sample_name} \
        --input_hdf5_path '{run_all_config.hdf5_path}' \
        --annotation '{run_all_config.annotation}' \
        --data_layer '{run_all_config.data_layer}'
    """
    with patch.object(sys, "argv", parse_bash_command(command)):
        main()

    # Verify Step 1
    assert run_all_config.hdf5_with_latent_path.exists(), "Latent representation h5ad file not created"
    assert run_all_config.hdf5_with_latent_path.stat().st_size > 0, "Latent representation h5ad file is empty"

    # Step 2: Latent to gene
    logger.info("Step 2: Mapping latent representations to genes")
    command = f"""
    gsmap run_latent_to_gene \
        --workdir '{run_all_config.workdir}' \
        --sample_name {run_all_config.sample_name} \
        --annotation '{run_all_config.annotation}' \
        --latent_representation 'latent_GVAE' \
        --num_neighbour {run_all_config.num_neighbour} \
        --num_neighbour_spatial {run_all_config.num_neighbour_spatial} \
        --homolog_file '{run_all_config.homolog_file}'
    """
    with patch.object(sys, "argv", parse_bash_command(command)):
        main()

    # Verify Step 2
    assert run_all_config.mkscore_feather_path.exists(), "Mkscore feather file not created"
    assert run_all_config.mkscore_feather_path.stat().st_size > 0, "Mkscore feather file is empty"

    # Step 3: Generate LDScores
    logger.info("Step 3: Generating LDScores")
    command = f"""
    gsmap run_generate_ldscore \
        --workdir '{run_all_config.workdir}' \
        --sample_name {run_all_config.sample_name} \
        --chrom all \
        --bfile_root '{run_all_config.bfile_root}' \
        --keep_snp_root '{run_all_config.keep_snp_root}' \
        --gtf_annotation_file '{run_all_config.gtffile}' \
        --gene_window_size 50000
    """
    with patch.object(sys, "argv", parse_bash_command(command)):
        main()

    # Verify Step 3
    # verify ldscore files
    ldscore_chunk1_file = run_all_config.ldscore_save_dir / f"{run_all_config.sample_name}_chunk1" / f"{run_all_config.sample_name}.1.l2.ldscore.feather"
    assert pd.read_feather(ldscore_chunk1_file).shape[0] > 0, "LDScore chunk1 file is empty"
    assert run_all_config.ldscore_save_dir.is_dir(), "LDScore directory not created"
    assert any(run_all_config.ldscore_save_dir.iterdir()), "LDScore directory is empty"

    # Step 4: Spatial LDSC
    logger.info("Step 4: Running spatial LDSC")
    trait_name = "IQ"
    command = f"""
    gsmap run_spatial_ldsc \
        --workdir '{run_all_config.workdir}' \
        --sample_name {run_all_config.sample_name} \
        --trait_name '{trait_name}' \
        --sumstats_file '{run_all_config.sumstats_file}' \
        --w_file '{run_all_config.w_file}' \
        --num_processes {run_all_config.max_processes}
    """
    with patch.object(sys, "argv", parse_bash_command(command)):
        main()

    # Verify Step 4
    spatial_ldsc_result = run_all_config.get_ldsc_result_file(trait_name)
    assert spatial_ldsc_result.exists(), "Spatial LDSC results not created"
    assert spatial_ldsc_result.stat().st_size > 0, "Spatial LDSC results file is empty"

    # Step 5: Cauchy combination
    logger.info("Step 5: Running Cauchy combination test")
    command = f"""
    gsmap run_cauchy_combination \
        --workdir '{run_all_config.workdir}' \
        --sample_name {run_all_config.sample_name} \
        --trait_name '{trait_name}' \
        --annotation '{run_all_config.annotation}'
    """
    with patch.object(sys, "argv", parse_bash_command(command)):
        main()

    # Verify Step 5
    cauchy_result = run_all_config.get_cauchy_result_file(trait_name)
    assert cauchy_result.exists(), "Cauchy combination results not created"
    assert cauchy_result.stat().st_size > 0, "Cauchy combination results file is empty"

    # Step 6: Generate report
    logger.info("Step 6: Generating final report")
    command = f"""
    gsmap run_report \
        --workdir '{run_all_config.workdir}' \
        --sample_name {run_all_config.sample_name} \
        --trait_name '{trait_name}' \
        --annotation '{run_all_config.annotation}' \
        --sumstats_file '{run_all_config.sumstats_file}' \
        --top_corr_genes 50
    """
    with patch.object(sys, "argv", parse_bash_command(command)):
        main()

    # Verify Step 6
    report_file = run_all_config.get_gsMap_report_file(trait_name)
    assert report_file.exists(), "Final report not created"
    assert report_file.stat().st_size > 0, "Final report file is empty"

    # Verify report directory structure
    report_dir = run_all_config.get_report_dir(trait_name)
    assert report_dir.is_dir(), "Report directory not created"
    assert any(report_dir.iterdir()), "Report directory is empty"

    logger.info("Pipeline test completed successfully")


@pytest.mark.real_data
def test_gsmap_quick_mode(example_data_dir, resource_dir, work_dir):
    """Test the gsMap quick_mode pipeline with real data"""
    logger = logging.getLogger("test_gsmap_quick_mode")
    logger.info("Starting quick_mode pipeline test")

    # Setup parameters
    sample_name = "quick_mode_E16.5_E1S1"
    annotation = "annotation"
    data_layer = "count"
    homolog_file = resource_dir / "homologs/mouse_human_homologs.txt"
    hdf5_path = example_data_dir / "ST/E16.5_E1S1.MOSTA.h5ad"
    trait_name = "IQ"
    sumstats_file = f'{example_data_dir}/GWAS/IQ_NG_2018.sumstats.gz'

    # Test the quick_mode command
    command = f"""
    gsmap quick_mode \
        --workdir '{work_dir}' \
        --homolog_file '{homolog_file}' \
        --sample_name {sample_name} \
        --gsMap_resource_dir '{resource_dir}' \
        --hdf5_path '{hdf5_path}' \
        --annotation '{annotation}' \
        --data_layer '{data_layer}' \
        --sumstats_file '{sumstats_file}' \
        --trait_name '{trait_name}' \
        --max_processes 4
    """

    with patch.object(sys, "argv", parse_bash_command(command)):
        main()

    # Verify output files and directories
    from gsMap.config import RunAllModeConfig

    run_all_config = RunAllModeConfig(
        workdir=work_dir,
        sample_name=sample_name,
        annotation=annotation,
        data_layer=data_layer,
        homolog_file=str(homolog_file),
        hdf5_path=str(hdf5_path),
        trait_name=trait_name,
        sumstats_file=str(sumstats_file),
        max_processes=4,
        gsMap_resource_dir=resource_dir,
    )

    # Verify find_latent_representations step
    latent_rep_file = run_all_config.hdf5_with_latent_path
    assert latent_rep_file.exists(), "Latent representation h5ad file not created"
    assert latent_rep_file.stat().st_size > 0, "Latent representation h5ad file is empty"

    # Verify latent_to_gene step
    mkscore_file = run_all_config.mkscore_feather_path
    assert mkscore_file.exists(), "Mkscore feather file not created"
    assert mkscore_file.stat().st_size > 0, "Mkscore feather file is empty"

    # Verify generate_ldscore step (in quick mode, it uses symbolic links)
    ldscore_dir = run_all_config.ldscore_save_dir
    assert ldscore_dir.exists(), "LDScore directory not created"
    assert (ldscore_dir / "baseline").exists(), "Baseline annotation directory not created"
    assert (ldscore_dir / "SNP_gene_pair").exists(), "SNP_gene_pair directory not created"

    # Verify spatial_ldsc step
    spatial_ldsc_result = run_all_config.get_ldsc_result_file(trait_name)
    assert spatial_ldsc_result.exists(), "Spatial LDSC results not created"
    assert spatial_ldsc_result.stat().st_size > 0, "Spatial LDSC results file is empty"

    # Verify cauchy_combination step
    cauchy_result = run_all_config.get_cauchy_result_file(trait_name)
    assert cauchy_result.exists(), "Cauchy combination results not created"
    assert cauchy_result.stat().st_size > 0, "Cauchy combination results file is empty"

    # Verify report generation
    report_file = run_all_config.get_gsMap_report_file(trait_name)
    assert report_file.exists(), "Final report not created"
    assert report_file.stat().st_size > 0, "Final report file is empty"

    # Verify report directory structure
    report_dir = run_all_config.get_report_dir(trait_name)
    assert report_dir.is_dir(), "Report directory not created"
    assert any(report_dir.iterdir()), "Report directory is empty"

    # Verify all key visualizations are present
    gsmap_plot_dir = run_all_config.get_gsMap_plot_save_dir(trait_name)
    assert gsmap_plot_dir.exists(), "gsMap plot directory not created"
    assert any(gsmap_plot_dir.iterdir()), "gsMap plot directory is empty"

    manhattan_plot_path = run_all_config.get_manhattan_html_plot_path(trait_name)
    assert manhattan_plot_path.exists(), "Manhattan plot not created"

    gss_plot_dir = run_all_config.get_GSS_plot_dir(trait_name)
    assert gss_plot_dir.exists(), "GSS plot directory not created"
    assert any(gss_plot_dir.iterdir()), "GSS plot directory is empty"

    logger.info("Quick mode pipeline test completed successfully")

    # Setup parameters
    # sample_name = "E16.5_E1S1.MOSTA"
    # annotation = "annotation"
    # data_layer = "count"
    # homolog_file = resource_dir / "homologs/mouse_human_homologs.txt"
    # hdf5_path = example_data_dir / "ST/E16.5_E1S1.MOSTA.h5ad"

    # Create a config file for multiple traits
    config_file = work_dir / "gwas_config.yaml"
    with open(config_file, "w") as f:
        f.write(f"""
    IQ: {example_data_dir}/GWAS/IQ_NG_2018.sumstats.gz
    Height: {example_data_dir}/GWAS/GIANT_EUR_Height_2022_Nature.sumstats.gz
    MCHC: {example_data_dir}/GWAS/BCX2_MCHC_EA_GWAMA.sumstats.gz
    """)

    # Test the quick_mode command with config
    command = f"""
        gsmap quick_mode \
            --workdir '{work_dir}' \
            --homolog_file '{homolog_file}' \
            --sample_name {sample_name} \
            --gsMap_resource_dir '{resource_dir}' \
            --hdf5_path '{hdf5_path}' \
            --annotation '{annotation}' \
            --data_layer '{data_layer}' \
            --sumstats_config_file '{config_file}' \
            --max_processes 4
        """

    with patch.object(sys, "argv", parse_bash_command(command)):
        main()

    # Create run_all_config for verification
    from gsMap.config import RunAllModeConfig

    run_all_config = RunAllModeConfig(
        workdir=f"{work_dir}",
        sample_name=sample_name,
        annotation=annotation,
        data_layer=data_layer,
        homolog_file=str(homolog_file),
        hdf5_path=str(hdf5_path),
        sumstats_config_file=str(config_file),
        max_processes=4,
        gsMap_resource_dir=resource_dir,
    )

    # Verify results for each trait
    for trait_name in ["IQ", "Height", "MCHC"]:
        spatial_ldsc_result = run_all_config.get_ldsc_result_file(trait_name)
        assert spatial_ldsc_result.exists(), f"Spatial LDSC results for {trait_name} not created"

        cauchy_result = run_all_config.get_cauchy_result_file(trait_name)
        assert cauchy_result.exists(), f"Cauchy combination results for {trait_name} not created"

        report_file = run_all_config.get_gsMap_report_file(trait_name)
        assert report_file.exists(), f"Final report for {trait_name} not created"

    logger.info("Quick mode pipeline test with config file completed successfully")