# test_gsmap.py
import logging
import shlex
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

from gsMap.main import main


def parse_bash_command(command: str) -> list[str]:
    cleaned_command = command.replace("\\\n", " ")
    cleaned_command = " ".join(cleaned_command.splitlines())
    cleaned_command = " ".join(cleaned_command.split())
    return shlex.split(cleaned_command)


@pytest.mark.real_data
def test_gsmap_pipeline(example_data_dir, resource_dir, work_dir):
    """Test complete gsMap pipeline with real data"""
    logger = logging.getLogger("test_gsmap_pipeline")
    sample_name = "E16.5_E1S1.MOSTA"
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
        --chrom 1 \
        --bfile_root '{run_all_config.bfile_root}' \
        --keep_snp_root '{run_all_config.keep_snp_root}' \
        --gtf_annotation_file '{run_all_config.gtffile}' \
        --gene_window_size 50000
    """
    with patch.object(sys, "argv", parse_bash_command(command)):
        main()

    # Verify Step 3
    ldscore_done_file = run_all_config.ldscore_save_dir / "generate_ldscore.done"
    assert ldscore_done_file.exists(), "LDScore generation not completed"
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
