import argparse
import glob
import os
import logging
import shutil
import subprocess

from sbmlpbkutils import run_config, load_config, plot_simulation_results

CONFIGS_PATH = './validation/scenarios/'
OUTPUT_PATH = './Output'
R_CONFIGS = []

# Configure logger for formatted console output
console_logger = logging.getLogger('create_simulation_reports')
console_logger.setLevel(logging.INFO)
_console_handler = logging.StreamHandler()
_console_handler.setLevel(logging.INFO)
_console_handler.setFormatter(logging.Formatter('[%(levelname)s] %(message)s'))
if not console_logger.handlers:
    console_logger.addHandler(_console_handler)

def main():
    parser = argparse.ArgumentParser(description="Run validation scripts arguments.")
    parser.add_argument(
        '-f',
        '--force_recompute',
        action='store_true',
        help="Force re-calculation of validation results instead of using cached results."
    )

    args = parser.parse_args()
    force_recompute = args.force_recompute

    # Configure logger for formatted console output
    logger = logging.getLogger('run_validation')
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter('[%(levelname)s] %(message)s'))
    if not logger.handlers:
        logger.addHandler(console_handler)

    # If force recompute then remove output
    if force_recompute:
        shutil.rmtree(OUTPUT_PATH)

    # Ensure output path
    os.makedirs(OUTPUT_PATH, exist_ok=True)

    # Run R validaton scenarios
    run_r_validation_scenarios(force_recompute, logger)

    # Glob scenario configs and run them
    configs = glob.glob(f'{CONFIGS_PATH}/**/*.yaml', recursive=True)
    for file in configs:
        file_dir = os.path.dirname(file)

        # Load config
        config = load_config(file)

        # Create output directory if it does not exist
        out_path = os.path.join(OUTPUT_PATH, os.path.relpath(file_dir, CONFIGS_PATH), config.id)

        # Ensure config output path
        os.makedirs(out_path, exist_ok=True)

        # Run simulations
        logger.info("Running simulation config %s", file)
        run_config(
            config = config,
            out_path = out_path,
            force_recompute = force_recompute,
            logger = logger
        )

        # Run simulations
        plot_simulation_results(
            config = config,
            out_path = out_path,
            combine_outputs = True
        )

def run_r_validation_scenarios(
    force_recompute: bool,
    logger: logging.Logger
):
    for r_config in R_CONFIGS:
        out_files = r_config['output_files']
        if not force_recompute and all([os.path.exists(o) for o in out_files]):
            logger.info(f"Skipping {r_config['id']} validation scenarios: results already available")
            return

        # Run R validation scenarios
        logger.info("Running R validation scenarios")
        subprocess.run(['Rscript', r_config['file_path']], check=False)

if __name__ == "__main__":
    main()
