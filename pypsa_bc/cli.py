#!/usr/bin/env python3

import argparse
import subprocess
import os
import sys
from pathlib import Path
from typing import Optional
from bc_combined_modelling.clews.runner import RunModel


class Clews:
    def __init__(self, clews_builder_cfg_path: Path, combined_model_cfg_path: Path):
        self.clews_builder_cfg_path = clews_builder_cfg_path.resolve()
        self.combined_model_cfg_path = combined_model_cfg_path.resolve()

    def update(self, keyword: str):
        base_dirs = [
            Path('workflow/scripts'),
            Path('models/BC_Nexus/scripts'),
        ]
        script_map = {
            'run': 'runner.py',
        }
        script_name = script_map.get(keyword)
        if not script_name:
            print(f"Unknown keyword. Supported keywords: {list(script_map.keys())}")
            return
        script_path = self.find_script(base_dirs, script_name)
        if script_path:
            self.run_script(script_path)
        else:
            print(f"Script for keyword '{keyword}' not found.")
    
    def run(self):
        # Instantiate and run the model
        args={
        'storage_algorithm': 'Kotzur',
        'scenario':'Base',
        'input_csvs':'data/clews_data/Model_Kotzur/inputs_csv',
        'combined_model_config_path':'config/config.yaml'
        }
        bcnexus = RunModel(**args)
        bcnexus.run(build=True)

class LinkingTool:
    def __init__(self, combined_model_cfg_path: Path):
        self.combined_model_cfg_path = combined_model_cfg_path.resolve()

    def update(self, keyword: str):
        base_dirs = [Path('models/Linking_tool/workflow/scripts')]
        script_map = {'run': 'run.py'}
        script_name = script_map.get(keyword)
        if not script_name:
            print(f"Unknown keyword. Supported keywords: {list(script_map.keys())}")
            return
        script_path = self.find_script(base_dirs, script_name)
        if script_path:
            self.run_script(script_path)
        else:
            print(f"Script for keyword '{keyword}' not found.")

    @staticmethod
    def find_script(base_dirs, script_name: str) -> Optional[Path]:
        for base_dir in base_dirs:
            script_path = base_dir / script_name
            if script_path.exists():
                return script_path.resolve()
        return None

    def run_script(self, script_path: Path):
        try:
            subprocess.run(['python', str(script_path), str(self.combined_model_cfg_path)], check=True)
            print(f">>>> Successfully executed {script_path}.")
        except subprocess.CalledProcessError as e:
            print(f"Error executing {script_path}: {e}")

@staticmethod
def find_project_root(marker="bc_combined_modelling"):
    """
    Traverse upwards to find the root directory containing a specific folder or file.
    """
    current_dir = Path.cwd()
    while current_dir != current_dir.parent:  # Stop at the filesystem root
        if (current_dir / marker).exists():  # Check for the marker
            return current_dir
        current_dir = current_dir.parent
    raise FileNotFoundError(f"Project root with marker '{marker}' not found.")

def main():
    try:
        # Find and set the project root dynamically
        root_dir = find_project_root("bc_combined_modelling")  # Adjust marker as needed
        os.chdir(root_dir)
        print(f"Running command from project root: {root_dir}")
        
    except FileNotFoundError as e:
        print(e)
        sys.exit(1)
        
    parser = argparse.ArgumentParser(
        description='BC Combined Model CLI: A tool for managing models and configurations.',
        epilog='Examples:\n  bccm clews run \n  bccm linkingtool run',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # CLEWS subcommand
    clews_parser = subparsers.add_parser('clews', help='Manage CLEWS models and configurations.')
    clews_parser.add_argument(
        'keyword',
        choices=['update_tech_schema', 'update_profiling_params', 'update_yearly_params', 'run'],
        help='Specify the operation to perform.',
    )
    clews_parser.add_argument(
        '--combined-model-config',
        type=Path,
        default=Path('config/config.yaml'),
        help='Path to the combined model configuration file.',
    )
    clews_parser.add_argument(
        '--scenario',
        type=str,
        default='Base',
        help='Scenario name harmonized with configurations',
    )
    clews_parser.add_argument(
        '--storage_algorithm',
        type=str,
        default='Kotzur',
        help='Storage algorithm name (e.g. Kotzur,Niet)',
    )
    clews_parser.add_argument(
        '--input_csvs',
        type=Path,
        default=Path('data/clews_data/Model_Kotzur/inputs_csv'),
        help='Path to the input csvs',
    )
    # LinkingTool subcommand
    linking_parser = subparsers.add_parser('linkingtool', help='Manage LinkingTool operations.')
    linking_parser.add_argument(
        'keyword',
        choices=['run'],
        help='Specify the operation to perform.',
    )
    linking_parser.add_argument(
        '--combined-model-config',
        type=Path,
        default=Path('config/config.yaml'),
        help='Path to the combined model configuration file.',
    )

    args = parser.parse_args()

    if args.command == 'clews':
        model = Clews(args.clews_builder_config, args.combined_model_config)
        
        # Apply the Method based on keyword
        if args.keyword!='run':
            model.run()
        else:
            pass
            
    elif args.command == 'linkingtool':
        model = LinkingTool(args.combined_model_config)
        model.update(args.keyword)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
