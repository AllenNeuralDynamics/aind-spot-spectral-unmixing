#!/usr/bin/env python3
"""
Main script for processing multichannel spot data.
This script demonstrates how to use the spot_analysis package to process
a complete round of multichannel data.
"""

import logging
from pathlib import Path
from typing import List, Optional
import argparse

import numpy as np
import pandas as pd

from spot_analysis.config import Config
from spot_analysis.data_loader import SpotDataLoader
from spot_analysis.spot_processor import SpotProcessor
from spot_analysis.ratio_calculator import RatioCalculator
from spot_analysis.unmixer import SpotUnmixer
from spot_analysis.cell_by_gene_table import cell_by_gene_processor


class SpotAnalysisPipeline:
    def __init__(
        self,
        spots_folder: Path = Path('/data/'),
        output_folder: Path = Path('/results/'),
        min_distances: Optional[List[float]] = None,
        config = None  # Lowercase parameter name for clarity
    ):
        # Use instance-only config
        if config is None: 
            self.config = Config()  # Create new instance
        else: 
            self.config = config  # Use passed instance
        
        # Update INSTANCE attributes only (no class modification!)
        self.config.SPOTS_FOLDER = spots_folder
        self.config.OUTPUT_FOLDER = output_folder
        
        # Create directories using instance config
        self.config.OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)
        self.config.SCRATCH_FOLDER.mkdir(parents=True, exist_ok=True)
        
        self.min_distances = min_distances or [3.0, 4.0, 5.0]
        
        # Initialize components
        self.data_loader = SpotDataLoader()
        self.processor = SpotProcessor()
        self.ratio_calculator = RatioCalculator()
        self.unmixer = SpotUnmixer(self.config)
        self.cell_by_gene_processor = cell_by_gene_processor()
        
        # Setup logging
        self.logger = self._setup_logging()
        
    def _setup_logging(self) -> logging.Logger:
        """Setup logging configuration"""
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.INFO)

        # Clear existing handlers to avoid duplicates when processing multiple tiles
        if logger.handlers:
            for handler in list(logger.handlers):
                logger.removeHandler(handler)
                handler.close()
        
        # Create handlers
        console_handler = logging.StreamHandler()
        tile_suffix = f'_tile_{self.config.CURRENT_TILE}' if self.config.CURRENT_TILE else ''
        log_path = self.config.OUTPUT_FOLDER / f'round_{self.config.ROUND_N}{tile_suffix}_processing.log'
        file_handler = logging.FileHandler(log_path)
        
        # Create formatters and add it to handlers
        log_format = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        console_handler.setFormatter(log_format)
        file_handler.setFormatter(log_format)
        
        # Add handlers to the logger
        logger.addHandler(console_handler)
        logger.addHandler(file_handler)
        
        return logger
    
    def run(self):
        """Run the complete spot analysis pipeline"""
        tile_suffix = f'_tile_{self.config.CURRENT_TILE}' if self.config.CURRENT_TILE else ''
        tile_label = self.config.CURRENT_TILE or 'N/A'
        self.logger.info(f"Starting analysis for Round {self.config.ROUND_N}, Tile {tile_label}")
        
        try:
            # 1. Load Data
            self.logger.info("Loading spot data...")
            spots_df = self.data_loader.load_all_spots()
            self.logger.info(f"Loaded {len(spots_df)} total spots")
            
            # 2. Calculate intensities and apply threshold filtering
            self.logger.info("Processing spots...")
            spots_df = self.processor.calculate_intensities(spots_df)
            spots_df, spots_over_thresh = self.processor.filter_by_threshold(spots_df)
            
            self.logger.info(
                f"Found {len(spots_over_thresh)} spots over threshold "
                f"({len(spots_over_thresh)/len(spots_df)*100:.1f}%)"
            )
            
            # Save intermediate results
            mixed_output_path = self.config.OUTPUT_FOLDER / f'mixed_spots_R{self.config.ROUND_N}{tile_suffix}.pkl'
            mixed_scratch_path = self.config.SCRATCH_FOLDER / f'mixed_spots_R{self.config.ROUND_N}{tile_suffix}.pkl'
            spots_df.to_pickle(mixed_output_path)
            spots_df.to_pickle(mixed_scratch_path)
            
            # 3. Calculate ratios
            self.logger.info("Calculating channel ratios...")
            ratio_path = self.config.OUTPUT_FOLDER / f'r{self.config.ROUND_N}{tile_suffix}_ratios.txt'
            intensity_cols = [
                f'chan_{ch}_intensity'
                for ch in self.config.get_round_spot_channels()
            ]
            # ratios = self.ratio_calculator.calculate_ratios( #tim's old way
            #     spots_df[intensity_cols].values,
            #     ratio_path
            # )
            ratios = self.ratio_calculator.calculate_ratios_by_channel( #trying a weighted avg
                spots_df[intensity_cols].values,
                ratio_path, 
                spots_df['chan'].to_numpy()
            )

            # 4. Calculate distances and create stats
            self.logger.info("Calculating distances between spots...")
            #thresh_spots_df = spots_df[spots_df['over_thresh']].copy()
            stats_df = self.unmixer.calculate_distances(spots_df, ratios)

            #save stats_df
            stats_df_csv_path = self.config.OUTPUT_FOLDER / f'spot_unmixing_stats{tile_suffix}.csv'
            stats_df.to_csv(stats_df_csv_path)
            
            # 5. Application of QC filters has been moved to interactive capsule
            self.logger.info("Applying QC filters...")
            spots_df = self.processor.apply_qc_filters(
                spots_df,
                stats_df
            )
            
            #self.logger.info(
            #    f"Kept {len(filtered_spots_df)} spots after QC "
            #    f"({len(filtered_spots_df)/len(thresh_spots_df)*100:.1f}%)"
            #)
            all_chans_filt_stats = self.unmixer.calculate_distances(spots_df, ratios)
            
            # 6. Process multiple minimum distances
            self.logger.info("Running unmixing with multiple minimum distances...")
            results = self.unmixer.process_multiple_distances(
                spots_df,
                all_chans_filt_stats,
                self.min_distances
            )
            
            # 7. Save summary statistics
            self._save_summary_statistics(results)
            
            

            # 8. Generate cleaned up tables for analysis 
            # processor = cell_by_gene_processor()
            # unmixed_results, mixed_results = processor.process_pipeline([Config.ROUND_N])
            # Save final results
            #unmixed_results.to_csv(
            #    Config.OUTPUT_FOLDER / 
            #    f'unmixed_spots_R{Config.ROUND_N}.csv'
            #)
            #unmixed_results.to_pickle(
            #    Config.OUTPUT_FOLDER / 
            #    f'unmixed_spots_R{Config.ROUND_N}.pkl'
            #)
            #mixed_results.to_csv(
            #    Config.OUTPUT_FOLDER / 
            #    f'mixed_spots_R{Config.ROUND_N}.csv'
            #)
            self.logger.info("Processing completed successfully")

            return results
            
        except Exception as e:
            self.logger.error(f"Error during processing: {str(e)}", exc_info=True)
            raise
            
    def _save_summary_statistics(self, results):
        """Save summary statistics for all processing runs"""
        summary_stats = []
        tile_suffix = f'_tile_{self.config.CURRENT_TILE}' if self.config.CURRENT_TILE else ''
        
        for min_dist, (unmixed_df, stats) in results.items():
            for channel_stat in stats:
                stat_dict = {
                    'min_dist': min_dist,
                    'round': self.config.ROUND_N,
                    'tile': self.config.CURRENT_TILE,
                    **channel_stat
                }
                summary_stats.append(stat_dict)
        
        # Create summary DataFrame and save
        summary_df = pd.DataFrame(summary_stats)
        summary_df.to_csv(
            self.config.OUTPUT_FOLDER / 
            f'round_{self.config.ROUND_N}{tile_suffix}_summary_stats.csv',
            index=False
        )


def main():
    """Main entry point"""
    # Example usage
    pipeline = SpotAnalysisPipeline(
       spots_folder=Path('/root/capsule/data/'),
       output_folder=Path('/root/capsule/results/'),
       min_distances=[3.0]
    )

    # parser = argparse.ArgumentParser()
    # parser.add_argument("round", type=str, help="Dataset round until procedures.json is finished")
    # args = parser.parse_args()
    # round = int(args.round)


    # pipeline = SpotAnalysisPipeline(round_number = round, min_distances = [3])
    
    results = pipeline.run()
    
    # Print final statistics
    print("\nFinal Results Summary:")
    for min_dist, (unmixed_df, stats) in results.items():
        print(f"\nResults for minimum distance {min_dist}:")
        for channel_stat in stats:
            print(
                f"Channel {channel_stat['channel']} ({channel_stat['gene']}): "
                f"Kept {channel_stat['kept_spots']} of {channel_stat['total_spots']} spots "
                f"({channel_stat['kept_spots']/channel_stat['total_spots']*100:.1f}%)"
            )


if __name__ == '__main__':
    main()