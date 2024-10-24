import numpy as np
import pandas as pd
from typing import Tuple
from .config import Config

class SpotProcessor:
    def __init__(self):
        self.config = Config
    
    def calculate_intensities(self, spots_df: pd.DataFrame) -> pd.DataFrame:
        """Calculate intensities for each channel"""
        channels = self.config.get_round_spot_channels()
        for channel in channels:
            spots_df[f'chan_{channel}_intensity'] = (
                spots_df[f'chan_{channel}_fg'] - spots_df[f'chan_{channel}_bg']
            )
        return spots_df
    
    def filter_by_threshold(self, spots_df: pd.DataFrame) -> Tuple[pd.DataFrame, np.ndarray]:
        """Filter spots based on intensity threshold"""
        intensity_cols = [f'chan_{ch}_intensity' for ch in self.config.get_round_spot_channels()]
        
        # Find spots over threshold
        spots_over_thresh = spots_df.loc[
            np.any(spots_df[intensity_cols] > np.percentile(
                spots_df[intensity_cols],
                self.config.PERCENTILE, 0
            ), 1)
        ]['spot_id'].values
        
        # Mark spots over threshold
        spots_df['over_thresh'] = False
        spots_df.loc[spots_df['spot_id'].isin(spots_over_thresh), 'over_thresh'] = True
        
        return spots_df, spots_over_thresh
    
    def apply_qc_filters(self, spots_df: pd.DataFrame, stats_df: pd.DataFrame) -> pd.DataFrame:
        """Apply quality control filters"""
        filters = [
            spots_df['dist'] < self.config.CENT_CUTOFF,
            spots_df['r'] < self.config.CORR_CUTOFF,
            stats_df['dist_r'] > self.config.DIST_CUTOFF
        ]
        
        all_filters = np.all(np.vstack(filters), 0)
        return spots_df[all_filters]