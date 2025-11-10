"""
Module for merging single-tile spot tables into combined multi-tile tables.
"""
import pandas as pd
from pathlib import Path
from typing import List, Tuple, Optional
import logging
from .config import Config


class TableMerger:
    """Merges single-tile spot tables (both mixed and unmixed) into combined tables."""
    
    def __init__(self):
        self.config = Config
        self.logger = logging.getLogger(__name__)
    
    def find_tile_tables(
        self, 
        folder: Path, 
        pattern: str,
        file_extension: str = '.pkl'
    ) -> List[Path]:
        """
        Find all tile-specific tables matching a pattern.
        
        Args:
            folder: Directory to search in
            pattern: Base filename pattern (e.g., 'unmixed_spots_R', 'unmixed_cell_by_gene')
            file_extension: File extension (.pkl or .csv)
            
        Returns:
            List of matching file paths
        """
        if not folder.exists():
            self.logger.warning(f"Folder does not exist: {folder}")
            return []
        
        # Find files matching pattern with tile suffix
        matching_files = []
        for file_path in folder.glob(f"{pattern}*_tile_*{file_extension}"):
            matching_files.append(file_path)
        
        self.logger.info(f"Found {len(matching_files)} files matching pattern '{pattern}*_tile_*{file_extension}'")
        return sorted(matching_files)
    
    def merge_pickle_tables(
        self, 
        file_paths: List[Path],
        output_path: Path
    ) -> Optional[pd.DataFrame]:
        """
        Merge multiple pickle files into a single dataframe.
        
        Args:
            file_paths: List of pickle file paths to merge
            output_path: Path where merged table should be saved
            
        Returns:
            Merged dataframe or None if no files to merge
        """
        if not file_paths:
            self.logger.warning("No files to merge")
            return None
        
        dataframes = []
        for file_path in file_paths:
            try:
                df = pd.read_pickle(file_path)
                tile_name = self._extract_tile_name(file_path)
                self.logger.info(f"Loaded {len(df)} rows from {file_path.name} (tile: {tile_name})")
                dataframes.append(df)
            except Exception as e:
                self.logger.error(f"Error loading {file_path}: {e}")
                continue
        
        if not dataframes:
            self.logger.warning("No dataframes successfully loaded")
            return None
        
        # Concatenate all dataframes
        merged_df = pd.concat(dataframes, ignore_index=True)
        self.logger.info(f"Merged {len(dataframes)} tables into {len(merged_df)} total rows")
        
        # Save merged table
        output_path.parent.mkdir(parents=True, exist_ok=True)
        merged_df.to_pickle(output_path)
        self.logger.info(f"Saved merged table to {output_path}")
        
        return merged_df
    
    def merge_csv_tables(
        self, 
        file_paths: List[Path],
        output_path: Path
    ) -> Optional[pd.DataFrame]:
        """
        Merge multiple CSV files into a single dataframe.
        
        Args:
            file_paths: List of CSV file paths to merge
            output_path: Path where merged table should be saved
            
        Returns:
            Merged dataframe or None if no files to merge
        """
        if not file_paths:
            self.logger.warning("No files to merge")
            return None
        
        dataframes = []
        for file_path in file_paths:
            try:
                df = pd.read_csv(file_path, index_col=0)
                tile_name = self._extract_tile_name(file_path)
                self.logger.info(f"Loaded {len(df)} rows from {file_path.name} (tile: {tile_name})")
                dataframes.append(df)
            except Exception as e:
                self.logger.error(f"Error loading {file_path}: {e}")
                continue
        
        if not dataframes:
            self.logger.warning("No dataframes successfully loaded")
            return None
        
        # Concatenate all dataframes
        merged_df = pd.concat(dataframes, ignore_index=True)
        self.logger.info(f"Merged {len(dataframes)} tables into {len(merged_df)} total rows")
        
        # Save merged table
        output_path.parent.mkdir(parents=True, exist_ok=True)
        merged_df.to_csv(output_path)
        self.logger.info(f"Saved merged table to {output_path}")
        
        return merged_df
    
    def _extract_tile_name(self, file_path: Path) -> str:
        """Extract tile name from filename."""
        # Assumes filename pattern: *_tile_<tile_name>_*.ext or *_tile_<tile_name>.ext
        parts = file_path.stem.split('_tile_')
        if len(parts) > 1:
            # Get everything after '_tile_' and remove any trailing suffixes
            tile_part = parts[1]
            # Remove trailing parts like '_minDist_3'
            if '_minDist_' in tile_part:
                tile_part = tile_part.split('_minDist_')[0]
            return tile_part
        return "unknown"
    
    def merge_all_spot_tables(
        self,
        round_n: Optional[int] = None,
        min_dist: Optional[float] = None
    ) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame], Optional[pd.DataFrame]]:
        """
        Merge all unmixed and mixed spot tables from individual tiles.
        
        Args:
            round_n: Round number (uses Config.ROUND_N if not specified)
            min_dist: Minimum distance value (uses default if not specified)
            
        Returns:
            Tuple of (merged unmixed output, merged mixed output, merged mixed scratch)
        """
        round_n = round_n or self.config.ROUND_N
        min_dist = min_dist or self.config.min_dist
        
        # UNMIXED SPOTS
        # Build pattern for unmixed spots
        unmixed_pattern = f'unmixed_spots_R{round_n}'
   
        unmixed_full_pattern = f'{unmixed_pattern}*_tile_*'
        
        # Merge unmixed from output folder - search recursively for tile subfolders
        self.logger.info(f"\nMerging unmixed spot tables for min_dist={min_dist}...")
        self.logger.info(f"searching in {self.config.OUTPUT_FOLDER} for pattern '**/{unmixed_full_pattern}.pkl'")
        unmixed_output_files = []
        for file_path in self.config.OUTPUT_FOLDER.glob(f"**/{unmixed_full_pattern}.pkl"):
            unmixed_output_files.append(file_path)
        
        self.logger.info(f"Found {len(unmixed_output_files)} unmixed spot files to merge")
        
        unmixed_output_merged = None
        if unmixed_output_files:
            output_path = Path(self.config.OUTPUT_FOLDER).parent / f'unmixed_spots_R{round_n}_merged_minDist_{min_dist}.pkl'
            unmixed_output_merged = self.merge_pickle_tables(unmixed_output_files, output_path)
        

        # MIXED SPOTS
        # Build pattern for mixed spots
        mixed_pattern = f'mixed_spots_R{round_n}'
        mixed_full_pattern = f'{mixed_pattern}*_tile_*'
        
        # Merge mixed from output folder - search recursively for tile subfolders
        mixed_output_files = []
        for file_path in self.config.OUTPUT_FOLDER.glob(f"**/{mixed_full_pattern}.pkl"):
            mixed_output_files.append(file_path)
        
        self.logger.info(f"Found {len(mixed_output_files)} mixed spot files to merge")
        
        mixed_output_merged = None
        if mixed_output_files:
            output_path = Path(self.config.OUTPUT_FOLDER).parent / f'mixed_spots_R{round_n}_merged.pkl'
            mixed_output_merged = self.merge_pickle_tables(mixed_output_files, output_path)
        
        # Merge mixed from scratch folder
        mixed_scratch_files = []
        for file_path in self.config.SCRATCH_FOLDER.glob(f"**/{mixed_full_pattern}.pkl"):
            mixed_scratch_files.append(file_path)
        
        mixed_scratch_merged = None
        if mixed_scratch_files:
            scratch_path = Path(self.config.OUTPUT_FOLDER).parent / f'mixed_spots_R{round_n}_merged.pkl'
            mixed_scratch_merged = self.merge_pickle_tables(mixed_scratch_files, scratch_path)
        
        return unmixed_output_merged, mixed_output_merged, mixed_scratch_merged
    
    def merge_all_tables(
        self,
        min_dist: Optional[float] = None
    ) -> dict:
        """
        Merge all types of tables (unmixed and mixed spot tables).
        
        Args:
            min_dist: Minimum distance value for unmixed spot tables
            
        Returns:
            Dictionary with merged dataframes
        """
        results = {}
        
        self.logger.info("Starting table merge process...")
        
        # Merge spot tables (both unmixed and mixed)
        self.logger.info("\nMerging spot tables...")
        unmixed_output, mixed_output, mixed_scratch = self.merge_all_spot_tables(min_dist=min_dist)
        results['unmixed_spots_output'] = unmixed_output
        results['mixed_spots_output'] = mixed_output
        results['mixed_spots_scratch'] = mixed_scratch
        
        # Summary
        self.logger.info("\n" + "="*80)
        self.logger.info("MERGE SUMMARY")
        self.logger.info("="*80)
        
        if unmixed_output is not None:
            self.logger.info(f"Unmixed spots (output): {len(unmixed_output)} rows")
        else:
            self.logger.info("Unmixed spots (output): No data merged")
        
        if mixed_output is not None:
            self.logger.info(f"Mixed spots (output): {len(mixed_output)} rows")
        else:
            self.logger.info("Mixed spots (output): No data merged")
            
        self.logger.info("="*80 + "\n")
        
        return results
