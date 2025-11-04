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
        min_dist: Optional[int] = None
    ) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
        """
        Merge all unmixed spot tables from individual tiles.
        
        Args:
            round_n: Round number (uses Config.ROUND_N if not specified)
            min_dist: Minimum distance value (uses default if not specified)
            
        Returns:
            Tuple of (merged output dataframe, merged scratch dataframe)
        """
        round_n = round_n or self.config.ROUND_N
        min_dist = min_dist or 3
        
        # Build pattern for unmixed spots
        pattern = f'unmixed_spots_R{round_n}'
        if min_dist:
            # Find files with specific minDist
            full_pattern = f'{pattern}*_tile_*_minDist_{min_dist}'
        else:
            full_pattern = f'{pattern}*_tile_*'
        
        # Merge from output folder
        output_files = []
        for file_path in self.config.OUTPUT_FOLDER.glob(f"{full_pattern}.pkl"):
            output_files.append(file_path)
        
        output_merged = None
        if output_files:
            output_path = self.config.OUTPUT_FOLDER / f'unmixed_spots_R{round_n}_merged_minDist_{min_dist}.pkl'
            output_merged = self.merge_pickle_tables(output_files, output_path)
        
        # Merge from scratch folder
        scratch_files = []
        for file_path in self.config.SCRATCH_FOLDER.glob(f"{full_pattern}.pkl"):
            scratch_files.append(file_path)
        
        scratch_merged = None
        if scratch_files:
            scratch_path = self.config.SCRATCH_FOLDER / f'unmixed_spots_R{round_n}_merged_minDist_{min_dist}.pkl'
            scratch_merged = self.merge_pickle_tables(scratch_files, scratch_path)
        
        return output_merged, scratch_merged
    
    def merge_all_cell_by_gene_tables(self) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
        """
        Merge all cell-by-gene tables (both unmixed and mixed) from individual tiles.
        
        Returns:
            Tuple of (merged unmixed dataframe, merged mixed dataframe)
        """
        # Merge unmixed cell-by-gene tables (CSV)
        unmixed_csv_files = self.find_tile_tables(
            self.config.OUTPUT_FOLDER,
            'unmixed_cell_by_gene',
            '.csv'
        )
        unmixed_merged = None
        if unmixed_csv_files:
            output_path = self.config.OUTPUT_FOLDER / 'unmixed_cell_by_gene_merged.csv'
            unmixed_merged = self.merge_csv_tables(unmixed_csv_files, output_path)
            
            # Also save pickle version from scratch folder
            unmixed_pkl_files = self.find_tile_tables(
                self.config.SCRATCH_FOLDER,
                'unmixed_cell_by_gene',
                '.pkl'
            )
            if unmixed_pkl_files:
                scratch_path = self.config.SCRATCH_FOLDER / 'unmixed_cell_by_gene_merged.pkl'
                self.merge_pickle_tables(unmixed_pkl_files, scratch_path)
        
        # Merge mixed cell-by-gene tables (CSV)
        mixed_csv_files = self.find_tile_tables(
            self.config.OUTPUT_FOLDER,
            'mixed_cell_by_gene',
            '.csv'
        )
        mixed_merged = None
        if mixed_csv_files:
            output_path = self.config.OUTPUT_FOLDER / 'mixed_cell_by_gene_merged.csv'
            mixed_merged = self.merge_csv_tables(mixed_csv_files, output_path)
            
            # Also save pickle version from scratch folder
            mixed_pkl_files = self.find_tile_tables(
                self.config.SCRATCH_FOLDER,
                'mixed_cell_by_gene',
                '.pkl'
            )
            if mixed_pkl_files:
                scratch_path = self.config.SCRATCH_FOLDER / 'mixed_cell_by_gene_merged.pkl'
                self.merge_pickle_tables(mixed_pkl_files, scratch_path)
        
        return unmixed_merged, mixed_merged
    
    def merge_all_tables(
        self,
        min_dist: Optional[int] = None
    ) -> dict:
        """
        Merge all types of tables (spot tables and cell-by-gene tables).
        
        Args:
            min_dist: Minimum distance value for spot tables
            
        Returns:
            Dictionary with merged dataframes
        """
        results = {}
        
        self.logger.info("Starting table merge process...")
        
        # Merge spot tables
        self.logger.info("\nMerging unmixed spot tables...")
        spot_output, spot_scratch = self.merge_all_spot_tables(min_dist=min_dist)
        results['unmixed_spots'] = spot_output
        
        # Merge cell-by-gene tables
        self.logger.info("\nMerging cell-by-gene tables...")
        unmixed_cbg, mixed_cbg = self.merge_all_cell_by_gene_tables()
        results['unmixed_cell_by_gene'] = unmixed_cbg
        results['mixed_cell_by_gene'] = mixed_cbg
        
        # Summary
        self.logger.info("\n" + "="*80)
        self.logger.info("MERGE SUMMARY")
        self.logger.info("="*80)
        for key, df in results.items():
            if df is not None:
                self.logger.info(f"{key}: {len(df)} rows")
            else:
                self.logger.info(f"{key}: No data merged")
        self.logger.info("="*80 + "\n")
        
        return results
