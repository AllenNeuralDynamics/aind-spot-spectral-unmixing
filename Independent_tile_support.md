# Independent Tile Support Implementation Plan

## Overview

This document outlines the implementation plan for adding multi-tile processing support to the aind-spot-spectral-unmixing package. The goal is to enable the package to process multiple microscopy data volumes (tiles) sequentially, with each tile processed independently and results saved separately.

## Background

### Current State
The package currently processes a single data volume with stats channels organized as:
```
ch_488_stats/
  image_data_488_versus_spots_488.csv
  image_data_488_versus_spots_561.csv
  ...
ch_561_stats/
  ...
```

### Target State
The package needs to handle multiple tiles with the following structure:
```
Tile_X_0000_Y_0000_Z_0000_ch_488_spots/
Tile_X_0000_Y_0000_Z_0000_ch_488_stats/
  image_data_Tile_X_0000_Y_0000_Z_0000_ch_488_versus_spots_488.csv
  image_data_Tile_X_0000_Y_0000_Z_0000_ch_488_versus_spots_561.csv
  image_data_Tile_X_0000_Y_0000_Z_0000_ch_488_versus_spots_594.csv
  image_data_Tile_X_0000_Y_0000_Z_0000_ch_488_versus_spots_638.csv
Tile_X_0000_Y_0000_Z_0000_ch_561_spots/
Tile_X_0000_Y_0000_Z_0000_ch_561_stats/
  ...
Tile_X_0000_Y_0001_Z_0000_ch_488_spots/
Tile_X_0000_Y_0001_Z_0000_ch_488_stats/
  ...
Tile_X_0001_Y_0000_Z_0000_ch_488_spots/
Tile_X_0001_Y_0000_Z_0000_ch_488_stats/
  ...
```

## Implementation Plan

---

## Phase 1: Update Config to Extract Tile Information

**File:** `code/spot_analysis/config.py`

### 1.1 Add Tile Tracking Class Variable

Add a new class variable to track the current tile being processed:

```python
class Config:
    # ... existing variables ...
    
    CURRENT_TILE: Optional[str] = None  # Track which tile is being processed
    folder_paths = None  # Will now be tile-aware
```

### 1.2 Update `get_folder_paths_pipeline()` Method

**Current regex:**
```python
spot_regex = r".*(\d{1,3})_stats\/image_data_.*_(\d{1,3})_versus_spots_(\d{1,3})\.csv"
```

**New regex to capture tile names:**
```python
spot_regex = r".*(Tile_X_\d{4}_Y_\d{4}_Z_\d{4}).*_ch_(\d{1,3})_stats/image_data_.*_ch_(\d{1,3})_versus_spots_(\d{1,3})\.csv"
```

**New return structure:**
```python
{
    'Tile_X_0000_Y_0000_Z_0000': {
        'spots_folders': {
            '488': 'Tile_X_0000_Y_0000_Z_0000_ch_488_stats/...',
            '561': 'Tile_X_0000_Y_0000_Z_0000_ch_561_stats/...',
            ...
        },
        'multichan_folders': {
            '488': {
                '561': 'Tile_X_0000_Y_0000_Z_0000_ch_488_stats/...versus_spots_561.csv',
                '594': '...',
                ...
            },
            '561': {...},
            ...
        }
    },
    'Tile_X_0000_Y_0001_Z_0000': {
        'spots_folders': {...},
        'multichan_folders': {...}
    },
    ...
}
```

**Implementation:**

```python
@classmethod
def get_folder_paths_pipeline(cls) -> Dict[str, Dict[str, Dict[str, str]]]:
    """Returns folder paths organized by tile from what is attached in /data/"""
    # Updated regex to capture tile name
    spot_regex = r".*(Tile_X_\d{4}_Y_\d{4}_Z_\d{4}).*_ch_(\d{1,3})_stats/image_data_.*_ch_(\d{1,3})_versus_spots_(\d{1,3})\.csv"
    exclude = set(['*.zarr'])
    
    # Organize by tile first
    tile_data = {}

    for root, dirs, files in os.walk(cls.DATA_FOLDER):
        # Exclude .zarr directories
        dirs[:] = [d for d in dirs if not d.endswith('.zarr')]
        for file in files:
            # Skip files within .zarr directories
            if '.zarr' in root:
                continue
            full_path = os.path.join(root, file)
            relative_path = os.path.relpath(full_path, cls.DATA_FOLDER)

            # Check for spot intensity files
            spot_match = re.match(spot_regex, relative_path)
            if spot_match:
                tile_name = spot_match.group(1)
                source_channel = spot_match.group(2)
                target_channel = spot_match.group(3)

                # Initialize tile structure if not exists
                if tile_name not in tile_data:
                    tile_data[tile_name] = {
                        'spots_folders': {},
                        'multichan_folders': {}
                    }

                if source_channel == target_channel: 
                    tile_data[tile_name]['spots_folders'][source_channel] = relative_path
                else:
                    if source_channel not in tile_data[tile_name]['multichan_folders']:
                        tile_data[tile_name]['multichan_folders'][source_channel] = {}
                    tile_data[tile_name]['multichan_folders'][source_channel][target_channel] = relative_path
    
    return tile_data
```

### 1.3 Add New Helper Methods

```python
@classmethod
def get_unique_tiles(cls) -> List[str]:
    """Extract list of unique tile names from folder structure"""
    if cls.folder_paths is None:
        cls.folder_paths = cls.get_folder_paths_pipeline()
    return list(cls.folder_paths.keys())

@classmethod
def get_folder_paths_for_tile(cls, tile_name: str) -> Dict[str, Dict[str, str]]:
    """Get folder paths filtered for a specific tile"""
    if cls.folder_paths is None:
        cls.folder_paths = cls.get_folder_paths_pipeline()
    
    if tile_name not in cls.folder_paths:
        raise ValueError(f"Tile {tile_name} not found in folder paths. "
                        f"Available tiles: {list(cls.folder_paths.keys())}")
    
    return cls.folder_paths[tile_name]

@classmethod
def set_current_tile(cls, tile_name: str) -> None:
    """Set the active tile context"""
    if cls.folder_paths is None:
        cls.folder_paths = cls.get_folder_paths_pipeline()
    
    if tile_name not in cls.folder_paths:
        raise ValueError(f"Tile {tile_name} not found. "
                        f"Available tiles: {list(cls.folder_paths.keys())}")
    
    cls.CURRENT_TILE = tile_name
    print(f"Set current tile to: {tile_name}")
```

### 1.4 Update `get_folder_paths()` Method

```python
@classmethod
def get_folder_paths(cls) -> Dict[str, Dict[str, str]]:
    """Get folder paths for the current tile"""
    if cls.CURRENT_TILE is None:
        # If no tile is set, return old behavior for backward compatibility
        # or raise an error
        tiles = cls.get_unique_tiles()
        if len(tiles) == 1:
            cls.CURRENT_TILE = tiles[0]
            print(f"Auto-selected single tile: {cls.CURRENT_TILE}")
        else:
            raise ValueError(f"Multiple tiles found but no current tile set. "
                           f"Call set_current_tile() first. Available tiles: {tiles}")
    
    return cls.get_folder_paths_for_tile(cls.CURRENT_TILE)
```

### 1.5 Update `validate_folder_paths()` Method

Update to work with tile-aware structure:

```python
@classmethod
def validate_folder_paths(cls, folder_paths: Dict[str, Dict[str, str]], tile_name: str = None) -> None:
    """Validates the generated folder paths for a specific tile"""
    expected_channels = set(cls.get_round_channels().keys())
    
    tile_label = f" for tile {tile_name}" if tile_name else ""
    
    # Validate spots folders
    spots_channels = set(folder_paths['spots_folders'].keys())
    if spots_channels != expected_channels:
        missing = expected_channels - spots_channels
        extra = spots_channels - expected_channels
        print(f"Warning{tile_label}: Mismatch in spots folders. Missing: {missing}, Extra: {extra}")

    # Validate multichannel folders
    multichan_channels = set(folder_paths['multichan_folders'].keys())
    if multichan_channels != expected_channels:
        missing = expected_channels - multichan_channels
        extra = multichan_channels - expected_channels
        print(f"Warning{tile_label}: Mismatch in multichannel folders. Missing: {missing}, Extra: {extra}")

    for source_channel, targets in folder_paths['multichan_folders'].items():
        expected_targets = expected_channels - {source_channel}
        if set(targets.keys()) != expected_targets:
            missing = expected_targets - set(targets.keys())
            extra = set(targets.keys()) - expected_targets
            print(f"Warning{tile_label}: Mismatch in multichannel targets for channel {source_channel}. Missing: {missing}, Extra: {extra}")
```

### 1.6 Update `get_and_validate_folder_paths()` Method

```python
@classmethod
def get_and_validate_folder_paths(cls) -> Dict[str, Dict[str, str]]:
    """Gets folder paths and validates them for the current tile"""
    if cls.folder_paths is None: 
        cls.folder_paths = cls.get_folder_paths_pipeline()
        print(f'Found {len(cls.folder_paths)} tiles: {list(cls.folder_paths.keys())}')
    
    # Get paths for current tile
    tile_paths = cls.get_folder_paths()
    
    # Validate
    cls.validate_folder_paths(tile_paths, cls.CURRENT_TILE)
    
    return tile_paths
```

---

## Phase 2: Update Data Loader for Tile-Specific Paths

**File:** `code/spot_analysis/data_loader.py`

### 2.1 Add Tile Name to Loaded DataFrames

Update all methods that load data to add a `tile_name` column:

```python
def load_detected_spots_for_channel(self, ch):
    round_n = self.config.ROUND_N
    spots_folders = self.config.get_folder_paths()['spots_folders']
    spot_col_order = ['spot_id','chan','chan_spot_id','cell_id','round','z','y','x','z_center','y_center','x_center','dist','r'] 
    
    spots = pd.read_csv(self.config.SPOTS_FOLDER.joinpath(spots_folders[str(ch)]))
    spots = spots.rename(columns = {
        'Z':'z', 'Y': 'y', 'X': 'x', 
        'Z_center': 'z_center', 'Y_center': 'y_center', 'X_center': 'x_center', 
        'SEG_ID': 'cell_id', 'FG': f'chan_{ch}_fg', 'BG': f'chan_{ch}_bg'
    })
    
    spots['round'] = str(round_n)
    spots['chan'] = str(ch)
    spots['tile_name'] = self.config.CURRENT_TILE  # Add tile tracking
    spots['spot_id'] = range(1, len(spots)+1)
    spots['chan_spot_id'] = range(1, len(spots)+1)
    
    spots = spots[list(spot_col_order)+['tile_name']+list(np.setdiff1d(spots.columns, spot_col_order+['tile_name']))]
    return spots
```

Apply similar changes to:
- `load_multichannel_data()`
- `load_all_spots()`

### 2.2 Update `load_all_spots()` Method

Ensure it properly uses the current tile context:

```python
def load_all_spots(self) -> pd.DataFrame:
    """Load and combine all spot data for the current tile"""
    
    if self.config.CURRENT_TILE is None:
        raise ValueError("No current tile set. Call Config.set_current_tile() first.")
    
    channels = self.config.get_round_spot_channels()
    channels = [str(i) for i in list(channels) if i!= '405']

    round_n = self.config.ROUND_N
    
    # Get tile-specific folder paths
    multichan_folders = self.config.get_folder_paths()['multichan_folders']

    # ... rest of the method remains the same ...
    # but add tile_name to the final dataframe
    
    for ch in channels:
        mixed_spots_df['tile_name'] = self.config.CURRENT_TILE
    
    return mixed_spots_df
```

---

## Phase 3: Update Output Paths to Include Tile Names

### 3.1 Update `process_round.py`

**File:** `code/process_round.py`

#### Update `_setup_logging()` method:

```python
def _setup_logging(self) -> logging.Logger:
    """Setup logging configuration"""
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    
    # Include tile name in log filename
    tile_suffix = f"_tile_{Config.CURRENT_TILE}" if Config.CURRENT_TILE else ""
    log_filename = f'round_{Config.ROUND_N}{tile_suffix}_processing.log'
    
    # Create handlers
    console_handler = logging.StreamHandler()
    file_handler = logging.FileHandler(Config.OUTPUT_FOLDER / log_filename)
    
    # ... rest of setup ...
    
    return logger
```

#### Update `run()` method to include tile names in outputs:

```python
def run(self):
    """Run the complete spot analysis pipeline"""
    tile_suffix = f"_tile_{Config.CURRENT_TILE}" if Config.CURRENT_TILE else ""
    
    self.logger.info(f"Starting analysis for Round {Config.ROUND_N}, Tile {Config.CURRENT_TILE}")
    
    try:
        # ... load and process data ...
        
        # Save intermediate results with tile name
        spots_df.to_pickle(
            Config.OUTPUT_FOLDER / f'mixed_spots_R{Config.ROUND_N}{tile_suffix}.pkl'
        )
        spots_df.to_pickle(
            Config.SCRATCH_FOLDER / f'mixed_spots_R{Config.ROUND_N}{tile_suffix}.pkl'
        )
        
        # Calculate ratios with tile name
        ratio_path = Config.OUTPUT_FOLDER / f'r{Config.ROUND_N}{tile_suffix}_ratios.txt'
        
        # ... more processing ...
        
        # Save stats with tile name
        stats_df_csv_name = f'/results/spot_unmixing_stats{tile_suffix}.csv'
        stats_df.to_csv(stats_df_csv_name)
        
        # ... rest of processing ...
```

#### Update `_save_summary_statistics()` method:

```python
def _save_summary_statistics(self, results):
    """Save summary statistics for all processing runs"""
    summary_stats = []
    tile_suffix = f"_tile_{Config.CURRENT_TILE}" if Config.CURRENT_TILE else ""
    
    for min_dist, (unmixed_df, stats) in results.items():
        for channel_stat in stats:
            stat_dict = {
                'min_dist': min_dist,
                'round': Config.ROUND_N,
                'tile': Config.CURRENT_TILE,
                **channel_stat
            }
            summary_stats.append(stat_dict)
    
    # Create summary DataFrame and save
    summary_df = pd.DataFrame(summary_stats)
    summary_df.to_csv(
        Config.OUTPUT_FOLDER / f'round_{Config.ROUND_N}{tile_suffix}_summary_stats.csv',
        index=False
    )
```

### 3.2 Update `unmixer.py`

**File:** `code/spot_analysis/unmixer.py`

#### Update `_save_results()` method:

```python
def _save_results(self, unmixed_df: pd.DataFrame, min_dist: float) -> None:
    """Save unmixed spots to file"""
    tile_suffix = f"_tile_{self.config.CURRENT_TILE}" if self.config.CURRENT_TILE else ""
    
    output_path = (
        self.config.OUTPUT_FOLDER /
        f'unmixed_spots_R{self.config.ROUND_N}{tile_suffix}_minDist_{int(min_dist)}.pkl'
    )
    scratch_path = (
        self.config.SCRATCH_FOLDER /
        f'unmixed_spots_R{self.config.ROUND_N}{tile_suffix}_minDist_{int(min_dist)}.pkl'
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    unmixed_df.to_pickle(output_path)
    unmixed_df.to_pickle(scratch_path)
```

### 3.3 Update `cell_by_gene_table.py`

**File:** `code/spot_analysis/cell_by_gene_table.py`

#### Update `load_spots()` method to handle tile-aware file patterns:

```python
def load_spots(self, rounds: List[int], unmixed: bool = True) -> pd.DataFrame:
    """Load and process spots data for given rounds"""
    spots_df = pd.DataFrame()
    tile_suffix = f"_tile_{self.config.CURRENT_TILE}" if self.config.CURRENT_TILE else "*"
    
    for rn in rounds:
        round_chans = list(self.config.GENE_DICT[str(rn)].keys())
        
        file_suffix = 'unmixed_spots' if unmixed else 'mixed_spots'
        if unmixed:
            file_name = f'{file_suffix}_R{rn}{tile_suffix}_minDist_{self.config.min_dist}.pkl'
        else:
            file_name = f'{file_suffix}_R{rn}{tile_suffix}.pkl'
        
        # Use glob pattern if tile not specified
        if '*' in file_name:
            file_locations = list(self.config.SCRATCH_FOLDER.glob(file_name))
            if not file_locations:
                raise FileNotFoundError(f"No files found matching {file_name}")
            file_location = file_locations[0]
        else:
            file_location = self.config.SCRATCH_FOLDER / file_name
        
        with open(file_location, 'rb') as file:
            ch_spots_df = pickle.load(file)
            
        # ... rest of method ...
```

#### Update `process_pipeline()` method to save with tile names:

```python
def process_pipeline(self, rounds: List[int]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Run the complete processing pipeline"""
    tile_suffix = f"_tile_{self.config.CURRENT_TILE}" if self.config.CURRENT_TILE else ""
    
    # Process unmixed spots
    unmixed_spots = self.load_spots(rounds, unmixed=True)
    unmixed_spots_filtered = self.apply_spot_filters(unmixed_spots)
    segmentation = self.load_segmentation(rounds)
    unmixed_annotations = self.process_cell_annotations(unmixed_spots_filtered, segmentation)
    
    # Process mixed spots
    mixed_spots = self.load_spots(rounds, unmixed=False)
    mixed_annotations = self.process_cell_annotations(mixed_spots, segmentation)
    
    # Save results with tile suffix
    unmixed_annotations.to_pickle(
        self.config.SCRATCH_FOLDER / f'unmixed_cell_by_gene{tile_suffix}.pkl'
    )
    mixed_annotations.to_pickle(
        self.config.SCRATCH_FOLDER / f'mixed_cell_by_gene{tile_suffix}.pkl'
    )
    unmixed_annotations.to_csv(
        self.config.OUTPUT_FOLDER / f'unmixed_cell_by_gene{tile_suffix}.csv'
    )
    mixed_annotations.to_csv(
        self.config.OUTPUT_FOLDER / f'mixed_cell_by_gene{tile_suffix}.csv'
    )
    
    return unmixed_annotations, mixed_annotations
```

---

## Phase 4: Create Tile Processing Orchestrator

**File:** `code/run_capsule.py`

### 4.1 Add `TileProcessor` Class

Replace or update the existing `run()` function with a tile processing orchestrator:

```python
""" top level run script """

from process_round import SpotAnalysisPipeline
from spot_analysis.config import Config
from pathlib import Path
from typing import Dict, List, Any
import traceback


class TileProcessor:
    """Orchestrates processing of multiple tiles"""
    
    def __init__(self, min_distances: List[float] = None):
        """
        Initialize the tile processor
        
        Args:
            min_distances: List of minimum distances for unmixing (default: [3])
        """
        self.config = Config()
        self.min_distances = min_distances or [3]
        self.tiles = self.config.get_unique_tiles()
        
        print(f"\nFound {len(self.tiles)} tiles to process:")
        for tile in self.tiles:
            print(f"  - {tile}")
    
    def process_single_tile(self, tile_name: str) -> Dict[float, Any]:
        """
        Process a single tile
        
        Args:
            tile_name: Name of the tile to process
            
        Returns:
            Dictionary of results by minimum distance
        """
        print(f"\n{'='*80}")
        print(f"Processing tile: {tile_name}")
        print(f"{'='*80}\n")
        
        try:
            # Set current tile context
            self.config.set_current_tile(tile_name)
            
            # Create and run pipeline for this tile
            pipeline = SpotAnalysisPipeline(min_distances=self.min_distances)
            tile_results = pipeline.run()
            
            print(f"\n{'='*80}")
            print(f"Successfully completed processing for tile: {tile_name}")
            print(f"{'='*80}\n")
            
            return tile_results
            
        except Exception as e:
            print(f"\n{'='*80}")
            print(f"ERROR processing tile {tile_name}:")
            print(f"{'='*80}")
            print(traceback.format_exc())
            print(f"{'='*80}\n")
            return None
    
    def process_all_tiles(self) -> Dict[str, Dict[float, Any]]:
        """
        Process all tiles sequentially
        
        Returns:
            Dictionary mapping tile names to their results
        """
        all_results = {}
        successful_tiles = []
        failed_tiles = []
        
        for i, tile_name in enumerate(self.tiles, 1):
            print(f"\n\n{'#'*80}")
            print(f"# Processing tile {i}/{len(self.tiles)}: {tile_name}")
            print(f"{'#'*80}\n")
            
            tile_results = self.process_single_tile(tile_name)
            
            if tile_results is not None:
                all_results[tile_name] = tile_results
                successful_tiles.append(tile_name)
            else:
                failed_tiles.append(tile_name)
        
        # Print final summary
        self._print_final_summary(all_results, successful_tiles, failed_tiles)
        
        return all_results
    
    def _print_final_summary(
        self, 
        all_results: Dict[str, Dict[float, Any]], 
        successful_tiles: List[str],
        failed_tiles: List[str]
    ) -> None:
        """Print summary across all tiles"""
        print("\n\n")
        print("=" * 80)
        print("FINAL SUMMARY - ALL TILES")
        print("=" * 80)
        print(f"\nTotal tiles found: {len(self.tiles)}")
        print(f"Successfully processed: {len(successful_tiles)}")
        print(f"Failed: {len(failed_tiles)}")
        
        if failed_tiles:
            print(f"\nFailed tiles:")
            for tile in failed_tiles:
                print(f"  - {tile}")
        
        print("\n" + "-" * 80)
        print("Results by Tile:")
        print("-" * 80)
        
        for tile_name, tile_results in all_results.items():
            print(f"\nTile: {tile_name}")
            
            for min_dist, (unmixed_df, stats) in tile_results.items():
                print(f"  Minimum distance: {min_dist}")
                
                for channel_stat in stats:
                    kept_pct = (channel_stat['kept_spots'] / channel_stat['total_spots'] * 100 
                               if channel_stat['total_spots'] > 0 else 0)
                    
                    print(f"    Channel {channel_stat['channel']} ({channel_stat['gene']}): "
                          f"{channel_stat['kept_spots']}/{channel_stat['total_spots']} spots "
                          f"({kept_pct:.1f}%)")
        
        print("\n" + "=" * 80)
        print("Processing complete!")
        print("=" * 80 + "\n")


def run():
    """Main entry point for tile processing"""
    # TODO: Add argparse for configurable parameters
    # parser = argparse.ArgumentParser()
    # parser.add_argument("--min-distances", nargs="+", type=float, default=[3.0],
    #                     help="Minimum distances for unmixing")
    # args = parser.parse_args()
    
    # Create tile processor
    processor = TileProcessor(min_distances=[3])
    
    # Process all tiles
    all_results = processor.process_all_tiles()
    
    return all_results


if __name__ == "__main__":
    run()
```

---

## Phase 5: Testing Strategy

### 5.1 Test with Single Tile
First, verify backward compatibility by testing with a single tile dataset:
- Should auto-detect and process the single tile
- Output filenames should include tile name
- Results should match previous single-tile behavior

### 5.2 Test with Multiple Tiles
Test with 2-3 tiles to verify:
- Each tile processes independently
- No data mixing between tiles
- All tiles complete successfully
- Output files correctly separated by tile name

### 5.3 Verify Output Structure
Check that all output files include tile identification:
```
/results/
  mixed_spots_R0_tile_Tile_X_0000_Y_0000_Z_0000.pkl
  mixed_spots_R0_tile_Tile_X_0000_Y_0001_Z_0000.pkl
  unmixed_spots_R0_tile_Tile_X_0000_Y_0000_Z_0000_minDist_3.pkl
  unmixed_spots_R0_tile_Tile_X_0000_Y_0001_Z_0000_minDist_3.pkl
  unmixed_cell_by_gene_tile_Tile_X_0000_Y_0000_Z_0000.csv
  unmixed_cell_by_gene_tile_Tile_X_0000_Y_0001_Z_0000.csv
  round_0_tile_Tile_X_0000_Y_0000_Z_0000_processing.log
  round_0_tile_Tile_X_0000_Y_0001_Z_0000_processing.log
  ...
```

### 5.4 Error Handling
Test error scenarios:
- Missing files for one tile (should fail gracefully for that tile only)
- Malformed tile names
- Empty tile directories

---

## Implementation Checklist

### Phase 1: Config Updates
- [ ] Add `CURRENT_TILE` class variable
- [ ] Update `get_folder_paths_pipeline()` with new regex
- [ ] Implement `get_unique_tiles()` method
- [ ] Implement `get_folder_paths_for_tile()` method
- [ ] Implement `set_current_tile()` method
- [ ] Update `get_folder_paths()` method
- [ ] Update `validate_folder_paths()` method
- [ ] Update `get_and_validate_folder_paths()` method

### Phase 2: Data Loader Updates
- [ ] Add `tile_name` column in `load_detected_spots_for_channel()`
- [ ] Add `tile_name` column in `load_multichannel_data()`
- [ ] Update `load_all_spots()` to use current tile context
- [ ] Verify no tile data mixing

### Phase 3: Output Path Updates
- [ ] Update `_setup_logging()` in `process_round.py`
- [ ] Update all output paths in `run()` method
- [ ] Update `_save_summary_statistics()` method
- [ ] Update `_save_results()` in `unmixer.py`
- [ ] Update `load_spots()` in `cell_by_gene_table.py`
- [ ] Update `process_pipeline()` output paths

### Phase 4: Orchestrator
- [ ] Create `TileProcessor` class in `run_capsule.py`
- [ ] Implement `process_single_tile()` method
- [ ] Implement `process_all_tiles()` method
- [ ] Implement `_print_final_summary()` method
- [ ] Update `run()` function

### Phase 5: Testing
- [ ] Test single tile dataset
- [ ] Test multiple tile dataset
- [ ] Verify output file naming
- [ ] Test error handling
- [ ] Verify no data contamination between tiles

---

## Key Design Decisions

### 1. Tile Context Management
Using a class variable `CURRENT_TILE` in Config to track which tile is being processed. This allows all downstream code to access the tile context without passing it through every function call.

### 2. Backward Compatibility
If only one tile is found, the system auto-selects it, maintaining compatibility with single-tile workflows.

### 3. Independent Processing
Each tile is processed completely before moving to the next. This ensures:
- Clear error boundaries (one tile failing doesn't affect others)
- Simpler memory management
- Easier debugging and logging

### 4. File Naming Convention
All output files include the tile name as a suffix: `{base_name}_tile_{tile_name}.{ext}`

### 5. Error Handling
The `TileProcessor` catches exceptions per tile and continues processing remaining tiles, providing a summary of successes and failures at the end.

---

## Potential Enhancements (Future Work)

1. **Parallel Processing**: Process multiple tiles in parallel for better performance
2. **Tile Aggregation**: Combine results from all tiles into summary files
3. **Configuration File**: Support for tile-specific parameters via config file
4. **Resume Capability**: Skip already-processed tiles on re-run
5. **Tile Selection**: CLI arguments to process specific tiles only

---

## Notes

- The regex pattern assumes tile names follow the format: `Tile_X_####_Y_####_Z_####`
- All intermediate and final outputs will be tile-specific
- Logging is per-tile for easier debugging
- The manifest file (`processing_manifest.json`) is assumed to be tile-independent (same channels/genes for all tiles)

---

## Questions for Clarification

Before implementation, consider:
1. Should there be a combined summary file aggregating all tiles?
2. Should tiles be processed in a specific order (e.g., sorted by coordinates)?
3. What should happen if a tile partially fails (e.g., some channels missing)?
4. Should there be a progress checkpoint system for very large numbers of tiles?
