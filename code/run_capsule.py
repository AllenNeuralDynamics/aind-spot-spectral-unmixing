""" top level run script """

from pathlib import Path
from typing import Any, Dict, List, Optional
import traceback

from process_round import SpotAnalysisPipeline
from spot_analysis.config import Config
from spot_analysis.table_merger import TableMerger


class TileProcessor:
    """Orchestrates sequential processing of tiles."""

    def __init__(
        self,
        spots_folder: Path = Path('/data/'),
        # spots_folder: Path = Path('/root/capsule/data/test_unmixing_independent_tiles'), # test folder /root/capsule/data/test_unmixing_independent_tiles
        output_folder: Path = Path('/results/'),
        min_distances: Optional[List[float]] = None
    ) -> None:
        self.spots_folder = spots_folder
        self.output_folder = output_folder
        self.min_distances = min_distances or [3.0]

        # Prime configuration and discover available tiles
        Config.SPOTS_FOLDER = spots_folder
        Config.DATA_FOLDER = Path('/data/')
        Config.OUTPUT_FOLDER = output_folder
        self.config = Config()
        self.tiles = self.config.get_unique_tiles()

        print(f"\nFound {len(self.tiles)} tile(s) to process:")
        for tile in self.tiles:
            print(f"  - {tile}")

    def process_single_tile(self, tile_name: str) -> Optional[Dict[float, Any]]:
        """Process a single tile and return results."""
        print(f"\n{'=' * 80}")
        print(f"Processing tile: {tile_name}")
        print(f"{'=' * 80}\n")

        try:
            # Configure tile context
            self.config.set_current_tile(tile_name)
            tile_spots_folder = self.spots_folder 

            # Run pipeline for tile
            pipeline = SpotAnalysisPipeline(
                spots_folder=tile_spots_folder,
                output_folder=self.output_folder.joinpath(tile_name),
                min_distances=self.min_distances, 
                Config = self.config
            )
            tile_results = pipeline.run()

            print(f"\n{'=' * 80}")
            print(f"Completed tile: {tile_name}")
            print(f"{'=' * 80}\n")
            return tile_results

        except Exception:
            print(f"\n{'=' * 80}")
            print(f"ERROR processing tile {tile_name}:")
            print(f"{'=' * 80}")
            print(traceback.format_exc())
            print(f"{'=' * 80}\n")
            return None

    def process_all_tiles(self) -> Dict[str, Dict[float, Any]]:
        """Process all discovered tiles sequentially."""
        all_results: Dict[str, Dict[float, Any]] = {}
        successful_tiles: List[str] = []
        failed_tiles: List[str] = []

        total_tiles = len(self.tiles)
        for index, tile_name in enumerate(self.tiles, start=1):
            print(f"\n{'#' * 80}")
            print(f"# Tile {index}/{total_tiles}: {tile_name}")
            print(f"{'#' * 80}\n")

            tile_results = self.process_single_tile(tile_name)
            if tile_results is None:
                failed_tiles.append(tile_name)
                continue

            all_results[tile_name] = tile_results
            successful_tiles.append(tile_name)

            self._print_tile_summary(tile_name, tile_results)

        self._print_final_summary(successful_tiles, failed_tiles, all_results)
        
        # Merge all tile tables into combined tables
        if len(successful_tiles) > 1:
            self._merge_tile_tables()
        
        return all_results
    
    def _merge_tile_tables(self) -> None:
        """Merge all single-tile tables into combined multi-tile tables."""
        print("\n" + "=" * 80)
        print("MERGING TILE TABLES")
        print("=" * 80)
        
        merger = TableMerger(self.config)
        
        # Merge all table types
        merged_results = merger.merge_all_tables(min_dist=self.min_distances[0])
        
        print("\nTable merge complete.")
        print("=" * 80 + "\n")

    def _print_tile_summary(self, tile_name: str, tile_results: Dict[float, Any]) -> None:
        """Print per-tile summary statistics."""
        print(f"\nTile summary: {tile_name}")
        for min_dist, (_, stats) in tile_results.items():
            print(f"  Minimum distance: {min_dist}")
            for channel_stat in stats:
                kept = channel_stat.get('kept_spots', 0)
                total = channel_stat.get('total_spots', 0)
                kept_pct = (kept / total * 100) if total else 0.0
                print(
                    f"    Channel {channel_stat.get('channel')} ({channel_stat.get('gene')}): "
                    f"{kept}/{total} spots ({kept_pct:.1f}%)"
                )

    def _print_final_summary(
        self,
        successful_tiles: List[str],
        failed_tiles: List[str],
        all_results: Dict[str, Dict[float, Any]]
    ) -> None:
        """Print summary across all processed tiles."""
        print("\n" + "=" * 80)
        print("FINAL SUMMARY")
        print("=" * 80)
        print(f"Total tiles discovered: {len(self.tiles)}")
        print(f"Successfully processed: {len(successful_tiles)}")
        print(f"Failed: {len(failed_tiles)}")

        if failed_tiles:
            print("\nFailed tiles:")
            for tile in failed_tiles:
                print(f"  - {tile}")

        if successful_tiles:
            print("\nSummary by tile:")
            for tile in successful_tiles:
                tile_results = all_results[tile]
                for min_dist, (_, stats) in tile_results.items():
                    kept_counts = [s.get('kept_spots', 0) for s in stats]
                    total_counts = [s.get('total_spots', 0) for s in stats]
                    total_kept = sum(kept_counts)
                    total_spots = sum(total_counts)
                    kept_pct = (total_kept / total_spots * 100) if total_spots else 0.0
                    print(
                        f"  - {tile} | min_dist={min_dist}: "
                        f"kept {total_kept}/{total_spots} spots ({kept_pct:.1f}%)"
                    )

        print("\nProcessing complete.\n")


def run() -> Dict[str, Dict[float, Any]]:
    """Entry point for processing all available tiles."""
    processor = TileProcessor(min_distances=[3.0])
    return processor.process_all_tiles()


if __name__ == "__main__":
    run()
