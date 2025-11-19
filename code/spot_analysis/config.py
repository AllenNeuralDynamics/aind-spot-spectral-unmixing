import pathlib
from typing import Dict, Any, List, Set, Optional
import os
import re
import json


class Config:

class Config:
    """Configuration for spot analysis pipeline.
    
    Each instance maintains isolated configuration state for independent processing.
    """
    
    # CONSTANTS: These are truly static defaults that should be shared (read-only)
    DEFAULT_GENE_DICT: Dict[str, Dict[str, str]] = {
        '0': {'1': 'Vip', '2': 'Sst', '4': 'Slc17a7'},
        '1': {'1': 'Cbln4', '2': 'Cdk18', '3': 'Kcnab1', '4': 'Nos1'},
        '2': {'1': 'Adcyap1', '2': 'Rorb', '3': 'Myh7', '4': 'Pdyn'},
        '3': {'1': 'Wfs1', '2': 'Npnt', '3': 'F2r12', '4': 'Trp53i11'},
        '4': {'1': 'Thsd7a', '2': 'Syt6', '3': 'Car4', '4': 'Tmem215'},
        '5': {'1': 'Pvalb', '2': 'Olig1', '3': 'Lypd1', '4': 'Synpr'},
        '6': {'1': 'Parm1', '2': 'Sfrp2', '3': 'Tnnc1', '4': 'Penk'},
        '7': {'1': 'Etv1', '2': 'Lsp1', '3': 'Slc18a3', '4': 'Calb1'},
        '8': {'1': 'Alcam', '2': 'Cidea', '3': 'Prss23', '4': 'Il1rap12'},
        '9': {'1': 'Cplx', '2': 'Ctss', '3': 'Npy'},
        '10': {'1': 'Slc18a8', '2': 'Tshz2', '3': 'Egln3', '4': 'Lpl'},
        '11': {'1': 'Gad2', '2': 'Ostn', '3': 'Lhx6', '4': 'Stk17b'},
        '12': {'1': 'Cck', '2': 'Crispld2', '3': 'Nmbr', '4': 'Anxa2'},
        '13': {'1': 'Snap25', '2': 'lgfbp4', '3': 'Chrm2', '4': 'Ndnf'}
    }
    
    # Default values for parameters (used as fallbacks during initialization)
    _DEFAULT_ROUND_N = 0
    _DEFAULT_MIN_DISTS = 5
    _DEFAULT_PERCENTILE = 95
    _DEFAULT_FRAC_SAMPLED = 0.1
    _DEFAULT_N_SUBSET = 100000
    _DEFAULT_EPOCHS = 10000
    _DEFAULT_RESAMPLE_ITER = 50
    _DEFAULT_L1 = 0
    _DEFAULT_LEARNING_RATE = 1e-9
    _DEFAULT_CENT_CUTOFF = 1
    _DEFAULT_CORR_CUTOFF = 0.5
    _DEFAULT_DIST_CUTOFF = 4
    _DEFAULT_MIN_DIST = 3
    _DEFAULT_VOLUME_QUANTILES = (0.08, 0.5, 0.95)
    _DEFAULT_OUTPUT_DATA_TYPE = 'zarr'
    _DEFAULT_APPLY_STITCHING_TRANSFORM = True

    def __init__(self):
        """Initialize a new Config instance with isolated state."""
        
        # Path configuration - INSTANCE ATTRIBUTES (mutable per tile)
        self.SPOTS_FOLDER = pathlib.Path('/data/')
        self.DATA_FOLDER = pathlib.Path('/data/')
        self.OUTPUT_FOLDER = pathlib.Path('/results/')
        self.SCRATCH_FOLDER = pathlib.Path('/scratch/')
        self.OUTPUT_DATA_TYPE = self._DEFAULT_OUTPUT_DATA_TYPE
        
        # Stitching configuration
        self.STITCHING_XML_PATH: Optional[str] = self.DATA_FOLDER.joinpath(
            'image_tile_alignment/combined_stitching_cam_alignment_all_channels.xml'
        ).as_posix()
        self.APPLY_STITCHING_TRANSFORM: bool = self._DEFAULT_APPLY_STITCHING_TRANSFORM
        
        # Processing parameters - INSTANCE ATTRIBUTES
        self.ROUND_N = self._DEFAULT_ROUND_N
        self.MIN_DISTS = self._DEFAULT_MIN_DISTS
        self.PERCENTILE = self._DEFAULT_PERCENTILE
        
        # Demixing parameters - INSTANCE ATTRIBUTES
        self.FRAC_SAMPLED = self._DEFAULT_FRAC_SAMPLED
        self.N_SUBSET = self._DEFAULT_N_SUBSET
        self.EPOCHS = self._DEFAULT_EPOCHS
        self.RESAMPLE_ITER = self._DEFAULT_RESAMPLE_ITER
        self.L1 = self._DEFAULT_L1
        self.LEARNING_RATE = self._DEFAULT_LEARNING_RATE
        
        # QC parameters - INSTANCE ATTRIBUTES
        self.CENT_CUTOFF = self._DEFAULT_CENT_CUTOFF
        self.CORR_CUTOFF = self._DEFAULT_CORR_CUTOFF
        self.DIST_CUTOFF = self._DEFAULT_DIST_CUTOFF
        
        # Cell by gene table parameters - INSTANCE ATTRIBUTES
        self.min_dist = self._DEFAULT_MIN_DIST
        self.volume_quantiles = self._DEFAULT_VOLUME_QUANTILES
        
        # Tile-specific state - INSTANCE ATTRIBUTES
        self.CURRENT_TILE: Optional[str] = None
        self.folder_paths = None
        self.folder_paths_by_tile: Optional[Dict[str, Dict[str, Any]]] = None
        
        # Initialize from manifest
        self._load_manifest()
        self._update_round_from_manifest()
        self._make_gene_dict_from_manifest()
        self._ensure_folder_paths_loaded()

    #@classmethod
    #def _load_manifest(cls):
    #    """Load the processing manifest JSON file"""

        # manifest_path = pathlib.Path(cls.dataset_name) / 'derived' / 'processing_manifest.json'
    #    manifest_path = list(pathlib.Path(cls.DATA_FOLDER).glob("derived/processing_manifest.json"))
        
    
    #    if not len(manifest_path):
    #        print(f'didnt find pipeline processing manifest')
            #raise FileNotFoundError("No processing_manifest.json was found!")
        
    #        manifest_path = list(pathlib.Path(cls.DATA_FOLDER).glob("*/derived/processing_manifest.json"))
    #        if not len(manifest_path):
    #            raise FileNotFoundError("No capsule processing_manifest.json was found!")

        
    #    print(f'Manifest_path {manifest_path}')

    #    try:
    #        with open(manifest_path[0], 'r') as f:
    #            cls.manifest = json.load(f)
    #    except FileNotFoundError:
    #        cls.manifest = None
    #        raise FileNotFoundError(f"Processing manifest not found at {manifest_path}")

    @classmethod
    def _update_round_from_manifest(cls):
        if not cls.manifest:
            cls.ROUND_N = cls._default_ROUND_N
            return
        round = cls.manifest['round']
        # if round != -1: 
        cls.ROUND_N = round
        # else:
            # cls.ROUND_N = cls._default_ROUND_N

        """ Processing Manifest Json Example
    {'segmentation_channels': {'background': '405', 'nuclear': None}, 'spot_channels': ['561', '488', '638'], 'round': 1, 'stitching_channels': ['561', '488', '638'], 'gene_dict': {'405': {'gene': 'Rn28s', 'barcode': '', 'fluorophore': '', 'wavelength': 'dtype:', 'round': 1}, '561': {'gene': 'Calb2', 'barcode': 'B7', 'fluorophore': '', 'wavelength': '561,', 'round': 1}, '488': {'gene': 'Npy', 'barcode': 'B1', 'fluorophore': '', 'wavelength': '488,', 'round': 1}, '638': {'gene': 'Tac1', 'barcode': 'B3', 'fluorophore': '', 'wavelength': '638,', 'round': 1}}}"""

    @classmethod
    def _make_gene_dict_from_manifest(cls):
        """Make a gene_dict from the processing manifest"""
        if not cls.manifest:
            cls.GENE_DICT = cls.DEFAULT_GENE_DICT
            return
        spot_channels = cls.manifest['spot_channels']
        round = cls.manifest['round']
        manifest_gene_dict = cls.manifest['gene_dict']

        #gene_dict is a dict of dicts with keys: round { channel: gene_name}
        temp_dict = {}
        
        for channel, gene in manifest_gene_dict.items():
            temp_dict[str(channel)] = str(gene['gene'])
        
        gene_dict= {}
        gene_dict[str(round)] = temp_dict
        cls.GENE_DICT = gene_dict
        

    @classmethod
    def get_round_channels(cls) -> Dict[str, str]:
        return cls.GENE_DICT[str(cls.ROUND_N)]
        
        
    @classmethod
    def get_round_spot_channels(cls) -> List[str]:
        if cls.manifest and cls.manifest.get('spot_channels'):
            return cls.manifest['spot_channels']
        return list(cls.get_round_channels().keys())

    @classmethod
    def _ensure_folder_paths_loaded(cls) -> None:
        if cls.folder_paths_by_tile is not None:
            return

        tile_data = cls.get_folder_paths_pipeline()
        if not tile_data:
            raise FileNotFoundError("No spot intensity files were found under the data folder.")

        cls.folder_paths_by_tile = tile_data
        tile_list = list(tile_data.keys())
        print(f"Found {len(tile_list)} tile(s): {tile_list}")

        if cls.CURRENT_TILE is None and len(tile_list) == 1:
            cls.CURRENT_TILE = tile_list[0]
            print(f"Auto-selected single tile: {cls.CURRENT_TILE}")

        if cls.CURRENT_TILE is not None:
            cls.folder_paths = cls.folder_paths_by_tile.get(cls.CURRENT_TILE)

    @classmethod
    def get_unique_tiles(cls) -> List[str]:
        cls._ensure_folder_paths_loaded()
        if cls.folder_paths_by_tile is None:
            return []
        return list(cls.folder_paths_by_tile.keys())

    @classmethod
    def get_folder_paths_for_tile(cls, tile_name: str) -> Dict[str, Any]:
        cls._ensure_folder_paths_loaded()
        available_tiles = list(cls.folder_paths_by_tile.keys()) if cls.folder_paths_by_tile else []
        if cls.folder_paths_by_tile is None or tile_name not in cls.folder_paths_by_tile:
            raise ValueError(
                f"Tile {tile_name} not found in folder paths. Available tiles: {available_tiles}"
            )
        return cls.folder_paths_by_tile[tile_name]

    @classmethod
    def set_current_tile(cls, tile_name: str) -> None:
        cls._ensure_folder_paths_loaded()
        available_tiles = list(cls.folder_paths_by_tile.keys()) if cls.folder_paths_by_tile else []
        if cls.folder_paths_by_tile is None or tile_name not in cls.folder_paths_by_tile:
            raise ValueError(
                f"Tile {tile_name} not found. Available tiles: {available_tiles}"
            )
        cls.CURRENT_TILE = tile_name
        cls.folder_paths = cls.folder_paths_by_tile[tile_name]
        print(f"Set current tile to: {tile_name}")

    @classmethod
    def get_folder_paths(cls) -> Dict[str, Any]:
        """
        Get folder paths for the current tile only.
        Ensures that only data from CURRENT_TILE is returned.
        """
        cls._ensure_folder_paths_loaded()
        if cls.CURRENT_TILE is None:
            tiles = cls.get_unique_tiles()
            if len(tiles) == 1:
                cls.CURRENT_TILE = tiles[0]
                if cls.folder_paths_by_tile is not None:
                    cls.folder_paths = cls.folder_paths_by_tile[cls.CURRENT_TILE]
                print(f"Auto-selected single tile: {cls.CURRENT_TILE}")
            else:
                raise ValueError(
                    "Multiple tiles detected but no current tile set. "
                    "Call Config.set_current_tile(tile_name) before accessing folder paths."
                )
        if cls.folder_paths_by_tile is None:
            raise ValueError("Folder paths not initialized")
        
        # Validate that we're only returning paths for the current tile
        current_tile_paths = cls.folder_paths_by_tile[cls.CURRENT_TILE]
        
        # Double-check that all paths contain the current tile name to prevent cross-contamination
        tile_name = cls.CURRENT_TILE
        for channel, path in current_tile_paths.get('spots_folders', {}).items():
            if tile_name not in str(path) and 'single_tile' not in str(path):
                print(f"Warning: spots_folders path for channel {channel} does not contain tile name {tile_name}: {path}")
        
        for source_ch, targets in current_tile_paths.get('multichan_folders', {}).items():
            for target_ch, path in targets.items():
                if tile_name not in str(path) and 'single_tile' not in str(path):
                    print(f"Warning: multichan_folders path for {source_ch}->{target_ch} does not contain tile name {tile_name}: {path}")
        
        return current_tile_paths

    
    @classmethod
    def _load_manifest(cls):
        """Load the processing manifest JSON file"""
        #pipeline path
        manifest_path = list(pathlib.Path(cls.DATA_FOLDER).glob("processing_manifest.json"))
        # manifest_path = list(pathlib.Path('/data').glob("processing_manifest.json"))
        
        if not len(manifest_path):
            print('Didn\'t find pipeline processing manifest')
            manifest_path = list(pathlib.Path(cls.DATA_FOLDER).glob("*/processing_manifest.json"))
            if not len(manifest_path):
                raise FileNotFoundError("No capsule processing_manifest.json was found!")

        print(f'Manifest_path {manifest_path}')

        try:
            with open(manifest_path[0], 'r') as f:
                cls.manifest = json.load(f)
                print(f"Loaded manifest with channels: {cls.manifest.get('spot_channels', [])}")
        except FileNotFoundError:
            cls.manifest = None
            raise FileNotFoundError(f"Processing manifest not found at {manifest_path}")

    @classmethod
    def get_folder_paths_pipeline(cls) -> Dict[str, Dict[str, Any]]: #get_folder_paths_pipeline
        """Returns folder paths from what is attached in /data/, organized by tile"""
        tile_regex = re.compile(
            r"(Tile_X_\d{4}_Y_\d{4}_Z_\d{4})_ch_(\d{1,3})_stats/"
            r"image_data_.*_ch_(\d{1,3})_versus_spots_(\d{1,3})\.csv"
        )
        legacy_regex = re.compile(
            r".*(\d{1,3})_stats/image_data_.*_(\d{1,3})_versus_spots_(\d{1,3})\.csv"
        )

        tile_data: Dict[str, Dict[str, Any]] = {}
        legacy_spots: Dict[str, str] = {}
        legacy_multichan: Dict[str, Dict[str, str]] = {}

        for root, dirs, files in os.walk(cls.DATA_FOLDER):
            # Exclude .zarr directories
            dirs[:] = [d for d in dirs if not d.endswith('.zarr')]
            for file in files:
                # Skip files within .zarr directories
                if '.zarr' in root:
                    continue
                full_path = os.path.join(root, file)
                relative_path = os.path.relpath(full_path, cls.DATA_FOLDER)

                # Check for tile-based spot intensity files
                tile_match = tile_regex.match(relative_path)
                if tile_match:
                    tile_name = tile_match.group(1)
                    source_channel = tile_match.group(2)
                    target_channel = tile_match.group(4)

                    if tile_name not in tile_data:
                        tile_data[tile_name] = {
                            'spots_folders': {},
                            'multichan_folders': {}
                        }

                    if source_channel == target_channel:
                        tile_data[tile_name]['spots_folders'][source_channel] = relative_path
                    else:
                        multichan = tile_data[tile_name]['multichan_folders']
                        if source_channel not in multichan:
                            multichan[source_channel] = {}
                        multichan[source_channel][target_channel] = relative_path
                    continue

                legacy_match = legacy_regex.match(relative_path)
                if legacy_match:
                    source_channel = legacy_match.group(2)
                    target_channel = legacy_match.group(3)

                    if source_channel == target_channel:
                        legacy_spots[source_channel] = relative_path
                    else:
                        if source_channel not in legacy_multichan:
                            legacy_multichan[source_channel] = {}
                        legacy_multichan[source_channel][target_channel] = relative_path

        if tile_data:
            return tile_data

        if legacy_spots or legacy_multichan:
            return {
                'single_tile': {
                    'spots_folders': legacy_spots,
                    'multichan_folders': legacy_multichan
                }
            }

        return {}

    @classmethod
    def validate_folder_paths(cls, folder_paths: Dict[str, Any], tile_name: Optional[str] = None) -> None:
        """Validates the generated folder paths"""
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
            if not isinstance(targets, dict):
                print(
                    f"Warning{tile_label}: Expected multichannel targets for channel {source_channel} to be a dict, "
                    f"but found {type(targets)}."
                )
                continue

            expected_targets = expected_channels - {source_channel}
            target_keys = set(targets.keys())
            if target_keys != expected_targets:
                missing = expected_targets - target_keys
                extra = target_keys - expected_targets
                print(
                    f"Warning{tile_label}: Mismatch in multichannel targets for channel {source_channel}. "
                    f"Missing: {missing}, Extra: {extra}"
                )

    @classmethod
    def get_and_validate_folder_paths(cls) -> Dict[str, Any]:
        """Gets folder paths and validates them"""
        cls._ensure_folder_paths_loaded()
        tile_paths = cls.get_folder_paths()
        cls.validate_folder_paths(tile_paths, cls.CURRENT_TILE)
        cls.folder_paths = tile_paths
        return tile_paths
        
