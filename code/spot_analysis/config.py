import pathlib
from typing import Dict, Any, List, Set
import os
import re
import json


class Config:

    SPOTS_FOLDER = pathlib.Path('/data/')
    DATA_FOLDER = pathlib.Path('/data/')
    OUTPUT_FOLDER = pathlib.Path('/results/')
    OUTPUT_DATA_TYPE = 'zarr'
    SCRATCH_FOLDER = pathlib.Path('/scratch/')
    
    # Processing parameters
    _default_ROUND_N = 0
    MIN_DISTS = 5
    PERCENTILE = 90
    
    # Demixing parameters
    FRAC_SAMPLED = 0.2
    N_SUBSET = 80000
    EPOCHS = 10000
    RESAMPLE_ITER = 50
    L1 = 0
    LEARNING_RATE = 1e-9
    
    # QC parameters --- these are getting moved to qc capsule
    CENT_CUTOFF = 1
    CORR_CUTOFF = 0.5
    DIST_CUTOFF = 4
    # cell by gene table parameters    
    min_dist = 3
    volume_quantiles = (0.08, 0.5, 0.95)


    folder_paths = None

    # Gene dictionary
    DEFAULT_GENE_DICT: Dict[str, Dict[str, str]] = {'0':{'1': 'Vip', '2': 'Sst', '4': 'Slc17a7'},
            '1':{'1': 'Cbln4', '2': 'Cdk18', '3': 'Kcnab1', '4': 'Nos1'},
            '2':{'1': 'Adcyap1', '2': 'Rorb', '3': 'Myh7', '4': 'Pdyn'},
            '3':{'1': 'Wfs1', '2': 'Npnt', '3': 'F2r12', '4': 'Trp53i11'},
            '4':{'1': 'Thsd7a', '2': 'Syt6', '3': 'Car4', '4': 'Tmem215'},
            '5':{'1': 'Pvalb', '2': 'Olig1', '3': 'Lypd1', '4': 'Synpr'},
            '6':{'1': 'Parm1', '2': 'Sfrp2', '3': 'Tnnc1', '4': 'Penk'},
            '7':{'1': 'Etv1', '2': 'Lsp1', '3': 'Slc18a3', '4': 'Calb1'},
            '8':{'1': 'Alcam', '2': 'Cidea', '3': 'Prss23', '4': 'Il1rap12'},
            '9':{'1': 'Cplx', '2': 'Ctss', '3': 'Npy',},
            '10':{'1': 'Slc18a8', '2': 'Tshz2', '3': 'Egln3', '4': 'Lpl'},
            '11':{'1': 'Gad2', '2': 'Ostn', '3': 'Lhx6', '4': 'Stk17b'},
            '12':{'1': 'Cck', '2': 'Crispld2', '3': 'Nmbr', '4': 'Anxa2'},
            '13':{'1': 'Snap25', '2': 'lgfbp4', '3': 'Chrm2', '4': 'Ndnf'}}

    def __init__(self):

        
        self._load_manifest()
        self._update_round_from_manifest()
        self._make_gene_dict_from_manifest()
        self.folder_paths = None
        self.folder_paths = self.get_and_validate_folder_paths()

        

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
    def get_round_spot_channels(cls) -> Dict[str, str]:
        spot_channels = cls.manifest['spot_channels']
        return spot_channels
    
    @classmethod
    def get_folder_paths(cls) -> Dict[str, Dict[str, str]]:
        return cls.get_and_validate_folder_paths()

    
    @classmethod
    def _load_manifest(cls):
        """Load the processing manifest JSON file"""
        manifest_path = list(pathlib.Path(cls.DATA_FOLDER).glob("derived/processing_manifest.json"))
        
        if not len(manifest_path):
            print('Didn\'t find pipeline processing manifest')
            manifest_path = list(pathlib.Path(cls.DATA_FOLDER).glob("*/derived/processing_manifest.json"))
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
    def _find_stats_files(cls, path: pathlib.Path) -> Dict[str, List[Dict[str, str]]]:
        """Find all stats files in the given path"""
        stats_files = {}
        expected_channels = set(str(ch) for ch in cls.manifest.get('spot_channels', []))
        print(f"\nDebug: Walking directory {path} looking for stats files")
        
        for root, dirs, files in os.walk(path):
            if '.zarr' in root:
                continue
                
            print(f"\nDebug: Examining directory: {root}")
            print(f"Debug: Found directories: {dirs}")
            print(f"Debug: Found files: {files}")
                
            # Look for CSV files
            for file in files:
                if file.endswith('.csv'):
                    # Try to extract channel information from the file path and name
                    file_path = os.path.relpath(os.path.join(root, file), path)
                    
                    # Extract channel numbers from the CSV filename
                    csv_match = re.search(r'channel_(\d+)_versus_spots_(\d+)\.csv', file)
                    if csv_match:
                        source_channel = csv_match.group(1)
                        target_channel = csv_match.group(2)
                        
                        if source_channel in expected_channels:
                            if source_channel not in stats_files:
                                stats_files[source_channel] = []
                            
                            stats_files[source_channel].append({
                                'path': file_path,
                                'source_wavelength': source_channel,
                                'target_wavelength': target_channel
                            })
        
        print(f"Found stats files: {stats_files}")
        return stats_files

    @classmethod
    def _find_spots_files(cls, path: pathlib.Path) -> Dict[str, str]:
        """Find all spots files (spots.csv or spots.npy) in *_spots folders"""
        spots_files = {}
        expected_channels = set(str(ch) for ch in cls.manifest.get('spot_channels', []))

        # Pattern for spots folders (both traditional and tile format)
        folder_patterns = [
            r'.*(\d{1,3})_spots$',  # Matches both ch_ and channel_
        ]
        
        for root, dirs, files in os.walk(path):
            if '.zarr' in root:
                continue
                
            current_dir = os.path.basename(root)
            
            # Check if this is a spots directory
            for pattern in folder_patterns:
                match = re.match(pattern, current_dir)
                if match:
                    channel = match.group(1)
                    if channel in expected_channels:
                        # Look for spots file
                        for file in files:
                            if file in ['spots.csv']:
                                spots_files[channel] = os.path.relpath(
                                    os.path.join(root, file), path)
                                break
        
        print(f"\nDebug: Found spots files: {spots_files}")
        return spots_files

    @classmethod
    def get_folder_paths_pipeline(cls) -> Dict[str, Dict[str, str]]:
        """Returns folder paths from what is attached in /data/"""
        # Ensure manifest is loaded
        if cls.manifest is None:
            cls._load_manifest()
        
        # Find spots files
        spots_folders = cls._find_spots_files(cls.DATA_FOLDER)
        
        # Find stats files
        stats_files = cls._find_stats_files(cls.DATA_FOLDER)
        
        # Process stats files into multichannel format
        multichan_folders = {}
        for source_channel, file_list in stats_files.items():
            if source_channel not in multichan_folders:
                multichan_folders[source_channel] = {}
                
            for file_info in file_list:
                if file_info['source_wavelength'] == source_channel:
                    target_channel = file_info['target_wavelength']
                    if target_channel != source_channel:
                        multichan_folders[source_channel][target_channel] = file_info['path']
        
        print(f"Final folder paths:")
        print(f"Spots folders: {spots_folders}")
        print(f"Multichannel folders: {multichan_folders}")
        
        return {
            'spots_folders': spots_folders,
            'multichan_folders': multichan_folders
        }

    @classmethod
    def validate_folder_paths(cls, folder_paths: Dict[str, Dict[str, str]]) -> None:
        """Validates the generated folder paths"""
        if cls.manifest is None:
            cls._load_manifest()
            
        expected_channels = set(str(ch) for ch in cls.manifest.get('spot_channels', []))
        print(f"Expected channels from manifest: {expected_channels}")
        
        # Validate spots folders
        spots_channels = set(folder_paths['spots_folders'].keys())
        if spots_channels != expected_channels:
            missing = expected_channels - spots_channels
            extra = spots_channels - expected_channels
            print(f"Warning: Mismatch in spots folders. Missing: {missing}, Extra: {extra}")

        # Validate multichannel folders
        multichan_channels = set(folder_paths['multichan_folders'].keys())
        if multichan_channels != expected_channels:
            missing = expected_channels - multichan_channels
            extra = multichan_channels - expected_channels
            print(f"Warning: Mismatch in multichannel folders. Missing: {missing}, Extra: {extra}")

        for source_channel, targets in folder_paths['multichan_folders'].items():
            expected_targets = expected_channels - {source_channel}
            if set(targets.keys()) != expected_targets:
                missing = expected_targets - set(targets.keys())
                extra = set(targets.keys()) - expected_targets
                print(f"Warning: Mismatch in multichannel targets for channel {source_channel}. Missing: {missing}, Extra: {extra}")

    @classmethod
    def get_and_validate_folder_paths(cls) -> Dict[str, Dict[str, str]]:
        """Gets folder paths and validates them"""
        if cls.folder_paths is None: 
            folder_paths = cls.get_folder_paths_pipeline()
            cls.validate_folder_paths(folder_paths)
            cls.folder_paths = folder_paths
        return cls.folder_paths