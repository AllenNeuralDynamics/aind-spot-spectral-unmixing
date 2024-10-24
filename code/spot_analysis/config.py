import pathlib
from typing import Dict, Any
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
    
    # QC parameters
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

        

    @classmethod
    def _load_manifest(cls):
        """Load the processing manifest JSON file"""

        # manifest_path = pathlib.Path(cls.dataset_name) / 'derived' / 'processing_manifest.json'
        manifest_path = list(pathlib.Path(cls.DATA_FOLDER).glob("*/derived/processing_manifest.json"))
        
    
        if not len(manifest_path):
            raise FileNotFoundError("No processing_manifest.json was found!")

        try:
            with open(manifest_path[0], 'r') as f:
                cls.manifest = json.load(f)
        except FileNotFoundError:
            cls.manifest = None
            raise FileNotFoundError(f"Processing manifest not found at {manifest_path}")

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
    def get_folder_paths_manual(cls) -> Dict[str, Dict[str, str]]:
        """Returns folder paths based on round number"""
        if cls.ROUND_N == 0:
            return {
                'spots_folders': {
                    '1': 'HCR_BL6-000-R0_ch1_spot_intensity/spots.npy',
                    '2': 'HCR_BL6-000-R0_ch2_spot_intensity/spots.npy',
                    '4': 'HCR_BL6-000-R0_ch4_spot_intensity/spots.npy'
                },
                'multichan_folders' : {'1':{'2':'HCR_BL6-000-R0_ch1_multichannel/ch_channel_4_spots_channel_2.npy',
                              '4':'HCR_BL6-000-R0_ch1_multichannel/ch_channel_4_spots_channel_4.npy'},
                         '2':{'1':'HCR_BL6-000-R0_ch2_multichannel/ch_channel_4_spots_channel_1.npy',
                              '4':'HCR_BL6-000-R0_ch2_multichannel/ch_channel_4_spots_channel_4.npy'},
                         '4':{'1':'HCR_BL6-000-R0_ch4_multichannel/ch_channel_4_spots_channel_1.npy',
                              '2':'HCR_BL6-000-R0_ch4_multichannel/ch_channel_4_spots_channel_2.npy'}}
            }

        if cls.ROUND_N == 13:
            return {
                    'spots_folders':{'1':'HCR_BL6-000-R13_ch1_spot_intensity/spots.npy',
                                        '2':'HCR_BL6-000-R13_ch2_spot_intensity/spots.npy',
                                        '3':'HCR_BL6-000-R13_ch3_spot_intensity/spots.npy',
                                        '4':'HCR_BL6-000-R13_ch4_spot_intensity/spots.npy'},
                    'multichan_folders': 
                        {'1':{'2':'HCR_BL6-000-R13_ch1_multichannel/ch_channel_4_spots_channel_2.npy',
                              '3':'HCR_BL6-000-R13_ch1_multichannel/ch_channel_4_spots_channel_3.npy',
                              '4':'HCR_BL6-000-R13_ch1_multichannel/ch_channel_4_spots_channel_4.npy'},
                         '2':{'1':'HCR_BL6-000-R13_ch2_multichannel/ch_channel_4_spots_channel_1.npy',
                              '3':'HCR_BL6-000-R13_ch2_multichannel/ch_channel_4_spots_channel_3.npy',
                              '4':'HCR_BL6-000-R13_ch2_multichannel/ch_channel_4_spots_channel_4.npy'},
                         '3':{'1':'HCR_BL6-000-R13_ch3_multichannel/ch_channel_4_spots_channel_1.npy',
                              '2':'HCR_BL6-000-R13_ch3_multichannel/ch_channel_4_spots_channel_2.npy',
                              '4':'HCR_BL6-000-R13_ch3_multichannel/ch_channel_4_spots_channel_4.npy'},
                         '4':{'1':'HCR_BL6-000-R13_ch4_multichannel/ch_channel_4_spots_channel_1.npy',
                              '2':'HCR_BL6-000-R13_ch4_multichannel/ch_channel_4_spots_channel_2.npy',
                              '3':'HCR_BL6-000-R13_ch4_multichannel/ch_channel_4_spots_channel_3.npy',}}
            }
    

    # def get_folder_paths_pipeline(cls) -> dict[str, Dict[str, str]]: 
    #     """Returns folder paths from what is attached in /data/"""

    #     #this happens a little differently in this function. 
    #     #first we find the channel names from the names 

    #     multichannel_regex = r"*_ch(\d{1,3})_multichannel/.*_channel_(\d{1,3})_spots_channel_(\d{1,3}).npy"
    #     spot_regex = r"*_ch(\d{1,3})_spot_intensity/spots.npy"

    #     data_paths = os.listdir(cls.DATA_FOLDER)
        

    #     for datapath in data_paths:
    #         groups = re.search(datapath, multichannel_regex)
            
    #         #find multichannel numbers

    #         #find spot channels

    #         #Make dictionary for spots_folders, multichan_folders with appropriate channels 

    # @classmethod
    # def get_folder_paths_pipeline(cls) -> Dict[str, Dict[str, str]]: #get_folder_paths_pipeline
    #     """Returns folder paths from what is attached in /data/"""
    #     multichannel_regex = r".*_ch(\d{1,3})_multichannel/.*_channel_(\d{1,3})_spots_channel_(\d{1,3}).npy"
    #     spot_regex = r"(\d{1,3}).spots.*\/spots.npy"
    #     # spot_regex = r"(\d{1,3})\/image_data_.*_(\d{1,3})_versus_spots_(\d{1,3})\.csv"

    #     exclude = set(['*.zarr'])
    #     spots_folders = {}
    #     multichan_folders = {}

    #     for root, dirs, files in os.walk(cls.DATA_FOLDER):
    #         # Exclude .zarr directories
    #         dirs[:] = [d for d in dirs if not d.endswith('.zarr')]

    #         for file in files:
    #             # Skip files within .zarr directories
    #             if '.zarr' in root:
    #                 continue
    #             full_path = os.path.join(root, file)
    #             relative_path = os.path.relpath(full_path, cls.DATA_FOLDER)

    #             # Check for spot intensity files
    #             spot_match = re.match(spot_regex, relative_path)
    #             if spot_match:
    #                 channel = spot_match.group(1)
    #                 spots_folders[channel] = relative_path

    #             # Check for multichannel files
    #             multichannel_match = re.match(multichannel_regex, relative_path)
    #             if multichannel_match:
    #                 source_channel = multichannel_match.group(1)
    #                 target_channel = multichannel_match.group(3)
                    
    #                 if source_channel not in multichan_folders:
    #                     multichan_folders[source_channel] = {}
                    
    #                 multichan_folders[source_channel][target_channel] = relative_path

    #     return {
    #         'spots_folders': spots_folders,
    #         'multichan_folders': multichan_folders
    #     }
    @classmethod
    def get_folder_paths_pipeline(cls) -> Dict[str, Dict[str, str]]: #get_folder_paths_pipeline
        """Returns folder paths from what is attached in /data/"""
        spot_regex = r".*(\d{1,3})\/image_data_.*_(\d{1,3})_versus_spots_(\d{1,3})\.csv"
        exclude = set(['*.zarr'])
        spots_folders = {}
        multichan_folders = {}

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
                    source_channel = spot_match.group(2)
                    target_channel = spot_match.group(3)

                    if source_channel == target_channel: 
                        spots_folders[source_channel] = relative_path
                    else:
                        if multichan_folders == {} or source_channel not in multichan_folders.keys():
                            multichan_folders[source_channel]= {target_channel: relative_path}
                        else: 
                            multichan_folders[source_channel][target_channel] = relative_path
        return {
            'spots_folders': spots_folders,
            'multichan_folders': multichan_folders
        }

    @classmethod
    def validate_folder_paths(cls, folder_paths: Dict[str, Dict[str, str]]) -> None:
        """Validates the generated folder paths"""
        expected_channels = set(cls.get_round_channels().keys())
        
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
        if cls.folder_paths == None: 
            folder_paths = cls.get_folder_paths_pipeline()
            cls.validate_folder_paths(folder_paths)
            cls.folder_paths = folder_paths
        else: 
            return cls.folder_paths
        
