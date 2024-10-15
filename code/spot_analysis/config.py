import pathlib
from typing import Dict, Any
import os
import re


class Config:
    # Path configurations
    SPOTS_FOLDER = pathlib.Path('/data/')
    DATA_FOLDER = pathlib.Path('/data/')
    OUTPUT_FOLDER = pathlib.Path('/results/')
    OUTPUT_DATA_TYPE = 'zarr'
    SCRATCH_FOLDER = pathlib.Path('/scratch/')
    
    # Processing parameters
    ROUND_N = 13
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

    # Gene dictionary
    GENE_DICT: Dict[str, Dict[str, str]] = {'0':{'1': 'Vip', '2': 'Sst', '4': 'Slc17a7'},
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
    
    @classmethod
    def get_round_channels(cls) -> Dict[str, str]:
        return cls.GENE_DICT[str(cls.ROUND_N)]
    
    
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

    @classmethod
    def get_folder_paths(cls) -> Dict[str, Dict[str, str]]: #get_folder_paths_pipeline
        """Returns folder paths from what is attached in /data/"""
        multichannel_regex = r".*_ch(\d{1,3})_multichannel/.*_channel_(\d{1,3})_spots_channel_(\d{1,3}).npy"
        spot_regex = r".*_ch(\d{1,3})_spot_intensity/spots.npy"

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
                    channel = spot_match.group(1)
                    spots_folders[channel] = relative_path

                # Check for multichannel files
                multichannel_match = re.match(multichannel_regex, relative_path)
                if multichannel_match:
                    source_channel = multichannel_match.group(1)
                    target_channel = multichannel_match.group(3)
                    
                    if source_channel not in multichan_folders:
                        multichan_folders[source_channel] = {}
                    
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
        folder_paths = cls.get_folder_paths_pipeline()
        cls.validate_folder_paths(folder_paths)
        return folder_paths
