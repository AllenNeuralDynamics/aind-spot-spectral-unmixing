import pathlib
from typing import Dict, Any

class Config:
    # Path configurations
    SPOTS_FOLDER = pathlib.Path('/data/')
    OUTPUT_FOLDER = pathlib.Path('/results/')
    OUTPUT_DATA_TYPE = 'zarr'
    
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
    def get_folder_paths(cls) -> Dict[str, Dict[str, str]]:
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
                