import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from .config import Config
from utils.transform_utils import apply_stitching_to_points


class SpotDataLoader:
    def __init__(self):
        self.config = Config
        self.spot_col_order = [
            'spot_id', 'chan', 'chan_spot_id', 'round',
            'z', 'y', 'x', 'z_center', 'y_center', 'x_center', 'dist', 'r'
        ]
        self.tile_col = 'tile_name'
    
    def apply_stitching_transform_to_spots(self, spots_df: pd.DataFrame) -> pd.DataFrame:
        """
        Apply stitching transforms to spot coordinates.
        Preserves original tile-space coordinates and adds stitched global coordinates.
        
        Parameters:
        -----------
        spots_df: pd.DataFrame
            DataFrame with spot data including x, y, z, x_center, y_center, z_center columns
            
        Returns:
        --------
        pd.DataFrame
            DataFrame with additional columns for stitched coordinates and original tile coordinates
        """
        if not self.config.APPLY_STITCHING_TRANSFORM:
            return spots_df
            
        if self.config.STITCHING_XML_PATH is None:
            print("Warning: APPLY_STITCHING_TRANSFORM is True but no STITCHING_XML_PATH provided. Skipping transform.")
            return spots_df
            
        if not Path(self.config.STITCHING_XML_PATH).exists():
            print(f"Warning: XML path does not exist: {self.config.STITCHING_XML_PATH}. Skipping transform.")
            return spots_df
        
        if spots_df.empty:
            return spots_df
            
        # Import transform utilities
        import sys
        utils_path = Path(__file__).parent.parent / 'utils'
        if str(utils_path) not in sys.path:
            sys.path.insert(0, str(utils_path))
        
        tile_name = self.config.CURRENT_TILE
        if tile_name is None:
            print("Warning: No current tile set. Skipping transform.")
            return spots_df
        
        # Store original tile-space coordinates
        spots_df['x_tile'] = spots_df['x'].copy()
        spots_df['y_tile'] = spots_df['y'].copy()
        spots_df['z_tile'] = spots_df['z'].copy()
        spots_df['x_center_tile'] = spots_df['x_center'].copy()
        spots_df['y_center_tile'] = spots_df['y_center'].copy()
        spots_df['z_center_tile'] = spots_df['z_center'].copy()
        
        # Apply transform to spot locations (x, y, z)
        points_xyz = spots_df[['x', 'y', 'z']].values  # N x 3
        transformed_xyz = apply_stitching_to_points(
            points_xyz, 
            tile_name, 
            self.config.STITCHING_XML_PATH
        )
        spots_df['x'] = transformed_xyz[:, 0]
        spots_df['y'] = transformed_xyz[:, 1]
        spots_df['z'] = transformed_xyz[:, 2]
        
        # Apply transform to spot centers (x_center, y_center, z_center)
        points_center = spots_df[['x_center', 'y_center', 'z_center']].values  # N x 3
        transformed_center = apply_stitching_to_points(
            points_center, 
            tile_name, 
            self.config.STITCHING_XML_PATH
        )
        spots_df['x_center'] = transformed_center[:, 0]
        spots_df['y_center'] = transformed_center[:, 1]
        spots_df['z_center'] = transformed_center[:, 2]
        
        print(f"Applied stitching transform to {len(spots_df)} spots from tile {tile_name}")
        
        return spots_df

    def load_multichannel_data(self, ch: str, m_ch: str) -> pd.DataFrame:
        
        
        round_n = self.config.ROUND_N
        multichan_folders = self.config.get_folder_paths()['multichan_folders']
        spot_col_order = ['spot_id','chan','chan_spot_id','round','z','y','x','z_center','y_center','x_center','dist','r'] 
        # spot_cols = ['Z','Y','X','Z_center','Y_center','X_center','dist','r','SEG_ID','FG','BG'] #what comes from CSV 
        #spot_cols = ['z','y','x','z_center','y_center','x_center','dist','r','chan_'+str(ch)+'_fg','chan_'+str(ch)+'_bg','cell_id']
        # spots = pd.DataFrame(np.load(self.config.SPOTS_FOLDER.joinpath(spots_folders[str(ch)])),columns=spot_cols)
        if self.config.CURRENT_TILE is None:
            raise ValueError("Config.CURRENT_TILE is not set. Call Config.set_current_tile before loading data.")

        # add support for case where one channel has NO spots. 
        if multichan_folders[str(ch)].get(str(m_ch)) is None: 
            spots = pd.DataFrame(columns = ['z',
                                 'y',
                                 'x', 
                                 'z_center', 
                                'y_center', 
                                'x_center', 
                                f'chan_{ch}_fg', 
                                f'chan_{ch}_bg',
                                'dist', 
                                'r'
                                ])
        else: 
            spots = pd.read_csv(self.config.SPOTS_FOLDER.joinpath(multichan_folders[str(ch)][str(m_ch)]))
            spots = spots.rename(columns = {'Z':'z',
                                    'Y': 'y',
                                    'X': 'x', 
                                    'Z_center': 'z_center', 
                                    'Y_center': 'y_center', 
                                    'X_center': 'x_center', 
                                    #  'SEG_ID': 'cell_id', 
                                    'FG': f'chan_{ch}_fg', 
                                    'BG': f'chan_{ch}_bg'})
        spots['round']=str(round_n)
        spots['chan']=str(ch)
        spots['spot_id']=range(1, len(spots)+1)
        spots['chan_spot_id']=range(1, len(spots)+1)
        spots[self.tile_col] = self.config.CURRENT_TILE
        base_cols = spot_col_order + [self.tile_col]
        remaining_cols = [c for c in spots.columns if c not in base_cols]
        spots = spots[base_cols + remaining_cols]
        return spots

    def load_detected_spots_for_channel(self, ch):
        round_n = self.config.ROUND_N
        spots_folders = self.config.get_folder_paths()['spots_folders']
        spot_col_order = ['spot_id','chan','chan_spot_id','round','z','y','x','z_center','y_center','x_center','dist','r'] 
        # spot_cols = ['Z','Y','X','Z_center','Y_center','X_center','dist','r','SEG_ID','FG','BG'] #what comes from CSV 
        # add support for case if no spots in channel are detected; this means that 
        # there will be no stats file for itself.... eg. 488 vs 488.csv 
        # in this case should we make an empty table?  
        if spots_folders.get(str(ch)) is None: 
            spots = pd.DataFrame(columns = ['z',
                                 'y',
                                 'x', 
                                 'z_center', 
                                'y_center', 
                                'x_center', 
                                f'chan_{ch}_fg', 
                                f'chan_{ch}_bg',
                                'dist', 
                                'r'
                                ])
        else: 
            spots = pd.read_csv(self.config.SPOTS_FOLDER.joinpath(spots_folders[str(ch)]))
            
            spots = spots.rename(columns = {'Z':'z',
                                    'Y': 'y',
                                    'X': 'x', 
                                    'Z_center': 'z_center', 
                                    'Y_center': 'y_center', 
                                    'X_center': 'x_center', 
                                    #  'SEG_ID': 'cell_id', 
                                    'FG': f'chan_{ch}_fg', 
                                    'BG': f'chan_{ch}_bg'})
        spots['round']=str(round_n)
        spots['chan']=str(ch)
        spots['spot_id']=range(1, len(spots)+1)
        spots['chan_spot_id']=range(1, len(spots)+1)
        if self.config.CURRENT_TILE is None:
            raise ValueError("Config.CURRENT_TILE is not set. Call Config.set_current_tile before loading data.")
        spots[self.tile_col] = self.config.CURRENT_TILE
        base_cols = spot_col_order + [self.tile_col]
        remaining_cols = [c for c in spots.columns if c not in base_cols]
        spots = spots[base_cols + remaining_cols]
        return spots


    def load_all_spots(self) -> pd.DataFrame:
        """Load and combine all spot data"""


        channels = self.config.get_round_spot_channels()
        channels = [str(i) for i in list(channels) if i!= '405']

        round_n = self.config.ROUND_N
        multichan_folders = self.config.get_folder_paths()['multichan_folders']
        if self.config.CURRENT_TILE is None:
            raise ValueError("Config.CURRENT_TILE is not set. Call Config.set_current_tile before loading data.")

        channel_spots = {}
        #get detected spots
        for ch in channels:
            channel_spots[ch] = self.load_detected_spots_for_channel(ch)


        for ch in channels:
            # Luminance from loop channel at channel 3's spot locations
            # if len(ch) ==1:
            #     # check if ch is single digit or 
            #     if str(ch) != '1':
            #         m_ch = 1
            #         multi_spot_cols = ['z','y','x','chan_'+str(ch)+'_fg','chan_'+str(ch)+'_bg']
            #         chan_multichan_df = pd.DataFrame(np.load(self.config.SPOTS_FOLDER.joinpath(multichan_folders[str(ch)][str(m_ch)])),columns=multi_spot_cols)
            #         channel_spots[str(m_ch)] = channel_spots[str(m_ch)].merge(chan_multichan_df, on=['z', 'y', 'x'], how='inner')
            #     # Luminance from loop channel at channel 2's spot locations
            #     if str(ch) != '2':
            #         m_ch = 2
            #         multi_spot_cols = ['z','y','x','chan_'+str(ch)+'_fg','chan_'+str(ch)+'_bg']
            #         chan_multichan_df = pd.DataFrame(np.load(self.config.SPOTS_FOLDER.joinpath(multichan_folders[str(ch)][str(m_ch)])),columns=multi_spot_cols)
            #         channel_spots[str(m_ch)] = channel_spots[str(m_ch)].merge(chan_multichan_df, on=['z', 'y', 'x'], how='inner')
            #     # Luminance from loop channel at channel 3's spot locations
            #     if (str(ch) != '3') & (round_n != 0):
            #         m_ch = 3
            #         multi_spot_cols = ['z','y','x','chan_'+str(ch)+'_fg','chan_'+str(ch)+'_bg']
            #         chan_multichan_df = pd.DataFrame(np.load(self.config.SPOTS_FOLDER.joinpath(multichan_folders[str(ch)][str(m_ch)])),columns=multi_spot_cols)
            #         channel_spots[str(m_ch)] = channel_spots[str(m_ch)].merge(chan_multichan_df, on=['z', 'y', 'x'], how='inner')
            #     # Luminance from loop channel at channel 4's spot locations
            #     if str(ch) != '4':
            #         m_ch = 4
            #         multi_spot_cols = ['z','y','x','chan_'+str(ch)+'_fg','chan_'+str(ch)+'_bg']
            #         chan_multichan_df = pd.DataFrame(np.load(self.config.SPOTS_FOLDER.joinpath(multichan_folders[str(ch)][str(m_ch)])),columns=multi_spot_cols)
            #         channel_spots[str(m_ch)] = channel_spots[str(m_ch)].merge(chan_multichan_df, on=['z', 'y', 'x'], how='inner')
            # else: 

            active_channel_list = channels
            active_channel_list_without_current_ch = [chan for chan in channels if chan!=ch]
            for m_ch in active_channel_list_without_current_ch: 
                
                try: 
                    multi_spot_cols = ['z','y','x','chan_'+str(ch)+'_fg','chan_'+str(ch)+'_bg']
                    chan_multichan_df = self.load_multichannel_data(ch, m_ch)
                    drop_cols = ['spot_id', 'chan', 'chan_spot_id', self.tile_col]
                    existing_drop_cols = [col for col in drop_cols if col in chan_multichan_df.columns]
                    test_multichan = chan_multichan_df.drop(existing_drop_cols, axis = 1)

                    channel_spots[str(m_ch)] = channel_spots[m_ch].merge(test_multichan, on = ['z', 'y', 'x', 'round', 'z_center', 'y_center', 'x_center', 'dist', 'r'], how = 'inner')
                except Exception as e: 
                    print(f"Error processing channel {ch} at measurement channel {m_ch}: {str(e)}")
                    continue

        mixed_spots_df = pd.DataFrame()
        #convert dict of channel_spots to list
        spots_list = list(channel_spots.values())
        for i, ch in enumerate(channels):
            mixed_spots_df = pd.concat([mixed_spots_df,spots_list[i] ])
        for ch in channels:
            mixed_spots_df['chan_'+str(ch)+'_intensity'] = mixed_spots_df['chan_'+str(ch)+'_fg']-mixed_spots_df['chan_'+str(ch)+'_bg']

        if not mixed_spots_df.empty:
            mixed_spots_df[self.tile_col] = self.config.CURRENT_TILE
            
        # Apply stitching transforms if enabled
        # mixed_spots_df = self.apply_stitching_transform_to_spots(mixed_spots_df)
            
        return mixed_spots_df    