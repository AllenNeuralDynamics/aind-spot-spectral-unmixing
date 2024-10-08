import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple
from .config import Config

class SpotDataLoader:
    def __init__(self):
        self.config = Config
        self.spot_col_order = [
            'spot_id', 'chan', 'chan_spot_id', 'cell_id', 'round',
            'z', 'y', 'x', 'z_center', 'y_center', 'x_center', 'dist', 'r'
        ]
    
    def load_channel_spots(self, channel: str) -> pd.DataFrame:
        """Load spots data for a specific channel"""
        spot_cols = [
            'z', 'y', 'x', 'z_center', 'y_center', 'x_center',
            'dist', 'r', f'chan_{channel}_fg', f'chan_{channel}_bg', 'cell_id'
        ]
        
        spots_path = self.config.SPOTS_FOLDER / self.config.get_folder_paths()['spots_folders'][channel]
        spots_data = pd.DataFrame(
            np.load(spots_path),
            columns=spot_cols
        )
        
        spots_data['round'] = str(self.config.ROUND_N)
        spots_data['chan'] = channel
        spots_data['spot_id'] = range(1, len(spots_data) + 1)
        spots_data['chan_spot_id'] = range(1, len(spots_data) + 1)
        
        return spots_data[list(self.spot_col_order) + list(np.setdiff1d(spots_data.columns, self.spot_col_order))]
    

    #currently broken... only returning a single item from each list? 
    def load_multichannel_data(self, spots_df: pd.DataFrame, channel: str) -> pd.DataFrame:
        """Load and merge multichannel data for a specific channel"""
        folder_paths = self.config.get_folder_paths()['multichan_folders']
        
        for other_chan in self.config.get_round_channels().keys():
            if other_chan != channel:
                multi_spot_cols = ['z', 'y', 'x', f'chan_{other_chan}_fg', f'chan_{other_chan}_bg']
                multichan_path = self.config.SPOTS_FOLDER / folder_paths[channel][other_chan]
                
                multichan_df = pd.DataFrame(
                    np.load(multichan_path),
                    columns=multi_spot_cols
                )
                # multichan_df['chan'] = other_chan
                spots_df = spots_df.merge(multichan_df, how='inner',  on=['z', 'y', 'x'] )
        
        return spots_df
    
    def load_detected_spots_for_channel(self, ch):
        round_n = self.config.ROUND_N
        spots_folders = self.config.get_folder_paths()['spots_folders']
        spot_col_order = ['spot_id','chan','chan_spot_id','cell_id','round','z','y','x','z_center','y_center','x_center','dist','r']
        spot_cols = ['z','y','x','z_center','y_center','x_center','dist','r','chan_'+str(ch)+'_fg','chan_'+str(ch)+'_bg','cell_id']
        spots = pd.DataFrame(np.load(self.config.SPOTS_FOLDER.joinpath(spots_folders[str(ch)])),columns=spot_cols)
        spots['round']=str(round_n)
        spots['chan']=str(ch)
        spots['spot_id']=range(1, len(spots)+1)
        spots['chan_spot_id']=range(1, len(spots)+1)
        spots = spots[list(spot_col_order)+list(np.setdiff1d(spots.columns,spot_col_order))]
        return spots


    def load_all_spots(self) -> pd.DataFrame:
        """Load and combine all spot data"""
        # all_spots = []
        # channels = self.config.get_round_channels().keys()
        
        # for channel in channels:
        #     channel_spots = self.load_channel_spots(channel)
        #     channel_spots = self.load_multichannel_data(channel_spots, channel)
        #     all_spots.append(channel_spots)
        
        # return pd.concat(all_spots, ignore_index=True)

        channels = self.config.get_round_channels().keys()
        round_n = self.config.ROUND_N
        multichan_folders = self.config.get_folder_paths()['multichan_folders']

        channel_spots = {}
        #get detected spots
        for ch in channels:
            channel_spots[ch] = self.load_detected_spots_for_channel(ch)

        for ch in channels:
            # Luminance from loop channel at channel 3's spot locations
            if str(ch) != '1':
                m_ch = 1
                multi_spot_cols = ['z','y','x','chan_'+str(ch)+'_fg','chan_'+str(ch)+'_bg']
                chan_multichan_df = pd.DataFrame(np.load(self.config.SPOTS_FOLDER.joinpath(multichan_folders[str(ch)][str(m_ch)])),columns=multi_spot_cols)
                channel_spots[str(m_ch)] = channel_spots[str(m_ch)].merge(chan_multichan_df, on=['z', 'y', 'x'], how='inner')
            # Luminance from loop channel at channel 2's spot locations
            if str(ch) != '2':
                m_ch = 2
                multi_spot_cols = ['z','y','x','chan_'+str(ch)+'_fg','chan_'+str(ch)+'_bg']
                chan_multichan_df = pd.DataFrame(np.load(self.config.SPOTS_FOLDER.joinpath(multichan_folders[str(ch)][str(m_ch)])),columns=multi_spot_cols)
                channel_spots[str(m_ch)] = channel_spots[str(m_ch)].merge(chan_multichan_df, on=['z', 'y', 'x'], how='inner')
            # Luminance from loop channel at channel 3's spot locations
            if (str(ch) != '3') & (round_n != 0):
                m_ch = 3
                multi_spot_cols = ['z','y','x','chan_'+str(ch)+'_fg','chan_'+str(ch)+'_bg']
                chan_multichan_df = pd.DataFrame(np.load(self.config.SPOTS_FOLDER.joinpath(multichan_folders[str(ch)][str(m_ch)])),columns=multi_spot_cols)
                channel_spots[str(m_ch)] = channel_spots[str(m_ch)].merge(chan_multichan_df, on=['z', 'y', 'x'], how='inner')
            # Luminance from loop channel at channel 4's spot locations
            if str(ch) != '4':
                m_ch = 4
                multi_spot_cols = ['z','y','x','chan_'+str(ch)+'_fg','chan_'+str(ch)+'_bg']
                chan_multichan_df = pd.DataFrame(np.load(self.config.SPOTS_FOLDER.joinpath(multichan_folders[str(ch)][str(m_ch)])),columns=multi_spot_cols)
                channel_spots[str(m_ch)] = channel_spots[str(m_ch)].merge(chan_multichan_df, on=['z', 'y', 'x'], how='inner')


        mixed_spots_df = pd.DataFrame()
        #convert dict of channel_spots to list
        spots_list = list(channel_spots.values())
        for i, ch in enumerate(channels):
            mixed_spots_df = pd.concat([mixed_spots_df,spots_list[i] ])
        for ch in channels:
            mixed_spots_df['chan_'+str(ch)+'_intensity'] = mixed_spots_df['chan_'+str(ch)+'_fg']-mixed_spots_df['chan_'+str(ch)+'_bg']
        return mixed_spots_df    