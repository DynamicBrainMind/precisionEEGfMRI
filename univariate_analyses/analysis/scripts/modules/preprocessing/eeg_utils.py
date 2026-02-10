import mne
from pathlib import Path
from pyprojroot import here
from collections import OrderedDict
import pandas as pd

def get_channel_coordinates(channels,
                            locations='analysis/scripts/Standard-10-20-Cap81.locs'):    
    '''
    Gives x, y, z coordinates for EEG channels for plotting
    
    --- PARAMETERS ---
    ------------------
    channels (list of str): 10-20 names of channels to get coordinates for
    locations (str): Path to Standard-10-20-Cap81.locs
    '''
    
    # Load montage at assumed location
    file = Path(locations)
    montage = mne.channels.read_custom_montage(here() / file)
    ch_pos = montage.get_positions()['ch_pos']
    ch_pos = OrderedDict((key, value*1000) for key, value in ch_pos.items())
    ch_pos = pd.DataFrame(ch_pos).transpose()
    ch_pos.columns = ['x', 'y', 'z']
    ch_pos.insert(0, 'channel', ch_pos.index)
    ch_pos = ch_pos[ch_pos['channel'].isin(channels)]
    
    return ch_pos

