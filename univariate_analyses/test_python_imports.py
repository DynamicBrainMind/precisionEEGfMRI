try:
    from pyprojroot import here
    import os
    from pathlib import Path
    from tqdm import tqdm
    import sys
    modules_root = here() / Path('analysis/scripts/modules')
    os.chdir(here() / Path('analysis'))
    sys.path.append(str(here()))
    sys.path.append(str(modules_root))
    from modeling.model_selection import train_test_across
    from modeling.dummy import BaselineModel
    from utility import get_subjects
    root = Path(modules_root / Path('../sandbox/run_ml'))
    import numpy as np
    import shap
    from concurrent.futures import ThreadPoolExecutor
    import pandas as pd
    from xgboost import XGBRegressor
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.metrics import r2_score
    import pickle
    from glob import glob

    print('ALL IMPORTS COMPLETED SUCCESSFULLY')

except Exception as e:
    print('Error installing library:')
    print('\n')
    print(e)

