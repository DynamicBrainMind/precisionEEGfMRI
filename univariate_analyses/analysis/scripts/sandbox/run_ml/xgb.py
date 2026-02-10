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



def rmse(y_true, y_pred):
    return np.sqrt(np.mean((y_true - y_pred)**2))

def get_p(surrogates, true_loss):
    return (np.sum(surrogates >= true_loss) + 1) / (len(surrogates) + 1)

def get_loss(subjects, network, surrogate):
	
    losses = []

    def prep_subject(subject, network, surrogate):
        X_train, y_train, X_test, y_test = train_test_across(network=network, held_out_subject=subject)
        
        if surrogate:
            np.random.shuffle(y_train)
        return X_train, y_train, X_test, y_test


    # CPU parallelize the data prep
    with ThreadPoolExecutor(max_workers=16) as pool:
        futures = [pool.submit(prep_subject, subject, network, surrogate) for subject in subjects]
        for fut in futures:
            X_train, y_train, X_test, y_test = fut.result()
            # XGBoost
            xgb = XGBRegressor(n_estimators=300,
                               max_depth=3,
                               learning_rate=.05,
                               subsample=.8,
                               colsample_bytree=.8,
                               tree_method='gpu_hist',
                               predictor='gpu_predictor',
                               verbosity=0)
            xgb.fit(X_train, y_train)
            y_pred = xgb.predict(X_test)
            losses.append(r2_score(y_test, y_pred))

    return np.array(losses), xgb


def get_shap(model, network):
    X, y = train_test_across(network=network)
    model.fit(X, y)

    explainer = shap.TreeExplainer(model)

    shap_values = explainer.shap_values(X)

    return shap_values
    

def get_surrogates(nsims, networks, subjects, 
            out_path=Path('scripts/sandbox/run_ml/surrogates.pkl')):

    if not os.path.exists(out_path):
        out = {x: [] for x in networks}
    else:
        with open(out_path, 'rb') as file:
            out = pickle.load(file)

    start = len(out[networks[0]]) + 1

    for sim in range(start, nsims+1):
        print(f'\n---SIMULATION {sim} of {nsims}---\n')

        for network in networks:
            print(f'Computing network {network}')

            loss, _ = get_loss(subjects, network, surrogate=True)

            out[network].append(loss.mean())

        # Write after each sim
        with open(out_path, 'wb') as file:
            pickle.dump(out, file)




# --- GET SUBJECTS AND NETWORKS --- #
dpath = Path('data/formatted/reduced')
subjects = get_subjects(dpath)

dpath = Path('data/formatted/reduced/sub-001.pkl')

with open(dpath, 'rb') as file:
    d = pickle.load(file)

networks = d['ses-001']['run-001']['y'].keys()
networks = [x for x in networks if x not in ['DAN', 'DMN']]



# --- GET AND WRITE SURROGATES --- #


get_surrogates(500, networks, subjects)




# --- GET AND WRITE TRUTH --- #
out = {}
for network in networks:
    # XGBoost
    xgb = XGBRegressor(n_estimators=300,
                       max_depth=3,
                       learning_rate=.05,
                       subsample=.8,
                       colsample_bytree=.8,
                       tree_method='gpu_hist',
                       predictor='gpu_predictor',
                       verbosity=0)

    out[network] = get_shap(xgb, network)
    

with open(root / Path('shap_networks.pkl'), 'wb') as file:
    pickle.dump(out, file)


# Get truth
truth = {}
for network in networks:
    print(f'Processing network {network}')
    truth[network], _ = get_loss(subjects, network, surrogate=False)

with open(root / Path('truth.pkl'), 'wb') as file:
    pickle.dump(truth, file)


out = {'truth': truth, 'shap': get_shap(xgb)}
with open(root / Path('truth.pkl'), 'wb') as file:
    pickle.dump(out, file)






