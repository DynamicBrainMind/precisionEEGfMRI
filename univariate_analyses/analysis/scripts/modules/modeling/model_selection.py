from pathlib import Path
import numpy as np
from glob import glob
import pickle

def train_test_across(held_out_subject):
    # Held out subject in form sub-\d\d\d

    # Get data
    droot = Path('data/formatted/reduced')
    subjects = [Path(x).stem for x in glob(str(droot / Path('*')))]
    subjects = sorted(subjects, key = lambda x: int(x.split('-')[1]))
    subjects_train = [x for x in subjects if x != held_out_subject]

    X_l = []
    y_l = []

    def concat(d):
        X_l = []
        y_l = []
        for session in d:
            runs = list(d[session].keys())
            X_ses = np.concatenate([d[session][x]['X'] for x in runs], axis=0)
            y_ses = np.concatenate([d[session][x]['y']['DAN'] for x in runs], axis=0)

            X_l.append(X_ses)
            y_l.append(y_ses)

        X_sub = np.concatenate(X_l, axis=0)
        y_sub = np.concatenate(y_l, axis=0)

        return X_sub, y_sub

    for subject in subjects_train:

        dpath = droot / Path(f'{subject}.pkl')
        with open(dpath, 'rb') as file:
            d = pickle.load(file)

        X_sub, y_sub = concat(d)

        X_l.append(X_sub)
        y_l.append(y_sub)

    X_train = np.concatenate(X_l, axis=0)
    y_train = np.concatenate(y_l, axis=0)

    dpath = droot / Path(f'{held_out_subject}.pkl')
    with open(dpath, 'rb') as file:
        d_test = pickle.load(file)

    X_test, y_test = concat(d_test)


    return X_train, y_train, X_test, y_test

