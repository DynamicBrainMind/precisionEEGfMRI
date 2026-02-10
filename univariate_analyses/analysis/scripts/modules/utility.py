from pathlib import Path
from glob import glob

def get_subjects(dpath):

    subjects = [Path(x).stem for x in glob(str(dpath / Path('*')))]
    subjects = sorted(subjects, key=lambda x: int(x.split('-')[1]))

    return subjects
