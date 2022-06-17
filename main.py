import json
import os
import shutil
import sys
import time
from argparse import ArgumentParser
import joblib
import pandas as pd
from ItvAnt.main import annotate
from pathlib import Path


def load_data(path) -> pd.DataFrame:
    data = pd.read_csv(path, index_col=0)
    with open(Path('models') / 'features.json') as handle:
        features = json.load(handle)['best_comb']
    columns = data.columns
    for feature in features:
        if feature not in columns:
            data[feature] = 0
    return data.drop(['start', 'end', 'chr'], axis=1)[features]


def predict(data: pd.DataFrame, threshold: float) -> pd.DataFrame:
    model = joblib.load(Path('models') / 'core.model')
    score = model.predict_proba(data)[:, 1]
    label = score > threshold
    return pd.DataFrame({'repeats_id': data.index, 'probability_score': score, 'pathogenic': label})


def get_args():
    ap = ArgumentParser(description='predict the pathogenic extension mutation')
    ap.add_argument('input', help='BED6 file of repeats to be predicted')
    ap.add_argument('--output', '-o', help='the path to save the result', required=True)
    return ap.parse_args()


def split_name(file) -> tuple:
    prefix, full_name = os.path.split(file)
    name, _ = os.path.splitext(full_name)
    return prefix, name


if __name__ == '__main__':
    args = get_args()
    path_name = args.input
    task_dir = str(time.time_ns()) + str(hash(path_name))
    try:
        os.makedirs(task_dir)
    except FileExistsError:
        print('a same task is running!!!')
        sys.exit(1)
    annotate(path_name, task_dir)
    tb = load_data(Path(task_dir) / 'annotated.csv')
    result = predict(tb, 0.04)
    shutil.rmtree(task_dir)
    _, input_name = split_name(path_name)
    os.makedirs(args.output, exist_ok=True)
    result.to_csv(Path(args.output) / '{}.csv'.format(input_name), index=False)
