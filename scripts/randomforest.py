"""
This script is a Real Bogus classifier.
"""
from collections import Counter
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.metrics import (
        accuracy_score,
        classification_report,
        ConfusionMatrixDisplay,
        confusion_matrix,
        RocCurveDisplay,
)

from imblearn.under_sampling import RandomUnderSampler

from config import Configuration
from libraries.utils import Utils

from pathlib import Path
from datetime import datetime

import os
import re
import joblib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# Random seed
SEED = 43

# Load the dataset
def load_rb_catalog(filepath):
    return pd.read_csv(filepath)


def make_big_dataset(rb_path):
    rb_catalogs = os.listdir(rb_path)
    base = Path("/Users/mcast/Research/data_directory/training/injected/rb_catalogs/2024-10-11/FIELD_36.007")

    for catalog_name in rb_catalogs:
        filepath = base / catalog_name
        rb_catalog  = load_rb_catalog(catalog)
        X, y = format_dataset(rb_catalog)


# Format dataset
def format_dataset(rb_catalog):
    X = rb_catalog.drop(columns=['REAL'])
    y = rb_catalog['REAL']
    return X, y

def conform(rb_catalog, train=False):
    # Using this column order to conform to O2 paper
    columns = ['FLUX_APER', 'FLUXERR_APER',
                'MAG_APER', 'MAGERR_APER',
                'FLUX_MAX',
                'ISOAREA_IMAGE', 'ISOAREAF_IMAGE',
                'X2_IMAGE', 'Y2_IMAGE', 'XY_IMAGE', 'ERRX2_IMAGE', 'ERRY2_IMAGE', 'ERRXY_IMAGE',
                'CXX_IMAGE', 'CYY_IMAGE', 'CXY_IMAGE', 'ERRCXX_IMAGE', 'ERRCYY_IMAGE', 'ERRCXY_IMAGE',
                'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE', 'ERRA_IMAGE', 'ERRB_IMAGE', 'ERRTHETA_IMAGE',
                'ISO0', 'ISO1', 'ISO2', 'ISO3', 'ISO4', 'ISO5', 'ISO6', 'ISO7',
                'FLAGS',
                'FWHM_IMAGE',
                'ELONGATION',
                'ELLIPTICITY',
                'POLAR_IMAGE',
                'VIGNET', 'VIGNET_SHIFT',
                'FLUX_GROWTHSTEP',
                'MAG_GROWTH', 'MAG_GROWTHSTEP',
                'FLUX_RADIUS',
               ]
    # number of columns used in O2 paper plus labels
    total_columns = 44

    if train:
        columns.append('REAL')
        total_columns += 1

    if len(columns) != total_columns:
        print("Columns did not conform")

    return rb_catalog[columns]

def train() -> None:
    base = Path("/Users/mcast/Research/data_directory/training/injected/rb_catalogs/2024-10-11/FIELD_36.007")
    filepath = base / "FIELD_36.007_300s_12000x10600_331_bkfcspiad_realbogus.csv"

    # Load dataset
    print("LOADING DATASET...")
    rb_catalog = load_rb_catalog(filepath)

    # Conform to O2 paper
    print("CONFORM...")
    rb_catalog = conform(rb_catalog, train=True)

    # Format dataset
    print("FORMATTING DATASET...")
    X, y = format_dataset(rb_catalog)

    # Make dataset test/train split
    print("MAKING TRAIN TEST SETS...")
    X_train, X_test, y_train, y_test = train_test_split(X,
                                                        y,
                                                        test_size=0.2,
                                                        stratify=y,
                                                        random_state=SEED)

    print('Original dataset shape %s' % Counter(y_train))

    # Balance the data set by undersampling
    # this matches the majority class to the minority class
    rus = RandomUnderSampler(random_state=SEED)
    print("Balance the dataset")

    X_train, y_train = rus.fit_resample(X_train, y_train)
    print('Resampled dataset shape %s' % Counter(y_train))

    # initialize RandomForestClassifier
    print("Initializing Random Forest Classifier.")
    TREES = 10
    MAX_FEATURES = 7
    MIN_SAMPLES = 20
    CRITERION = "gini"
    CPU_LIMIT = 1

    rf = RandomForestClassifier(
            n_estimators=TREES,
            max_features=MAX_FEATURES,
            min_samples_leaf=MIN_SAMPLES,
            criterion=CRITERION,
            random_state=SEED,
            n_jobs=CPU_LIMIT)

    # Train the model
    print("Training the model!")
    rf.fit(X_train, y_train)

    # Evaluate
    print("EVALUATE")
    y_pred = rf.predict(X_test)
    print(f"Accuracy: {accuracy_score(y_test, y_pred):.4f}\n")
    report = classification_report(y_test, y_pred, target_names=['bogus', 'real'],
                                   output_dict=True)
    print(classification_report(y_test, y_pred, target_names=['bogus', 'real']))

    # Confusion Matrix
    cm = confusion_matrix(y_test, y_pred)
    cm_display = ConfusionMatrixDisplay(cm).plot()
    plt.title(f"Confusion Matrix of Test Set\n bogus: {report['bogus']['support']} real: {report['real']['support']}\n Accuracy Score: {accuracy_score(y_test, y_pred):.4f} ")
    plt.savefig("ConfusionMatrix.png", dpi=300, bbox_inches="tight")
    plt.show()

    # ROC - Receiver Operating Curve
    rf.fit(X_train, y_train)
    ax = plt.gca()
    plt.grid(True)
    rf_disp = RocCurveDisplay.from_estimator(
        rf, X_test, y_test, ax=ax, curve_kwargs=dict(alpha=0.8)
    )
    ax.set_xlim(0.0, 0.5)
    ax.set_ylim(0.5, 1.0)
    ax.set_box_aspect(1)
    plt.title("Receiver Operating Curve")
    plt.savefig("ROC.png", dpi=300, bbox_inches="tight")
    plt.show()

    # Feature importances radial plot
    importances = rf.feature_importances_
    feature_names = X_train.columns
    categories = list(feature_names)
    N = len(categories)

    # Compute angle for each axis
    angles = [n / float(N) * 2 * np.pi for n in range(N)]
    angles += angles[:1]

    # Values
    values = list(importances)
    values += values[:1]

    # Radial plot
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(polar=True))
    ax.fill(angles, values, 'b', alpha=0.1)
    ax.set_rorigin(-0.02)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, fontsize=8)
    ax.spines['polar'].set_visible(False)
    ax.plot(angles, values, linewidth=1, linestyle='solid', label='Importance')
    plt.tick_params(axis='x', which='major', pad=10)
    plt.yticks([0.02, 0.04, 0.06, 0.08, 0.10, 0.12],
               ha='center', va='center')
    plt.title("Feature Importance")
    plt.savefig("feature_importance.png", dpi=300, bbox_inches="tight")
    plt.show()


    # Save trained model
    basename = 'toros_rf_classifier'
    timestamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    extension = 'model'
    filename = f'{basename}_{timestamp}.{extension}'
    joblib.dump(rf, filename)
    # run env = randomforest.main()
    # then set globals().update(env) to get output of main()
    return rf

def predict(rf_model, filepath):
    # Load data to process
    data = load_rb_catalog(filepath)

    # conform data to rf model
    X = conform(data)

    # predict
    predictions = rf_model.predict(X)

    # add precition to data
    data['CLASS_REAL'] = predictions

    # update filepath to save classification
    filepath = filepath.replace('/rb_catalogs/', '/class_rf/')
    filepath = filepath.replace('_realbogus.csv','_classification.csv')
    # save classification as new file
    data.to_csv(filepath)
    print(f"{filepath=}")
    return predictions

def get_rb_catalog_list():
    files, date_dirs = Utils.get_all_files_per_field(
            Configuration.REALBOGUS_CATALOG_DIRECTORY,
            Configuration.FIELD,
            'realbogus',
            '.csv')

    return files, date_dirs

def load_model():
    # look to see if model exists
    directory = Path()
    basename = "toros_rf_classifier"
    extension = "model"
    models = list(directory.glob(f"{basename}*.{extension}"))

    if models:
        print("Models Found!")
        newest_model = None
        latest_timestamp = datetime.min
        for model in models:
            timestamp = model.stem[-15:]
            timestamp = datetime.strptime(timestamp, "%Y%m%dT%H%M%S")
            if timestamp > latest_timestamp:
                latest_timestamp = timestamp
                newest_model = model
        print(f"Loading {newest_model}...")
        rf_model = joblib.load(newest_model)
    else:
        print("NO Models Found... training new model!")
        rf_model = train()

    return rf_model

def generate_rf_classification_directories():
    # get the file list for all dates the FIELD was observed
    Utils.log("Generating directories for RANDOM FOREST classifications", "info")
    files, date_dirs = get_rb_catalog_list()

    class_rf_dir = Configuration.CLASSIFICATION_RF_DIRECTORY

    # make the output directories
    output_dirs = []
    for date in date_dirs:
        output_dirs.append(class_rf_dir)
        output_dirs.append(class_rf_dir + date)
        output_dirs.append(class_rf_dir + date + "/" + Configuration.FIELD)

    Utils.create_directories(output_dirs)

    return files, date_dirs


def main():
    # load rf_model
    rf_model = load_model()

    # get rb catalogs
    files, date_dirs = get_rb_catalog_list()

    # generate classification dirs to store output
    generate_rf_classification_directories()

    for file in files:
        predictions = predict(rf_model, file)

    return locals()

if __name__ == "__main__":
    main()
