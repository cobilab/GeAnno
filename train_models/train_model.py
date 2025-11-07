from concurrent.futures import ThreadPoolExecutor
from sklearn.preprocessing import StandardScaler, LabelEncoder, label_binarize
from sklearn.pipeline import Pipeline
from imblearn.over_sampling import SMOTE
from imblearn.under_sampling import RandomUnderSampler
from sklearn.model_selection import RandomizedSearchCV, train_test_split
from sklearn.metrics import average_precision_score, classification_report, roc_auc_score
from sklearn.ensemble import RandomForestClassifier
from collections import Counter
from matplotlib import pyplot as plt
from xgboost import XGBClassifier

import os, logging, joblib, argparse
import seaborn as sea
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format='[PRED] %(message)s')

MOTIFS = "A,C,T,G,ATG,TAA,TAG,TGA,TATAAA,CAAT,CGCG,GT,AG,AAA,TTT,GGG,CGG,CAG"

plant_dir = os.environ.get("PLANT_DIR")
if not plant_dir:
    raise EnvironmentError("PLANT_DIR environment variable is not set")

def import_data_everything(dataset):
    data = pd.read_csv(dataset)
    logging.info(f"Read CSV: {data.shape}")

    data = data.sample(frac=1, random_state=42).reset_index(drop=True)  # Shuffle
    logging.info(f"After shuffle: {data.shape}")

    data = data.drop(columns=["species"])

    data = data.dropna()
    logging.info(f"After dropna: {data.shape}")

    print("Original size:", len(data))
    logging.info(f"Original size: {len(data)}")

    plt.figure(figsize=(6, 6)) 
    sea.countplot(x="label", hue="label", data=data, palette="Set3", legend=False)

    plt.title("Before filtering", fontsize=12)
    plt.xticks(rotation=45, ha='right') 
    plt.tight_layout()
    plt.show()
    
    return data

    
def import_split_scale_data(dataset):
    le = LabelEncoder()

    dataset = dataset.sample(frac=1, random_state=42).reset_index(drop=True)

    features = dataset.drop(columns=["label"])
    dataset["_features_tuple"] = features.apply(lambda row: tuple(row), axis=1)
    dataset = dataset.drop_duplicates(subset="_features_tuple", keep="first")
    dataset = dataset.drop(columns=["_features_tuple"])

    dataset["label"] = le.fit_transform(dataset["label"])

    X = dataset.drop(['label'], axis=1)
    y = dataset['label'].values

    features_names = X.columns.tolist()

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42, stratify=y
    )

    return X_train, X_test, y_train, y_test, le, features_names


def undersample_data(X_train, y_train):
    rs = RandomUnderSampler(random_state=42)
    X_train, y_train = rs.fit_resample(X_train, y_train)

    print('Resampled dataset shape %s' % Counter(y_train))
    logging.info('Resampled dataset shape %s' % Counter(y_train))

    return X_train, y_train

def oversample_data(X_train, y_train):
    sm = SMOTE(random_state=42)
    X_train, y_train = sm.fit_resample(X_train, y_train)

    print('Resampled dataset shape %s' % Counter(y_train))
    logging.info('Resampled dataset shape %s' % Counter(y_train))

    return X_train, y_train

def tune_and_evaluate(name, model, params, sampling, 
    X_train, y_train, X_test, y_test, le, y_test_bin,
    feature_names, output_file=None, output_dir=".",
    cv=3, n_jobs=4, n_iter=10):
    
    pipe = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", model)
    ])

    grid = RandomizedSearchCV(
        pipe, params, n_iter=n_iter, cv=cv, n_jobs=n_jobs, random_state=42
    )
    grid.fit(X_train, y_train)
    best_model = grid.best_estimator_

    print(f"Best params for {name}: {grid.best_params_}")

    y_pred = best_model.predict(X_test)
    report = classification_report(y_test, y_pred, output_dict=True)
    y_pred_bin = label_binarize(y_pred, classes=list(range(len(le.classes_))))

    auc_roc = roc_auc_score(y_test_bin, y_pred_bin, average="macro", multi_class="ovr")
    auc_pr = average_precision_score(y_test_bin, y_pred_bin, average="macro")

    result = {
        "Model": name,
        "Precision": report["weighted avg"]["precision"],
        "Recall": report["weighted avg"]["recall"],
        "F1-score": report["weighted avg"]["f1-score"],
        "AUC-ROC": auc_roc,
        "AUC-PRC": auc_pr
    }

    if output_file is None:
        os.makedirs(output_dir, exist_ok=True)
        safe_name = name.replace(" ", "_")
        output_file = os.path.join(output_dir, f"model_{sampling}_{safe_name}.pkl")

    joblib.dump({
        "pipeline": best_model,
        "label_encoder": le,
        "results": result,
        "feature_names": feature_names
    }, output_file)

    print(f"Saved model to {output_file}")
    return result, best_model


def main(args):
    data = import_data_everything(args.dataset_file)
    X_train, X_test, y_train, y_test, le, features_names = import_split_scale_data(data)
    y_test_bin = label_binarize(y_test, classes=list(range(len(le.classes_))))

    all_tasks = {
        "RF": (
            RandomForestClassifier(random_state=42, n_jobs=args.n_jobs),
            {
                'clf__n_estimators': [100, 200, 300],
                'clf__max_depth': [None, 10, 20, 30],
                'clf__min_samples_split': [2, 5, 10],
                'clf__max_features': ['sqrt', 'log2']
            }
        ),
        "XGBoost": (
            XGBClassifier(eval_metric='mlogloss', random_state=42, n_jobs=args.n_jobs),
            {
                'clf__n_estimators': [100, 200],
                'clf__max_depth': [3, 5, 7],
                'clf__learning_rate': [0.01, 0.1, 0.2],
                'clf__subsample': [0.7, 0.8, 1.0],
                'clf__colsample_bytree': [0.7, 0.8, 1.0]
            }
        )
    }

    if args.oversampling:
        X_train_, y_train_ = oversample_data(X_train, y_train)
        sampling_type = "oversampling"
    else:
        X_train_, y_train_ = undersample_data(X_train, y_train)
        sampling_type = "undersampling"

    selected_tasks = []
    for model in args.models:
        if model not in all_tasks:
            raise ValueError(f"Invalid model name: {model}. Choose from {list(all_tasks.keys())}.")
        selected_tasks.append((model, *all_tasks[model]))

    if args.output_files:
        if len(args.output_files) != len(selected_tasks):
            raise ValueError("Number of output files must match number of selected models.")
        output_files = args.output_files
    else:
        os.makedirs(args.output_dir, exist_ok=True)
        output_files = [
            os.path.join(args.output_dir, f"model_{sampling_type}_{name}.pkl")
            for name, _, _ in selected_tasks
        ]

    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        for (name, model, params), output_file in zip(selected_tasks, output_files):
            if not output_file.endswith(".pkl"):
                raise ValueError(f"Output file '{output_file}' must have a .pkl extension.")
            executor.submit(
                tune_and_evaluate,
                name, model, params, sampling_type,
                X_train_, y_train_, X_test, y_test, le, y_test_bin, features_names,
                output_file,
                args.output_dir,
                args.cv,
                args.n_jobs,
                args.n_iter
            )

if __name__ == "__main__":
    plant_dir = os.environ.get("PLANT_DIR", ".")
    OUTPUT_DIR_NAME = os.path.join(plant_dir, "models", "trained_models")

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dataset_file", required=True)
    parser.add_argument("-m", "--models", nargs="+", choices=["RF", "XGBoost"], default=["RF", "XGBoost"],
                        help="Which models to train: RF, XGBoost, or both")
    
    parser.add_argument("-o", "--output_files", nargs="+", default=None,
                        help="Optional .pkl filenames (one per model). If omitted, defaults are auto-generated.")
    parser.add_argument("--output_dir", default=OUTPUT_DIR_NAME,
                        help="Directory where models will be saved (default: models/trained_models)")
    
    parser.add_argument("-ov", "--oversampling", action="store_true", default=False,
                        help="Use oversampling (default: oversampling)")

    parser.add_argument("--cv", type=int, default=3, help="Number of CV folds for RandomizedSearchCV")
    parser.add_argument("--n_jobs", type=int, default=4, help="Number of parallel jobs for training")
    parser.add_argument("--max_workers", type=int, default=8, help="Max workers for parallel model training")
    parser.add_argument("--n_iter", type=int, default=10, help="Number of iterations for RandomizedSearchCV")

    args = parser.parse_args()
    main(args)
