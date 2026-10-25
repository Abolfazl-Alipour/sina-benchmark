"""Phase 4: nested model selection + HPO on v5-all features, chr folds.

Pre-registered protocol (fixed before seeing any HPO number):
- Outer: 5 chromosome-disjoint folds (gene folds consistently worse; chr is
  the primary scheme. Documented, not hidden).
- Inner: 3-fold random split of each outer train pool, seed 1337.
- Objective: mean inner AUPRC. Budget: 50 Optuna trials each for LightGBM
  and XGBoost (fixed grids below, recorded in output); L1-logistic C-path
  (20 values) as linear floor.
- Winner per outer fold retrained on the full outer train pool, scored once
  on the outer test fold. OOF scores pooled for covered PR AUC diff vs OT.
- Samplers seeded; n_jobs=1 (deterministic).
"""

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts" / "m4b"))
sys.path.insert(0, str(REPO_ROOT / "m5b" / "scripts"))
from run_v5 import BASE, load_matrix, prauc_diff_ci  # noqa: E402
from tabular_baseline import (  # noqa: E402
    auprc,
    chr_folds,
    paired_top1,
)

M5B_DIR = REPO_ROOT / "m5b"
SEED = 1337
N_TRIALS = 50

LGBM_SPACE = {"num_leaves": (15, 63), "learning_rate": (0.01, 0.1),
              "min_child_samples": (10, 100), "feature_fraction": (0.6, 1.0),
              "bagging_fraction": (0.6, 1.0)}
XGB_SPACE = {"max_depth": (3, 8), "eta": (0.01, 0.1),
             "subsample": (0.6, 1.0), "colsample_bytree": (0.6, 1.0),
             "min_child_weight": (1, 10)}


def train_lgbm_params(tr, va, cols, seed, **kw):
    import lightgbm as lgb

    pos = float(tr["label"].sum())
    params = {"objective": "binary", "metric": "average_precision",
              "seed": seed, "deterministic": True, "num_threads": 4,
              "verbosity": -1, "bagging_freq": 1,
              "scale_pos_weight": float(len(tr) - pos) / max(pos, 1.0)}
    params.update(kw)
    dtr = lgb.Dataset(tr[cols], label=tr["label"])
    dva = lgb.Dataset(va[cols], label=va["label"], reference=dtr)
    model = lgb.train(params, dtr, num_boost_round=2000,
                      valid_sets=[dva],
                      callbacks=[lgb.early_stopping(100, verbose=False)])
    return model


def train_xgb_params(tr, va, cols, seed, **kw):
    import xgboost as xgb

    pos = float(tr["label"].sum())
    params = {"objective": "binary:logistic", "eval_metric": "aucpr",
              "seed": seed, "nthread": 4, "verbosity": 0,
              "scale_pos_weight": float(len(tr) - pos) / max(pos, 1.0)}
    params.update(kw)
    dtr = xgb.DMatrix(tr[cols], label=tr["label"])
    dva = xgb.DMatrix(va[cols], label=va["label"])
    return xgb.train(params, dtr, num_boost_round=2000,
                     evals=[(dva, "va")], early_stopping_rounds=100,
                     verbose_eval=False)


def train_logistic(tr, cols, seed, c_val):
    from sklearn.linear_model import LogisticRegression

    model = LogisticRegression(penalty="l1", solver="liblinear", C=c_val,
                               class_weight="balanced",
                               random_state=seed, max_iter=2000)
    model.fit(tr[cols].fillna(0.0), tr["label"])
    return model


def predict(model, df, cols, kind):
    if kind == "xgb":
        import xgboost as xgb

        return model.predict(xgb.DMatrix(df[cols]))
    if kind == "logistic":
        return model.predict_proba(df[cols].fillna(0.0))[:, 1]
    return model.predict(df[cols])


def inner_splits(n, seed=SEED):
    import numpy as np

    rng = np.random.default_rng(seed)
    idx = rng.permutation(n)
    return [idx[i::3] for i in range(3)]


def objective_for(kind, df_pool, cols, seed):
    import numpy as np

    pos = np.arange(len(df_pool))
    folds = inner_splits(len(df_pool), seed)

    def objective(trial):
        ap = []
        for k in range(3):
            va_idx = folds[k]
            tr_idx = np.concatenate([folds[j] for j in range(3) if j != k])
            tr = df_pool.iloc[tr_idx]
            va = df_pool.iloc[va_idx]
            cut = int(len(tr) * 0.9)
            if kind == "lgbm":
                kw = {"num_leaves": trial.suggest_int("num_leaves", 15, 63),
                      "learning_rate": trial.suggest_float("learning_rate",
                                                           0.01, 0.1),
                      "min_child_samples": trial.suggest_int(
                          "min_child_samples", 10, 100),
                      "feature_fraction": trial.suggest_float(
                          "feature_fraction", 0.6, 1.0),
                      "bagging_fraction": trial.suggest_float(
                          "bagging_fraction", 0.6, 1.0)}
                m = train_lgbm_params(tr.iloc[:cut], tr.iloc[cut:], cols,
                                      seed, **kw)
            elif kind == "xgb":
                kw = {"max_depth": trial.suggest_int("max_depth", 3, 8),
                      "eta": trial.suggest_float("eta", 0.01, 0.1),
                      "subsample": trial.suggest_float("subsample", 0.6, 1.0),
                      "colsample_bytree": trial.suggest_float(
                          "colsample_bytree", 0.6, 1.0),
                      "min_child_weight": trial.suggest_int(
                          "min_child_weight", 1, 10)}
                m = train_xgb_params(tr.iloc[:cut], tr.iloc[cut:], cols,
                                     seed, **kw)
            else:
                c_v = trial.suggest_float("C", 1e-3, 1e2, log=True)
                m = train_logistic(tr, cols, seed, c_v)
            s = predict(m, va, cols, kind)
            ap.append(auprc(va["label"].to_numpy(),
                            np.nan_to_num(np.asarray(s), nan=-1e18)))
        return float(np.mean(ap))

    _ = pos
    return objective


def main() -> int:
    import numpy as np
    import optuna

    optuna.logging.set_verbosity(optuna.logging.WARNING)
    df, ot_scores, avail = load_matrix()
    cols = BASE + [c for t in avail.values() for c in t]
    print(f"hpo cols ({len(cols)}): {cols}", flush=True)
    folds = chr_folds(df)
    oof = np.zeros(len(df))
    winners = []
    for fold in range(5):
        te = folds == fold
        pool = df[~te].sample(frac=1.0, random_state=SEED).reset_index(
            drop=True)
        best = ("none", -1.0, {})
        for kind in ("lgbm", "xgb", "logistic"):
            trials = 20 if kind == "logistic" else N_TRIALS
            study = optuna.create_study(
                direction="maximize",
                sampler=optuna.samplers.TPESampler(seed=SEED))
            study.optimize(objective_for(kind, pool, cols, SEED),
                           n_trials=trials, n_jobs=1)
            if study.best_value > best[1]:
                best = (kind, study.best_value, study.best_params)
        kind, val, params = best
        cut = int(len(pool) * 0.9)
        if kind == "logistic":
            model = train_logistic(pool, cols, SEED, params["C"])
        elif kind == "xgb":
            model = train_xgb_params(pool.iloc[:cut], pool.iloc[cut:],
                                     cols, SEED, **params)
        else:
            model = train_lgbm_params(pool.iloc[:cut], pool.iloc[cut:],
                                      cols, SEED, **params)
        oof[te] = predict(model, df[te], cols, kind)
        winners.append({"fold": fold, "kind": kind,
                        "inner_auprc": round(float(val), 4), "params": params})
        print(f"fold {fold}: {kind} inner {val:.4f} {params}", flush=True)
    df = df.assign(hpo=oof)
    df["negdist"] = -df["dist_kb"]
    cov = df[df["locus_id"].isin(set(ot_scores))].copy()
    from tabular_baseline import summarize_arm

    out = {"cols": cols, "winners": winners,
           "covered_n_loci": int(cov["locus_id"].nunique()),
           "covered_hpo": summarize_arm(cov.assign(gbm=cov["hpo"]), "gbm"),
           "covered_ot": summarize_arm(cov.assign(ot=cov["ot"]), "ot"),
           "prauc_diff_hpo_minus_ot": prauc_diff_ci(
               cov.assign(gbm=cov["hpo"], ot=cov["ot"]), "gbm", "ot"),
           "top1_diff_hpo_minus_ot": paired_top1(
               cov.assign(gbm=cov["hpo"], ot=cov["ot"]), "gbm", "ot")}
    (M5B_DIR / "reports" / "hpo_results.json").write_text(
        json.dumps(out, indent=1), encoding="utf-8")
    d = out["prauc_diff_hpo_minus_ot"]
    print(f"HPO covered PR AUC diff {d['mean']} {d['ci']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
