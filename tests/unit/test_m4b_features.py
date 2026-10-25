"""Focused tests for m5b Track F overlap math (synthetic only)."""

import importlib.util
import sys
from pathlib import Path

SCRIPT = (Path(__file__).resolve().parents[2] / "m5b" / "scripts" /
          "features_trackf.py")


def _load():
    spec = importlib.util.spec_from_file_location("features_trackf", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


mod = _load()


def test_overlap_picks_max_and_counts():
    f = mod.overlap_features([(1000, 0.5), (2000, 0.05), (10**7, 0.99)],
                             1500)
    assert f["fm_max_pp"] == 0.5
    assert f["fm_sum_pp"] == 0.55
    assert f["fm_n_vars"] == 2


def test_overlap_window_boundary():
    assert mod.overlap_features([(100_000, 0.3)], 0)["fm_n_vars"] == 1
    assert mod.overlap_features([(100_001, 0.3)], 0)["fm_n_vars"] == 0


def test_overlap_empty_is_zero():
    f = mod.overlap_features([], 500)
    assert (f["fm_max_pp"], f["fm_sum_pp"], f["fm_n_vars"]) == (0.0, 0.0, 0)


def _load_trackc():
    script = (Path(__file__).resolve().parents[2] / "m5b" / "scripts" /
              "features_trackc.py")
    spec = importlib.util.spec_from_file_location("features_trackc", script)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _load_tracke():
    script = (Path(__file__).resolve().parents[2] / "m5b" / "scripts" /
              "features_tracke.py")
    spec = importlib.util.spec_from_file_location("features_tracke", script)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _load_trackh():
    script = (Path(__file__).resolve().parents[2] / "m5b" / "scripts" /
              "features_trackh.py")
    spec = importlib.util.spec_from_file_location("features_trackh", script)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _load_hpo():
    script = (Path(__file__).resolve().parents[2] / "m5b" / "scripts" /
              "hpo_select.py")
    spec = importlib.util.spec_from_file_location("hpo_select", script)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_hpo_deterministic_same_seed():
    import pandas as pd

    hpo = _load_hpo()
    rng = pd.DataFrame({"label": [1, 0] * 30,
                        "f1": sorted([float(i % 7) for i in range(60)]),
                        "f2": [float(i % 3) for i in range(60)]})
    import optuna

    optuna.logging.set_verbosity(optuna.logging.WARNING)

    def run():
        study = optuna.create_study(
            direction="maximize",
            sampler=optuna.samplers.TPESampler(seed=1337))
        study.optimize(hpo.objective_for("logistic", rng, ["f1", "f2"],
                                         1337),
                       n_trials=3, n_jobs=1)
        return study.best_params, study.best_value

    assert run() == run()


def _load_v5():
    script = (Path(__file__).resolve().parents[2] / "m5b" / "scripts" /
              "run_v5.py")
    spec = importlib.util.spec_from_file_location("run_v5", script)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_prauc_diff_ci_deterministic_and_sane():
    import pandas as pd

    v5 = _load_v5()
    df = pd.DataFrame({
        "locus_id": ["l1"] * 3 + ["l2"] * 3,
        "label": [1, 0, 0, 0, 1, 0],
        "a": [0.9, 0.1, 0.2, 0.8, 0.9, 0.1],
        "b": [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]})
    r1 = v5.prauc_diff_ci(df, "a", "b", n_boot=50, seed=7)
    r2 = v5.prauc_diff_ci(df, "a", "b", n_boot=50, seed=7)
    assert r1 == r2
    assert r1["n"] == 2
    assert r1["mean"] > 0 and r1["ci"][0] <= r1["mean"] <= r1["ci"][1]
    r0 = v5.prauc_diff_ci(df, "a", "a", n_boot=50, seed=7)
    assert r0["mean"] == 0.0 and r0["ci"] == [0.0, 0.0]


def _load_curate():
    script = (Path(__file__).resolve().parents[2] / "m5b" / "scripts" /
              "curate_prospective.py")
    spec = importlib.util.spec_from_file_location("curate_prospective",
                                                  script)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_cluster_locus_greedy_window():
    cur = _load_curate()
    rows = [(1000, 1e-9, "a"), (1200, 1e-8, "b"), (10**7, 1e-10, "c")]
    loci = cur.cluster_locus(rows)
    assert len(loci) == 2
    assert [m[2] for m in loci[0]] == ["c"]
    assert {m[2] for m in loci[1]} == {"a", "b"}
    assert cur.parse_pval("2E-9") == 2e-9
    assert cur.parse_pval("NR") is None
    assert cur.split_genes("NYAP2 - MIR5702, APOE") == \
        ["NYAP2", "MIR5702", "APOE"]
    assert cur.split_genes(" - ") == []


def test_ccre_overlap_flags_sorted_scan():
    trackh = _load_trackh()
    elems = [(100, 200, False), (500, 600, True)]
    assert trackh.overlap_flags([150, 550, 9999], elems) == (2, True)
    assert trackh.overlap_flags([150], elems) == (1, False)
    assert trackh.overlap_flags([], elems) == (0, False)


def test_eqtl_overlap_stats_exact_key():
    tracke = _load_tracke()
    dis = [(("1", 100), 0.4), (("1", 200), 0.1), (("2", 100), 0.9)]
    eqt = [(("1", 100), "Liver"), (("1", 100), "Blood"), (("9", 9), "Brain")]
    f = tracke.overlap_stats(dis, eqt)
    assert f["eq_best_pp"] == 0.4
    assert f["eq_n_tissues"] == 2
    assert f["eq_any"] == 1.0
    f2 = tracke.overlap_stats(dis, [(("8", 8), "Brain")])
    assert (f2["eq_best_pp"], f2["eq_n_tissues"], f2["eq_any"]) == \
        (0.0, 0, 0.0)


def test_clinvar_gene_features_counts_and_window():
    import pandas as pd

    trackc = _load_trackc()
    patho = pd.DataFrame(
        [("G1", "1", 1000, 1000, 3),
         ("G1", "1", 500_000, 500_000, 1),
         ("G2", "1", 1200, 1200, 2)],
        columns=["symbol", "chr", "start", "stop", "nsub"])
    f = trackc.gene_features(patho, "G1", "1", 1100)
    assert f["cv_gene_count"] == 2
    assert f["cv_gene_maxsub"] == 3
    assert f["cv_win_count"] == 2
    assert f["cv_win_maxsub"] == 3
    f2 = trackc.gene_features(patho, "ZZZ", "2", 1100)
    assert (f2["cv_gene_count"], f2["cv_win_count"]) == (0, 0)
