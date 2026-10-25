"""Track F fetch: per-variant posteriorProbability rows for every cached
OT credible set. Cache-first + resume-safe into m5b/cache/finemap/."""

import json
import time
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC = REPO_ROOT / "data" / "processed" / "m4b" / "ot_cache" / "credible_sets"
DST = REPO_ROOT / "m5b" / "cache" / "finemap"
API = "https://api.platform.opentargets.org/api/v4/graphql"
PAGE = 100

Q = (f"query($slid:String!,$idx:Int!){{credibleSet(studyLocusId:$slid){{"
     f"studyLocusId locus(page:{{index:$idx,size:{PAGE}}})"
     "{count rows{posteriorProbability "
     "variant{id chromosome position}}}}}")


def post(payload: dict) -> dict:
    last = None
    for attempt in range(4):
        try:
            req = urllib.request.Request(
                API, data=json.dumps(payload).encode(),
                headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=90) as resp:
                return json.loads(resp.read().decode())
        except Exception as exc:
            last = exc
            time.sleep(2 + attempt * 2)
    raise RuntimeError(f"GraphQL failed after retries: {last}")


def fetch_set(slid: str) -> dict:
    rows, idx, total = [], 0, None
    while True:
        d = post({"query": Q, "variables": {"slid": slid, "idx": idx}})
        try:
            cs = d["data"]["credibleSet"]
        except (KeyError, TypeError) as err:
            raise RuntimeError(f"bad locus response for {slid}: "
                               f"{json.dumps(d)[:300]}") from err
        if cs is None:
            return {"studyLocusId": slid, "count": 0, "rows": [],
                    "null_upstream": True}
        loc = cs["locus"]
        total = loc["count"]
        rows.extend(loc["rows"])
        if len(rows) >= total or not loc["rows"]:
            break
        idx += 1
    assert len(rows) == total, f"{slid}: {len(rows)} != {total}"
    return {"studyLocusId": slid, "count": total, "rows": rows}


def main() -> int:
    DST.mkdir(parents=True, exist_ok=True)
    sids = sorted(p.stem for p in SRC.glob("*.json"))
    done, pulled = 0, 0
    for slid in sids:
        fpath = DST / f"{slid}.json"
        if fpath.exists():
            done += 1
            continue
        fpath.write_text(json.dumps(fetch_set(slid)) + "\n", encoding="utf-8")
        pulled += 1
        done += 1
        time.sleep(0.2)
        if done % 50 == 0:
            print(f"finemap: {done}/{len(sids)}", flush=True)
    mani = {"n_sets": len(sids), "pulled_this_run": pulled,
            "rows_total": sum(json.loads(
                (DST / f"{s}.json").read_text(
                    encoding="utf-8"))["count"] for s in sids)}
    (DST / "manifest.json").write_text(json.dumps(mani, indent=1),
                                       encoding="utf-8")
    print(f"finemap done: {len(sids)} sets, {mani['rows_total']} rows",
          flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
