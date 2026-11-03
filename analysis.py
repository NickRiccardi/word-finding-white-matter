import argparse

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

WFD = "Visual - TOT"
RT = "Visual - Mean RT"
SPEED = "proc_speed_osd"
TRACTS = ["AF_L", "AF_R", "IFOF_L", "IFOF_R", "ILF_L", "ILF_R", "SLF_L", "SLF_R", "UF_L", "UF_R"]
METRICS = {"FA": "dti_fa", "MD": "md", "RD": "rd", "AD": "ad"}
GM_REGIONS = [
    "Opercular inf. frontal gyrus", "Orbital inf. frontal gyrus", "Triangular inf. frontal gyrus",
    "Inf. temporal gyrus", "Middle temporal gyrus", "Sup. temporal gyrus", "Transverse temporal gyrus",
    "Fusiform gyrus", "Temporal pole", "Angular gyrus", "Supramarginal gyrus", "Anterior insula",
    "Posterior insula", "Precentral gyrus",
]
WM = [(tract, metric, f"{tract}_{suffix}") for tract in TRACTS for metric, suffix in METRICS.items()]
GM = [(f"{region} {side}", "volume", f"{region} {side} volume %")
      for region in GM_REGIONS for side in ("left", "right")]
TAILS = {"FA": "less", "MD": "greater", "RD": "greater", "AD": "two-sided", "volume": "less"}
WINSOR_CAP = 2
N_PERM = 5000
SEED = 20260623


def residual_ranks(d, cols, covars):
    ranks = d.rank()
    z = np.column_stack([np.ones(len(d))] + [ranks[c].to_numpy() for c in covars])
    return [ranks[c].to_numpy() - z @ np.linalg.lstsq(z, ranks[c].to_numpy(), rcond=None)[0] for c in cols]


def p_value(t, dof, tail):
    if tail == "less":
        return stats.t.cdf(t, dof)
    if tail == "greater":
        return stats.t.sf(t, dof)
    return 2 * stats.t.sf(abs(t), dof)


def partial_spearman(df, x, y, covars, tail):
    d = df[[x, y] + covars].dropna()
    ex, ey = residual_ranks(d, [x, y], covars)
    rho = np.corrcoef(ex, ey)[0, 1]
    dof = len(d) - 2 - len(covars)
    t = rho * np.sqrt(dof / (1 - rho**2))
    return {"n": len(d), "rho": rho, "p": p_value(t, dof, tail)}


def permutation_test(df, x, y, covars, tail):
    d = df[[x, y] + covars].dropna()
    ex, ey = residual_ranks(d, [x, y], covars)
    ex, ey = ex - ex.mean(), ey - ey.mean()
    denom = np.sqrt((ex**2).sum() * (ey**2).sum())
    rho = ex @ ey / denom
    rng = np.random.default_rng(SEED)
    null = np.array([ex @ rng.permutation(ey) for _ in range(N_PERM)]) / denom
    if tail == "less":
        hits = np.sum(null <= rho)
    elif tail == "greater":
        hits = np.sum(null >= rho)
    else:
        hits = np.sum(np.abs(null) >= abs(rho))
    return {"n": len(d), "rho": rho, "p": (hits + 1) / (N_PERM + 1)}


def run(analysis, df, measures, y, covars, test=partial_spearman, fdr=False, p_digits=4):
    rows = pd.DataFrame([{"analysis": analysis, "region": region, "metric": metric, "tail": TAILS[metric],
                          **test(df, col, y, covars, TAILS[metric])} for region, metric, col in measures])
    if fdr:
        rows["q"] = multipletests(rows["p"], method="fdr_bh")[1]
    return rows.round({"rho": 3, "p": p_digits, "q": p_digits})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv")
    parser.add_argument("output_csv")
    args = parser.parse_args()

    columns = ["Age", "sex", "education", WFD, RT, SPEED] + [col for _, _, col in WM + GM]
    df = pd.read_csv(args.input_csv, usecols=columns)
    df["male"] = (df["sex"] == "Male").astype(int)
    df["wfd_winsorized"] = df[WFD].clip(upper=WINSOR_CAP)
    rt_recorded = df[df[RT] > 0]
    low_scorers = df[df[WFD] <= WINSOR_CAP]

    age_wfd = pd.DataFrame([{"analysis": "age_x_wfd", "region": "", "metric": "", "tail": "two-sided",
                             **partial_spearman(df, "Age", WFD, [], "two-sided")}]).round({"rho": 3, "p": 3})

    wm_wfd = run("white_matter_x_wfd", df, WM, WFD, ["Age"], fdr=True)
    survivors = [m for m, q in zip(WM, wm_wfd["q"]) if q < 0.05]

    results = pd.concat([
        age_wfd,
        run("age_x_grey_matter", df, GM, "Age", [], p_digits=3),
        run("age_x_white_matter", df, WM, "Age", []),
        wm_wfd,
        run("grey_matter_x_wfd", df, GM, WFD, ["Age"], fdr=True, p_digits=3),
        run("white_matter_x_wfd_full_covariates", df, WM, WFD, ["Age", SPEED, "male", "education"]),
        run("white_matter_x_wfd_naming_rt", rt_recorded, WM, WFD, ["Age", RT]),
        run("white_matter_x_naming_rt", rt_recorded, survivors, RT, ["Age"]),
        run("white_matter_x_wfd_permutation", df, survivors, WFD, ["Age"], permutation_test),
        run("white_matter_x_wfd_winsorized", df, survivors, "wfd_winsorized", ["Age"]),
        run("white_matter_x_wfd_high_scorers_removed", low_scorers, survivors, WFD, ["Age"]),
    ])
    results.to_csv(args.output_csv, index=False)


if __name__ == "__main__":
    main()
