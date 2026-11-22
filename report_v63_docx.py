"""Map the 45-table template to saved M01, S01 and external analysis results."""
import argparse
import csv
import hashlib
import json
import re
import sys
from pathlib import Path
from zipfile import ZipFile

from lxml import etree as E
sys.dont_write_bytecode = True
from report_v51_docx import NS, grid, put, text


GROUPS = ("unchanged_untreated", "newly_intervention", "intervention_removed",
          "treatment_intensified", "unchanged_intensity", "treatment_deintensified")
STRATEGIES = ("S0", "SN", "ST", "SNT")
EXTERNAL_MATRICES = {38, 39, 42, 43}
EXTERNAL_REDUCTIONS = {40, 41, 44, 45}


def external_design(table):
    cohort = "S03_UKHLS" if table % 2 == 0 else "S04_ELSA"
    return cohort, "strict" if table <= 41 else "expanded"


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def select(rows, **keys):
    matches = [row for row in rows if all(row[k] == str(v) for k, v in keys.items())]
    if len(matches) != 1:
        raise ValueError(f"Expected one source row for {keys}, found {len(matches)}")
    return matches[0]


def count(n, pct):
    return f"{float(n):,.0f} ({float(pct):.1f})"


def summary(row):
    if row["n_positive"]:
        return count(row["n_positive"], row["positive_pct"])
    return f'{float(row["median"]):.1f} ({float(row["q25"]):.1f}–{float(row["q75"]):.1f})'


class Report:
    def __init__(self, analysis, intervention, legacy, blueprint):
        self.analysis = analysis
        self.intervention = intervention
        self.legacy = legacy
        self.blueprint = read_csv(blueprint)
        self.cache = {}
        self.used = set()
        self.hashes = {}
        self.verify_legacy()

    def verify_legacy(self):
        manifest = self.legacy / "quality_control/source_manifest.csv"
        for row in self.read(manifest):
            path = Path(row["path"])
            if not path.is_file() or hashlib.md5(path.read_bytes()).hexdigest() != row["md5"]:
                raise ValueError(f"Stale v6.2 statistics. Refresh M01 step 13 legacy sources: {path}")

    def read(self, path):
        self.used.add(path)
        if path not in self.cache:
            self.cache[path] = read_csv(path)
        return self.cache[path]

    def legacy_cell(self, n, r, c, labels):
        old = {16: 17, 17: 16, **{i: i-1 for i in range(19, 23)},
               **{i: i-2 for i in range(24, 30)}}.get(n, n)
        if n == 12:
            # The sample has two identical baseline LDL-C columns.
            c = {3: 2, 4: 5, 5: 6}.get(c, c)
        if n >= 24 and n <= 28:
            if n == 26:
                r -= 1
            col = int(n == 24)
            label = labels[col]
            for new, prior in (("HMR-1", "HMR-A"), ("HMR-2", "HMR-B"),
                               ("STRICT-1", "STRICT"), ("STRICT-2", "REMOVE2"),
                               ("STRICT-3", "DEP-P90"), ("STRICT-4", "CRP-2–10")):
                if label.startswith(new):
                    label = prior + label[len(new):]
                    break
            prior = select(self.blueprint, table=old, row=r, col=col)["text"]
            if prior != label:
                # The labels for the hs-CRP threshold differ between releases.
                if not (labels[col].startswith("STRICT-4") and
                        prior.endswith(labels[col].split(" ", 1)[1]) and "CRP" in prior):
                    raise ValueError(f"Scenario mismatch S{n}/{r}: {prior!r}, {label!r}")
        path = self.legacy / "tables" / f"Supplementary_Table_{old:02}_cells.csv"
        item = select(self.read(path), row=r, col=c)
        return item["value"]

    def external(self, cohort, filename):
        return self.read(self.intervention / cohort / "external_RM_analyses/03_pooled" / filename)

    def matrix(self, n, r, c):
        panel, row = divmod(r-2, 8)
        strategy = STRATEGIES[panel+1]
        if n == 18:
            rows = self.read(self.analysis / "treatment_transition_matrix_S0_vs_RM/12A_matrix_cells_MI_mean.csv")
            selected = [x for x in rows if x["scenario"] == strategy and
                        (row == 6 or x["row_order"] == str(row+1)) and
                        (c == 7 or x["column_order"] == str(c))]
            field = "mean_n"
        else:
            cohort, rm_set = external_design(n)
            rows = self.external(cohort, "05_six_group_matrix_MI_mean.csv")
            selected = [x for x in rows if x["RM_set"] == rm_set and x["strategy"] == strategy and
                        (row == 6 or x["score_based_group"] == GROUPS[row]) and
                        (c == 7 or x["RM_based_group"] == GROUPS[c-1])]
            field = "n"
        if len(selected) != (6 if row == 6 else 1) * (6 if c == 7 else 1):
            raise ValueError(f"Incomplete six-group matrix S{n}/{r}/{c}")
        return count(sum(float(x[field]) for x in selected),
                     sum(float(x["pct_total"]) for x in selected))

    def fixed_contribution(self, r, c):
        rows = self.read(self.analysis / "fixed_S0_treatment_group_RD/12B_pooled_RD_effects.csv")
        group = "Overall" if r == 7 else f"G{r}"
        matches = [x for x in rows if x["fixed_group"].lower() == group.lower()]
        strategy = STRATEGIES[c-1 if c <= 4 else c-4]
        estimand = "weighted_rd_contribution" if c <= 4 else "weighted_delta_rd_contribution"
        x = select(matches, scenario=strategy, time=15, estimand=estimand)
        if x["status"] != "estimable":
            raise ValueError(f"Unusable fixed group contribution: {x}")
        a, lo, hi = (100 * float(x[k]) for k in ("estimate", "lcl", "ucl"))
        return f"{a:.2f} ({lo:.2f} to {hi:.2f})"

    def external_cell(self, n, r, c, labels):
        if n in (32, 33):
            cohort = "S03_UKHLS" if c <= (2 if n == 32 else 1) else "S04_ELSA"
        elif n == 34:
            cohort = "S03_UKHLS" if r < 12 else "S04_ELSA"
        elif n in (35, 36, 37):
            cohort = "S03_UKHLS" if c <= (3 if n == 35 else 2) else "S04_ELSA"
        else:
            cohort, rm_set = external_design(n)
        elsa = cohort == "S04_ELSA"
        if n == 32:
            rows = self.external(cohort, "S32_selection_flow.csv")
            x = select(rows, step=min(8, max(1, r-1)))
            field = "n_before" if r == 1 else ("n_excluded" if c % 2 else "n_remaining")
            return f"{int(x[field]):,}"
        if n == 33:
            rows = self.external(cohort, "S33_baseline.csv")
            fields = ("age", "age", "female", "current_smoking", "sbp", "dbp",
                      "tc", "ldl", "hdl", "trig" if elsa else "tg", "moderate_2019", "moderate_2025")
            x = select(rows, variable=fields[r])
            if r == 0:
                return f'{"ELSA" if elsa else "UKHLS"} (N = {int(x["n_total"]):,})'
            return summary(x)
        if n == 34:
            uk = ("rm_south_asian", "bmi", "waist_cm", "rm_obesity_central", "hscrp",
                  "rm_hscrp_gt2", "rm_social_deprivation", "rm_physical_inactivity_sport_proxy",
                  "rm_chronic_inflammatory_arthritis_proxy")
            el = ("bmi", "waist", "rm_obesity_central", "hscrp_w6", "rm_hscrp_gt2",
                  "rm_social_deprivation_proxy", "rm_physical_inactivity_proxy", "rm_rheumatoid_proxy",
                  "rm_major_psychiatric_proxy", "rm_natural_menopause_lt40_proxy")
            variable = el[r-14] if elsa else uk[r-2]
            rows = self.external(cohort, "S34_imputation.csv")
            x = select(rows, variable=variable, stage="before" if c <= 2 else "after")
            if c == 2:
                return count(x["n_missing"], x["missing_pct"])
            if c == 4:
                return f'{float(x["n_missing"]):,.0f}'
            return summary(x)
        if n == 35:
            modifier = {2:"obesity", 3:"hscrp", 4:"SN", 5:"ST", 7:"deprivation",
                        8:"psychiatric", 9:"south_asian", 10:"obesity", 11:"inactivity",
                        12:"arthritis", 13:"menopause", 14:"hscrp", 15:"SN", 16:"ST", 17:"SNT"}[r]
            x = select(self.external(cohort, "S35_RM_prevalence.csv"), modifier=modifier,
                       RM_set="strict" if r < 6 else "expanded", guideline=2019 if c % 2 == 0 else 2025)
            return count(x["n"], x["pct"])
        if n in (36, 37):
            year = 2019 if c % 2 else 2025
            levels = self.external(cohort, "01_risk_levels_MI_mean.csv")
            if n == 36 and r < 6:
                x = select(levels, RM_set="expanded", strategy="S0", guideline=year,
                           category=("L", "M", "H", "VH")[r-2])
                return count(x["n"], x["pct_total"])
            strategy = re.search(r"\b(S0|SN|ST|SNT)\b", labels[0]).group()
            rm_set = "strict" if (r < 10 if n == 36 else 3 < r < 7) else "expanded"
            x = select(self.external(cohort, "02_hvh_eligibility_MI_mean.csv"),
                       RM_set=rm_set, strategy=strategy, guideline=year)
            if n == 37:
                return count(x["n_eligible"], x["pct_eligible"])
            base = select(levels, RM_set=rm_set, strategy="S0", guideline=year, category="M")
            return count(x["n_upgraded"], 100*float(x["n_upgraded"])/float(base["n"]))
        if n not in EXTERNAL_REDUCTIONS:
            raise ValueError(f"No source mapping for table {n}")
        panel, offset = divmod(r-2, 8)
        year = (2019, 2025)[panel]
        strategy = STRATEGIES[c-1]
        if offset == 6:
            if year == 2019:
                strategy = {"ST": "S0", "SNT": "SN"}.get(strategy, strategy)
            x = select(self.external(cohort, "02_hvh_eligibility_MI_mean.csv"),
                       RM_set=rm_set, guideline=year, strategy=strategy)
            return f'{float(x["mean_minimum_D"]):.2f}'
        x = select(self.external(cohort, "03_six_transitions_MI_mean.csv"),
                   RM_set=rm_set, strategy=strategy, transition=GROUPS[offset])
        return f'{float(x[f"mean_D_{year}"]):.2f}'

    def value(self, n, r, c, labels):
        if n == 18 or n in EXTERNAL_MATRICES:
            return self.matrix(n, r, c)
        if n == 23:
            return self.fixed_contribution(r, c)
        if n >= 32:
            return self.external_cell(n, r, c, labels)
        return self.legacy_cell(n, r, c, labels)


def build(args):
    spec = json.loads(args.lock.read_text(encoding="utf-8-sig"))
    if hashlib.sha256(args.template.read_bytes()).hexdigest() != spec["sha256"]:
        print("Template SHA-256 differs from the lock; checking the actual cell layout.")
    report = Report(args.analysis, args.intervention, args.legacy, args.blueprint)
    audit = []
    with ZipFile(args.template) as source:
        doc = E.fromstring(source.read("word/document.xml"))
        tables = doc.findall("w:body/w:tbl", NS)
        if len(tables) != len(spec["tables"]) or len(tables) != 45:
            raise ValueError("Expected 45 tables")
        for n, (table, info) in enumerate(zip(tables, spec["tables"]), 1):
            rows = table.findall("w:tr", NS)
            if len(rows) != info["rows"]:
                raise ValueError(f"Row count mismatch S{n}")
            for r, (row, row_spec) in enumerate(zip(rows, info["rowdata"])):
                cells = grid(row)
                specs = row_spec["cells"]
                labels = [s["text"] for s in specs]
                if len(cells) != len(specs):
                    raise ValueError(f"Grid mismatch S{n}/{r}: {len(cells)} versus {len(specs)}")
                for c, slot in enumerate(specs):
                    if c == 0 and text(cells[c]).strip() and slot["text"] and not slot.get("dynamic", False) and text(cells[c]).strip() != slot["text"].strip():
                        raise ValueError(f"Row label mismatch S{n}/{r}: {text(cells[c])!r} versus {slot['text']!r}")

                    if not slot.get("dynamic", False):
                        continue
                    where = f"S{n}/{r}/{c}"
                    try:
                        cell = cells[c]
                        report.used.clear()
                        value = report.value(n, r, c, labels)
                        if not value or re.search(r"XX|\[n\]|TBC", value):
                            raise ValueError(f"Missing value: {value!r}")
                        if c and cell is cells[c-1] and text(cell) != value:
                            raise ValueError("Merged slots require different values")
                        put(cell, value)
                        paths = sorted(report.used)
                        for path in paths:
                            if path not in report.hashes:
                                report.hashes[path] = hashlib.sha256(path.read_bytes()).hexdigest()
                        audit.append((n, r, c, value, "|".join(map(str, paths)),
                                      "|".join(report.hashes[path] for path in paths)))
                    except (IndexError, KeyError, ValueError) as exc:
                        raise ValueError(f"{where}: {exc}") from exc
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_suffix(".tmp.docx")
        try:
            with ZipFile(temporary, "w") as target:
                for entry in source.infolist():
                    content = E.tostring(doc, encoding="UTF-8", xml_declaration=True) if entry.filename == "word/document.xml" else source.read(entry)
                    target.writestr(entry, content)
            temporary.replace(args.output)
        finally:
            temporary.unlink(missing_ok=True)
    with args.output.with_suffix(".cells.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("table", "row_zero_based", "column_zero_based", "value", "source_files", "source_sha256"))
        writer.writerows(audit)
    print(f"Wrote {len(audit)} cells to {args.output}")


def build_main_table2(template, intervention, output):
    annual = []
    paired = []
    for cohort in ("S03_UKHLS", "S04_ELSA"):
        sources = intervention / cohort / "external_RM_analyses/03_pooled"
        annual.append(read_csv(sources / "Table_2_guideline_results.csv"))
        paired.append(read_csv(sources / "Table_2_transition_net_changes.csv"))
    with ZipFile(template) as src:
        doc = E.fromstring(src.read("word/document.xml"))
        table, = doc.findall("w:body/w:tbl", NS)
        rows = table.findall("w:tr", NS)
        if len(rows) != 26:
            raise ValueError(f"Table 2 expected 26 rows, found {len(rows)}")
        for section, start in enumerate((3, 11, 19)):
            rm_set = "strict" if section < 2 else "expanded"
            strategy = "S0" if section == 0 else "SNT"
            for offset in range(7):
                cells = grid(rows[start + offset])
                if offset in (0, 1, 6):
                    for cohort in range(2):
                        for year_index, year in enumerate((2019, 2025)):
                            item = select(annual[cohort], RM_set=rm_set, strategy=strategy, guideline=year)
                            if offset == 6:
                                value = (f'{float(item["mean_minimum_D"]):.2f}' if not section
                                         else f'{float(item["delta_mean_D_vs_S0"]):+.2f}')
                            else:
                                suffix = "hvh" if offset == 0 else "eligible"
                                value = count(item[f"n_{suffix}"], item[f"pct_{suffix}"])
                            put(cells[1 + 2*cohort + year_index], value)
                else:
                    for cohort in range(2):
                        item = select(paired[cohort], RM_set=rm_set, transition=GROUPS[(1,2,3,5)[offset-2]])
                        value = (count(item["n_S0"], item["pct_S0"]) if not section else
                                 f'{float(item["delta_n"]):+,.0f} ({float(item["delta_pp"]):+.1f})')
                        put(cells[1+2*cohort], value)
        note, = [p for p in doc.findall("w:body/w:p", NS) if text(p).startswith("Note.")]
        put(note, "Note. The strict RM set includes obesity or central obesity in both frameworks and "
            "elevated high-sensitivity C-reactive protein (hs-CRP >2 mg/L) in the 2025 framework only. "
            "The expanded RM set adds proxies for social deprivation, physical inactivity, and chronic "
            "inflammatory disease (arthritis) in both cohorts. It also adds a proxy for major psychiatric "
            "disorder in ELSA under both frameworks, high-risk South Asian ethnicity in UKHLS under the "
            "2025 framework, and premature menopause in ELSA under the 2025 framework. "
            "S0 applies the 2019 and 2025 guidelines without risk modifiers (RMs). "
            "Each RM set pairs 2019 SN with 2025 SNT. The categories are defined by comparing treatment "
            "decisions between the two frameworks. Category changes are net changes from S0 in the full "
            "analytical cohort, with percentages expressed as percentage-point changes. High/very-high risk "
            "and treatment eligibility are reported as n (%). Mean minimum required LDL-C reduction is "
            "reported in mmol/L under S0 and as the within-guideline change from S0 under each RM set.")
        tmp = output.with_suffix(".tmp.docx")
        try:
            with ZipFile(tmp, "w") as dst:
                for entry in src.infolist():
                    dst.writestr(entry, E.tostring(doc, encoding="UTF-8", xml_declaration=True)
                                 if entry.filename == "word/document.xml" else src.read(entry))
            tmp.replace(output)
        finally:
            tmp.unlink(missing_ok=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("template", "lock", "analysis", "intervention", "legacy", "blueprint", "output"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--table2-template", type=Path)
    args = parser.parse_args()
    build(args)
    if args.table2_template:
        build_main_table2(args.table2_template, args.intervention,
                          args.output.with_name("M01_Main_Table_2_v6.3.docx"))








