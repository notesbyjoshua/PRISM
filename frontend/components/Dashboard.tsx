"use client";

import { useEffect, useMemo, useState } from "react";
import {
  ArrowRight,
  ChartScatter,
  CheckCircle,
  CirclesFour,
  Database,
  Dna,
  Flask,
  MagnifyingGlass,
  PulseIcon,
  SlidersHorizontal,
  Sparkle,
  TrendUp,
} from "@phosphor-icons/react";
import { EffectScatter, SpecificityBars } from "./Charts";
import MetricCard from "./MetricCard";
import { getJson, Result, Specificity, Summary } from "@/lib/api";

const formatNumber = (value: number | null, digits = 2) =>
  value === null || value === undefined ? "—" : value.toFixed(digits);

const formatQ = (value: number | null) => {
  if (value === null || value === undefined) return "—";
  return value < 0.001 ? value.toExponential(1) : value.toFixed(3);
};

export default function Dashboard() {
  const [summary, setSummary] = useState<Summary | null>(null);
  const [results, setResults] = useState<Result[]>([]);
  const [specificity, setSpecificity] = useState<Specificity[]>([]);
  const [diseases, setDiseases] = useState<string[]>([]);
  const [disease, setDisease] = useState("");
  const [search, setSearch] = useState("");
  const [priorityOnly, setPriorityOnly] = useState(false);
  const [showAll, setShowAll] = useState(false);
  const [error, setError] = useState("");

  useEffect(() => {
    Promise.all([
      getJson<Summary>("/api/summary"),
      getJson<string[]>("/api/diseases"),
      getJson<Specificity[]>("/api/specificity?limit=20"),
    ]).then(([summaryData, diseaseData, specificityData]) => {
      setSummary(summaryData);
      setDiseases(diseaseData);
      setSpecificity(specificityData);
    }).catch(() => setError("The PRISM API is not reachable. Start the FastAPI server on port 8000."));
  }, []);

  useEffect(() => {
    const params = new URLSearchParams({ analysis_type: "pairwise", limit: "100", sort_by: "q_value" });
    if (disease) params.set("disease", disease);
    if (search) params.set("search", search);
    if (priorityOnly) params.set("high_priority", "true");
    getJson<{ items: Result[] }>(`/api/results?${params}`).then((data) => setResults(data.items)).catch(() => undefined);
  }, [disease, search, priorityOnly]);

  const visibleResults = useMemo(
    () => showAll ? results : results.slice(0, 12),
    [results, showAll],
  );

  return (
    <div className="app-shell">
      <aside className="sidebar">
        <div className="brand"><span className="brand-mark"><Dna size={22} weight="bold" /></span><span>PRISM</span></div>
        <nav aria-label="Primary navigation">
          <a className="nav-item active" href="#overview"><CirclesFour size={19} />Overview</a>
          <a className="nav-item" href="#explorer"><ChartScatter size={19} />Feature explorer</a>
          <a className="nav-item" href="#specificity"><Sparkle size={19} />Specificity</a>
          <a className="nav-item" href="http://localhost:8000/docs"><Database size={19} />API documentation</a>
        </nav>
        <div className="sidebar-note">
          <Flask size={20} weight="duotone" />
          <div><strong>Research preview</strong><span>For exploratory analysis, not clinical diagnosis.</span></div>
        </div>
      </aside>

      <main>
        <header className="topbar">
          <div><span className="status-dot" />Local analysis workspace</div>
          <button className="avatar" aria-label="User profile">JP</button>
        </header>

        <div className="content">
          <section className="hero" id="overview">
            <div>
              <p className="kicker">Phenotype intelligence</p>
              <h1>See which facial measurements <em>matter.</em></h1>
              <p>Explore effect sizes, reliability, specificity, and discriminative power across rare-disease cohorts.</p>
            </div>
            <div className="hero-orbit" aria-hidden="true"><Dna size={68} weight="duotone" /></div>
          </section>

          {error && <div className="alert">{error}</div>}
          {summary?.is_demo && <div className="demo-banner"><Sparkle size={17} />Showing illustrative demo data. Run the analysis to populate this dashboard with your results.</div>}

          <section className="metric-grid" aria-label="Dataset summary">
            <MetricCard label="Features analyzed" value={(summary?.feature_count ?? 0).toLocaleString()} detail="Geometric phenotypes" icon={PulseIcon} />
            <MetricCard label="Disease cohorts" value={(summary?.disease_count ?? 0).toLocaleString()} detail="Compared with healthy" icon={Dna} tone="violet" />
            <MetricCard label="High priority" value={(summary?.high_priority_count ?? 0).toLocaleString()} detail="FDR + effect threshold" icon={CheckCircle} tone="amber" />
            <MetricCard label="Mean AUROC" value={formatNumber(summary?.average_auc ?? null)} detail="Univariate discrimination" icon={TrendUp} />
          </section>

          <section className="filter-panel" id="explorer">
            <div className="filter-title"><SlidersHorizontal size={19} /><span>Explore comparisons</span></div>
            <label className="search-field"><MagnifyingGlass size={18} /><input value={search} onChange={(event) => setSearch(event.target.value)} placeholder="Search a feature…" /></label>
            <label><span>Disease cohort</span><select value={disease} onChange={(event) => setDisease(event.target.value)}><option value="">All diseases</option>{diseases.map((item) => <option key={item}>{item}</option>)}</select></label>
            <label className="check-field"><input type="checkbox" checked={priorityOnly} onChange={(event) => setPriorityOnly(event.target.checked)} /><span>High priority only</span></label>
          </section>

          <section className="chart-grid">
            <article className="panel">
              <div className="panel-heading"><div><p className="eyebrow">Discriminative signal</p><h2>Effect size vs. AUROC</h2></div><span className="legend"><i />High priority</span></div>
              <EffectScatter results={results} />
            </article>
            <article className="panel" id="specificity">
              <div className="panel-heading"><div><p className="eyebrow">Distinctive markers</p><h2>Feature specificity</h2></div><span className="panel-tag">Top weighted</span></div>
              <SpecificityBars data={specificity} />
            </article>
          </section>

          <section className="panel results-panel">
            <div className="panel-heading"><div><p className="eyebrow">Evidence table</p><h2>Feature comparisons</h2></div><span className="result-count">{results.length} results</span></div>
            <div className="table-wrap">
              <table>
                <thead><tr><th>Feature</th><th>Disease</th><th>Hedges’ g</th><th>95% CI</th><th>AUROC</th><th>FDR q</th><th>Stability</th><th>Status</th></tr></thead>
                <tbody>{visibleResults.map((row) => (
                  <tr key={row.id}>
                    <td><strong>{row.feature.replaceAll("_", " ")}</strong></td>
                    <td>{row.disease}</td>
                    <td className={(row.hedges_g ?? 0) < 0 ? "negative" : "positive"}>{formatNumber(row.hedges_g)}</td>
                    <td>{formatNumber(row.hedges_g_ci_low)} · {formatNumber(row.hedges_g_ci_high)}</td>
                    <td><span className="auc-pill">{formatNumber(row.roc_auc)}</span></td>
                    <td>{formatQ(row.q_value)}</td>
                    <td>{row.rank_stability === null ? "—" : `${Math.round(row.rank_stability * 100)}%`}</td>
                    <td>{row.high_priority ? <span className="priority"><CheckCircle size={14} weight="fill" />Priority</span> : <span className="muted-status">Monitor</span>}</td>
                  </tr>
                ))}</tbody>
              </table>
            </div>
            {results.length > 12 && (
              <button className="text-button" onClick={() => setShowAll((value) => !value)}>
                {showAll ? "Show top results" : "View all comparisons"} <ArrowRight size={16} />
              </button>
            )}
          </section>
        </div>
      </main>
    </div>
  );
}
