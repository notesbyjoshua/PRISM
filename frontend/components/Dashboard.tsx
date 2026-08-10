"use client";

import { useEffect, useMemo, useRef, useState } from "react";
import {
  ArrowRight, ChartScatter, CheckCircle, CirclesFour, Database, Dna, Flask,
  Images, Info, MagnifyingGlass, PulseIcon, SlidersHorizontal, Sparkle, TrendDown,
  TrendUp, UsersThree, X,
} from "@phosphor-icons/react";
import { EffectScatter, SpecificityBars } from "./Charts";
import MetricCard from "./MetricCard";
import { getJson, Result, Specificity, Summary } from "@/lib/api";
import { displayDisease, displayFeature, featureRegion, REGIONS, Region } from "@/lib/labels";

const formatNumber = (value: number | null, digits = 2) => value == null ? "—" : value.toFixed(digits);
const formatQ = (value: number | null) => value == null ? "—" : value < .001 ? "<0.001" : value.toPrecision(2);

type Distribution = { feature: string; disease: string; healthy: number[]; disease_values: number[] };

function median(values: number[]) {
  if (!values.length) return null;
  const sorted = [...values].sort((a, b) => a - b);
  const middle = Math.floor(sorted.length / 2);
  return sorted.length % 2 ? sorted[middle] : (sorted[middle - 1] + sorted[middle]) / 2;
}

function mean(values: number[]) {
  return values.length ? values.reduce((sum, value) => sum + value, 0) / values.length : null;
}

function Histograms({ data }: { data: Distribution }) {
  const all = [...data.healthy, ...data.disease_values];
  const low = Math.min(...all), high = Math.max(...all), bins = 18;
  const width = high - low || 1;
  const histogram = (values: number[]) => {
    const counts = Array(bins).fill(0) as number[];
    values.forEach((value) => counts[Math.min(bins - 1, Math.floor(((value - low) / width) * bins))]++);
    // Relative frequency keeps groups with different sample sizes comparable.
    return counts.map((count) => count / Math.max(values.length, 1));
  };
  const healthy = histogram(data.healthy), affected = histogram(data.disease_values);
  const sharedMaximum = Math.max(...healthy, ...affected, .01);
  return <div className="distribution-grid">
    {[["Healthy controls", healthy, "healthy"], [displayDisease(data.disease), affected, "affected"]].map(([label, values, tone]) =>
      <div className="distribution-group" key={String(label)}>
        <div className="histogram" aria-label={`${label} relative-frequency histogram`}>
          <span className="histogram-scale"><b>{(sharedMaximum * 100).toFixed(0)}%</b><b>0%</b></span>
          <span className="histogram-midline" aria-hidden="true" />
          {(values as number[]).map((value, index) => <i key={index} className={String(tone)} style={{ height: `${Math.max(2, value / sharedMaximum * 100)}%` }} title={`${(value * 100).toFixed(1)}%`} />)}
        </div>
        <strong>{label}</strong>
        <span>Mean {formatNumber(mean(String(tone) === "healthy" ? data.healthy : data.disease_values), 3)} · Median {formatNumber(median(String(tone) === "healthy" ? data.healthy : data.disease_values), 3)} · n = {String(tone) === "healthy" ? data.healthy.length : data.disease_values.length}</span>
      </div>)}
    <div className="axis-labels"><span>{low.toFixed(3)}</span><span>Measurement value · shared y-axis shows relative frequency</span><span>{high.toFixed(3)}</span></div>
  </div>;
}

function DistributionModal({ row, onClose }: { row: Result; onClose: () => void }) {
  const [data, setData] = useState<Distribution | null>(null);
  const [error, setError] = useState("");
  useEffect(() => {
    const params = new URLSearchParams({ feature: row.feature, disease: row.disease ?? "" });
    getJson<Distribution>(`/api/distributions?${params}`).then(setData).catch(() => setError("The underlying measurements could not be loaded."));
  }, [row]);
  return <div className="modal-backdrop" role="presentation" onMouseDown={onClose}>
    <section className="modal" role="dialog" aria-modal="true" aria-labelledby="distribution-title" onMouseDown={(event) => event.stopPropagation()}>
      <button className="close-button" onClick={onClose} aria-label="Close"><X size={20} /></button>
      <p className="eyebrow">Underlying measurements</p>
      <h2 id="distribution-title">{displayFeature(row.feature)}</h2>
      <p className="modal-subtitle">Healthy controls compared with {displayDisease(row.disease)}</p>
      {error && <div className="alert">{error}</div>}
      {!data && !error && <div className="loading">Loading distributions…</div>}
      {data && <Histograms data={data} />}
      <div className="stat-strip">
        <span><small>Hedges’ g</small><b>{formatNumber(row.hedges_g)}</b></span>
        <span><small>95% CI</small><b>{formatNumber(row.hedges_g_ci_low)} to {formatNumber(row.hedges_g_ci_high)}</b></span>
        <span><small>FDR q</small><b>{formatQ(row.q_value)}</b></span>
        <span><small>AUROC</small><b>{formatNumber(row.roc_auc)}</b></span>
      </div>
    </section>
  </div>;
}

const pipeline = ["Facial images", "FaceKit", "125 geometric measurements", "219 disease cohorts vs healthy controls", "Statistical testing + FDR", "Effect size + stability + AUROC", "Quantitative phenotype atlas"];

export default function Dashboard() {
  const [summary, setSummary] = useState<Summary | null>(null);
  const [results, setResults] = useState<Result[]>([]);
  const [specificity, setSpecificity] = useState<Specificity[]>([]);
  const [diseases, setDiseases] = useState<string[]>([]);
  const [disease, setDisease] = useState("");
  const [diseaseQuery, setDiseaseQuery] = useState("");
  const [search, setSearch] = useState("");
  const [region, setRegion] = useState<Region>("All regions");
  const [priorityOnly, setPriorityOnly] = useState(false);
  const [showAll, setShowAll] = useState(false);
  const [selected, setSelected] = useState<Result | null>(null);
  const [error, setError] = useState("");
  const diseaseSelect = useRef<HTMLInputElement>(null);
  const featureInput = useRef<HTMLInputElement>(null);

  useEffect(() => {
    Promise.all([getJson<Summary>("/api/summary"), getJson<string[]>("/api/diseases"), getJson<Specificity[]>("/api/specificity?limit=200")])
      .then(([summaryData, diseaseData, specificityData]) => { setSummary(summaryData); setDiseases(diseaseData); setSpecificity(specificityData); })
      .catch(() => setError("The PRISM API is not reachable. Start the FastAPI server on port 8000."));
  }, []);

  useEffect(() => {
    if (!disease && !search) {
      return;
    }
    const params = new URLSearchParams({ analysis_type: "pairwise", limit: disease ? "200" : "500", sort_by: "pds_score", descending: "true" });
    if (disease) params.set("disease", disease);
    if (search) params.set("search", search);
    if (priorityOnly) params.set("high_priority", "true");
    getJson<{ items: Result[] }>(`/api/results?${params}`).then((data) => setResults(data.items)).catch(() => undefined);
  }, [disease, search, priorityOnly]);

  const filtered = useMemo(() => region === "All regions" ? results : results.filter((row) => featureRegion(row.feature) === region), [results, region]);
  const hasSelection = Boolean(disease || search.trim());
  const visibleSpecificity = useMemo(() => {
    const byFeature = new Map(specificity.map((row) => [row.feature, row]));
    if (disease) {
      // Keep the disease's PDS ordering, then display the global specificity of
      // those same features. This makes the card follow the selected cohort.
      return filtered.flatMap((row) => {
        const match = byFeature.get(row.feature);
        return match ? [match] : [];
      }).filter((row, index, rows) => rows.findIndex((item) => item.feature === row.feature) === index);
    }
    return specificity.filter((row) => {
      const matchesSearch = !search || row.feature.toLowerCase().includes(search.toLowerCase()) || displayFeature(row.feature).toLowerCase().includes(search.toLowerCase());
      const matchesRegion = region === "All regions" || featureRegion(row.feature) === region;
      return matchesSearch && matchesRegion;
    });
  }, [specificity, disease, filtered, search, region]);
  const visibleResults = useMemo(() => showAll ? filtered : filtered.slice(0, 12), [filtered, showAll]);
  const profile = useMemo(() => {
    if (!disease) return null;
    const significant = results.filter((row) => row.q_value != null && row.q_value < .05);
    return {
      n: Math.max(...results.map((row) => row.n_disease ?? 0)), significant: significant.length,
      priority: results.filter((row) => row.high_priority).length,
      increases: [...results].filter((row) => (row.hedges_g ?? 0) > 0).sort((a, b) => (b.hedges_g ?? 0) - (a.hedges_g ?? 0)).slice(0, 3),
      decreases: [...results].filter((row) => (row.hedges_g ?? 0) < 0).sort((a, b) => (a.hedges_g ?? 0) - (b.hedges_g ?? 0)).slice(0, 3),
    };
  }, [disease, results]);
  const regionScores = useMemo(() => REGIONS.slice(1).map((name) => {
    const rows = results.filter((row) => featureRegion(row.feature) === name);
    const mean = rows.reduce((sum, row) => sum + Math.min(Math.abs(row.hedges_g ?? 0), 1.5), 0) / Math.max(rows.length, 1);
    return { name, dots: Math.max(1, Math.min(5, Math.round(mean / 1.5 * 5))) };
  }), [results]);

  const focusControl = (ref: React.RefObject<HTMLInputElement | HTMLSelectElement | null>) => {
    document.querySelector("#explorer")?.scrollIntoView({ behavior: "smooth" });
    setTimeout(() => ref.current?.focus(), 450);
  };

  return <div className="app-shell">
    <aside className="sidebar">
      <div className="brand"><span className="brand-mark"><Dna size={22} weight="bold" /></span><span>PRISM</span></div>
      <nav aria-label="Primary navigation">
        <a className="nav-item active" href="#overview"><CirclesFour size={19} />Overview</a>
        <a className="nav-item" href="#explorer"><ChartScatter size={19} />Feature explorer</a>
        <a className="nav-item" href="#specificity"><Sparkle size={19} />Specificity</a>
        <a className="nav-item" href="/ancestry"><UsersThree size={19} />Ancestry analysis</a>
        <a className="nav-item" href="/distributions"><Images size={19} />Data distributions</a>
        <a className="nav-item" href="http://localhost:8000/docs"><Database size={19} />API documentation</a>
      </nav>
      <div className="sidebar-note"><Flask size={20} weight="duotone" /><div><strong>Research preview</strong><span>For exploratory analysis, not clinical diagnosis.</span></div></div>
    </aside>

    <main>
      <header className="topbar"><div><span className="status-dot" />Local analysis workspace</div><span className="version-tag">Quantitative phenotype atlas</span></header>
      <div className="content">
        <section className="hero" id="overview">
          <div><p className="kicker">Rare-disease phenotype intelligence</p><h1>Quantifying facial phenotypes across rare genetic disorders.</h1>
            <p>Explore statistically significant differences in {summary?.feature_count ?? 125} geometric facial measurements between healthy controls and {summary?.disease_count ?? 219} rare-disease cohorts.</p>
            <div className="hero-actions"><button onClick={() => focusControl(diseaseSelect)}>Explore a disease</button><button onClick={() => focusControl(featureInput)}>Explore a feature</button><a href="/distributions">View dataset distributions</a></div>
          </div><div className="hero-orbit" aria-hidden="true"><Dna size={68} weight="duotone" /></div>
        </section>

        {error && <div className="alert">{error}</div>}
        {summary?.is_demo && <div className="demo-banner"><Sparkle size={17} />Showing illustrative demo data. Run the analysis to populate this dashboard with your results.</div>}

        <section className="study-overview" aria-labelledby="study-heading"><div className="section-intro"><p className="eyebrow">Study overview</p><h2 id="study-heading">From facial images to an interpretable atlas</h2></div><div className="pipeline">{pipeline.map((step, index) => <div className="pipeline-step" key={step}><span>{index === 0 ? <Images size={18} /> : index + 1}</span><b>{step}</b>{index < pipeline.length - 1 && <ArrowRight className="pipeline-arrow" size={15} />}</div>)}</div></section>

        <section className="metric-grid" aria-label="Dataset summary">
          <MetricCard label="Features analyzed" value={(summary?.feature_count ?? 0).toLocaleString()} detail="Geometric phenotypes" icon={PulseIcon} />
          <MetricCard label="Disease cohorts" value={(summary?.disease_count ?? 0).toLocaleString()} detail="Compared with healthy" icon={Dna} tone="violet" />
          <MetricCard label="High priority" value={(summary?.high_priority_count ?? 0).toLocaleString()} detail="FDR + effect threshold" icon={CheckCircle} tone="amber" />
          <MetricCard label="Mean AUROC" value={formatNumber(summary?.average_auc ?? null)} detail="Univariate discrimination" icon={TrendUp} />
        </section>

        <section className="filter-panel" id="explorer">
          <div className="filter-title"><SlidersHorizontal size={19} /><span>Explore comparisons</span></div>
          <label className="search-field"><span>Feature measurement</span><MagnifyingGlass size={18} /><input ref={featureInput} value={search} onChange={(event) => setSearch(event.target.value)} placeholder="Search features…" /></label>
          <label className="search-field disease-search"><span>Disease cohort</span><MagnifyingGlass size={18} /><input ref={diseaseSelect} list="disease-options" value={diseaseQuery} placeholder="Search diseases…" autoComplete="off" onChange={(event) => {
            const value = event.target.value;
            setDiseaseQuery(value);
            const match = diseases.find((item) => displayDisease(item).toLowerCase() === value.trim().toLowerCase());
            setDisease(match ?? "");
            setShowAll(false);
          }} /><datalist id="disease-options">{diseases.map((item) => <option key={item} value={displayDisease(item)} />)}</datalist></label>
          <label><span>Facial region</span><select value={region} onChange={(event) => setRegion(event.target.value as Region)}>{REGIONS.map((item) => <option key={item}>{item}</option>)}</select></label>
          <label className="check-field"><input type="checkbox" checked={priorityOnly} onChange={(event) => setPriorityOnly(event.target.checked)} /><span>High priority only</span></label>
        </section>

        {!hasSelection && <section className="explorer-introduction">
          <div className="intro-copy"><p className="eyebrow">How to explore PRISM</p><h2>Start with a scientific question—not an undirected list.</h2><p>PRISM compares quantitative facial measurements from rare-disease cohorts with healthy controls. Select a disease to build its phenotype profile, or search for a measurement to see how it differs across disorders.</p><button onClick={() => focusControl(diseaseSelect)}>Choose a disease <ArrowRight size={16} /></button></div>
          <div className="intro-options">
            <article><span>01</span><div><h3>Select a disease</h3><p>See sample size, significant measurements, strongest increases and decreases, and anatomical-region patterns.</p></div></article>
            <article><span>02</span><div><h3>Inspect the evidence</h3><p>Compare Hedges’ g, confidence intervals, FDR q-values, AUROC, bootstrap stability, and PDS rankings.</p></div></article>
            <article><span>03</span><div><h3>View underlying values</h3><p>Click a feature or scatterplot point to compare its healthy and disease measurement distributions.</p></div></article>
          </div>
        </section>}

        {hasSelection && <>
        {profile && <section className="disease-profile">
          <div className="profile-summary"><p className="eyebrow">Disease profile</p><h2>{displayDisease(disease)}</h2><div className="profile-stats"><span><b>{profile.n}</b> images</span><span><b>{profile.significant} / {results.length}</b> measurements with q &lt; 0.05</span><span><b>{profile.priority}</b> high-priority phenotypes</span></div></div>
          <div className="signal-list"><h3><TrendUp size={16} /> Strongest increases</h3>{profile.increases.map((row) => <button key={row.id} onClick={() => setSelected(row)}><span>{displayFeature(row.feature)}</span><b>+{formatNumber(row.hedges_g)}</b></button>)}</div>
          <div className="signal-list"><h3><TrendDown size={16} /> Strongest decreases</h3>{profile.decreases.map((row) => <button key={row.id} onClick={() => setSelected(row)}><span>{displayFeature(row.feature)}</span><b>{formatNumber(row.hedges_g)}</b></button>)}</div>
        </section>}

        {disease && <section className="region-panel"><div><p className="eyebrow">Anatomical summary</p><h2>Signal by facial region</h2><p>Strength summarizes absolute standardized effects within each region. Select a region to filter the atlas.</p></div><div className="region-grid">{regionScores.map((item) => <button key={item.name} className={region === item.name ? "active" : ""} onClick={() => setRegion(item.name)}><span>{item.name}</span><i>{Array.from({ length: 5 }, (_, index) => <b className={index < item.dots ? "filled" : ""} key={index} />)}</i><small>{item.dots >= 4 ? "Strong" : item.dots >= 2 ? "Moderate" : "Limited"}</small></button>)}</div></section>}

        <section className="chart-grid">
          <article className="panel"><div className="panel-heading"><div><p className="eyebrow">Discriminative signal</p><h2>Effect size vs. AUROC</h2><p className="panel-description">Hover for statistical evidence; click a point to inspect its underlying distributions.</p></div><span className="legend"><i />High priority</span></div><EffectScatter results={filtered} onSelect={setSelected} /></article>
          <article className="panel" id="specificity"><div className="panel-heading align-start"><div><p className="eyebrow">Distinctive markers</p><h2>Feature specificity <span className="info-tip" tabIndex={0}><Info size={15} /><span>Specificity = 1 − (diseases with a significant moderate-or-large effect ÷ diseases tested). Higher values indicate a more disease-selective feature.</span></span></h2><p className="panel-description">Measures how selectively a facial measurement is associated with particular disorders rather than broadly altered across many disorders.</p></div><span className="panel-tag">{disease ? "Selected cohort" : region !== "All regions" ? region : "Top weighted"}</span></div>{visibleSpecificity.length ? <SpecificityBars data={visibleSpecificity} /> : <div className="empty-chart">No specificity results match the current filters.</div>}</article>
        </section>

        <section className="panel results-panel" id="evidence">
          <div className="panel-heading"><div><p className="eyebrow">Evidence table</p><h2>Feature comparisons</h2><p className="panel-description">Ranked by Phenotype Difference Score (PDS), highest first.</p></div><span className="result-count">{filtered.length} results</span></div>
          <div className="table-wrap"><table><thead><tr><th>Feature</th>{!disease && <th>Disease</th>}<th>Direction</th><th>Hedges’ g</th><th>95% CI</th><th>FDR q</th><th>AUROC</th><th>PDS <span className="info-tip" tabIndex={0}><Info size={12} /><span><b>Composite ranking index.</b> PDS integrates standardized effect magnitude, statistical confidence, bootstrap direction stability, and sample reliability. It prioritizes features and does not replace significance testing.</span></span></th><th>Stability</th><th>Status</th></tr></thead>
            <tbody>{visibleResults.map((row) => <tr key={row.id}><td><button className="feature-link" onClick={() => setSelected(row)}>{displayFeature(row.feature)}</button></td>{!disease && <td>{displayDisease(row.disease)}</td>}<td><span className={(row.hedges_g ?? 0) < 0 ? "direction down" : "direction up"}>{(row.hedges_g ?? 0) < 0 ? "↓ Smaller" : "↑ Larger"}</span></td><td className={(row.hedges_g ?? 0) < 0 ? "negative" : "positive"}>{formatNumber(row.hedges_g)}</td><td>{formatNumber(row.hedges_g_ci_low)} to {formatNumber(row.hedges_g_ci_high)}</td><td>{formatQ(row.q_value)}</td><td><span className="auc-pill">{formatNumber(row.roc_auc)}</span></td><td><span className="pds-pill">{row.pds_score == null ? "—" : Math.round(row.pds_score)}</span></td><td>{row.rank_stability == null ? "—" : `${Math.round(row.rank_stability * 100)}%`}</td><td>{row.high_priority ? <span className="priority"><CheckCircle size={14} weight="fill" />Priority</span> : <span className="muted-status">Monitor</span>}</td></tr>)}</tbody></table></div>
          {filtered.length > 12 && <button className="text-button" onClick={() => setShowAll((value) => !value)}>{showAll ? "Show top results" : "View all comparisons"} <ArrowRight size={16} /></button>}
        </section>
        </>}
      </div>
    </main>
    {selected && <DistributionModal row={selected} onClose={() => setSelected(null)} />}
  </div>;
}
