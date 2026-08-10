"use client";

import { useEffect, useMemo, useState } from "react";
import Link from "next/link";
import { ArrowLeft, ChartBar, CheckCircle, Database, Dna, Flask, Info, MagnifyingGlass, UsersThree } from "@phosphor-icons/react";
import { Bar, BarChart, CartesianGrid, ResponsiveContainer, Tooltip, XAxis, YAxis } from "recharts";
import { displayDisease, displayFeature } from "@/lib/labels";
import { getJson } from "@/lib/api";

type Layer = "broad" | "subcategory";
type Overview = {
  group_counts: { group: string; healthy_images: number }[];
  healthy_significant_features: number;
  healthy_pairwise_findings: number;
  matched_high_priority_findings: number;
  syndrome_group_comparisons: number;
  groups: string[];
  diseases: string[];
  interaction_counts: Record<string, number>;
  harmonization: { dataset: string; layer: string; total_rows: number; kept_rows: number; excluded_rows: number; percent_kept: number }[];
};
type MatchedResult = { syndrome_name: string; group: string; feature: string; n_disease: number; n_healthy_matched: number; evidence_level: string; disease_mean: number; healthy_mean: number; hedges_g: number; hedges_g_ci_lower: number; hedges_g_ci_upper: number; global_fdr_q_value: number; high_priority: boolean };
type HealthyEffect = { group_1: string; group_2: string; feature: string; n_group_1: number; n_group_2: number; mean_group_1: number; mean_group_2: number; hedges_g: number; effect_category: string; direction: string; within_pair_fdr_q_value: number };
type Interaction = { syndrome_name: string; feature: string; group_1: string; group_2: string; n_disease_group_1: number; n_disease_group_2: number; evidence_level: string; hedges_g_group_1: number; hedges_g_group_2: number; delta_hedges_g_group2_minus_group1: number; interaction_beta_group2_minus_group1: number; interaction_ci_lower: number; interaction_ci_upper: number; global_fdr_q_value: number; direction_reversal: boolean; high_priority_interaction: boolean; high_confidence_interaction: boolean };
type InteractionSummary = { syndrome_name: string; comparisons_tested: number; significant_interactions: number; global_significant_interactions: number; high_priority_interactions: number; high_confidence_interactions: number; direction_reversals: number; max_abs_delta_g: number; median_abs_delta_g: number };
type InteractionSet = "all" | "within_syndrome_fdr" | "global_fdr" | "high_priority" | "high_confidence" | "direction_reversals";

const number = (value: number | null, digits = 2) => value == null ? "—" : Number(value).toFixed(digits);
const qValue = (value: number | null) => value == null ? "—" : value < .001 ? "<0.001" : Number(value).toPrecision(2);

export default function AncestryDashboard() {
  const [layer, setLayer] = useState<Layer>("broad");
  const [overview, setOverview] = useState<Overview | null>(null);
  const [matched, setMatched] = useState<MatchedResult[]>([]);
  const [healthy, setHealthy] = useState<HealthyEffect[]>([]);
  const [interactions, setInteractions] = useState<Interaction[]>([]);
  const [interactionSummaries, setInteractionSummaries] = useState<InteractionSummary[]>([]);
  const [interactionSet, setInteractionSet] = useState<InteractionSet>("high_priority");
  const [disease, setDisease] = useState("");
  const [diseaseQuery, setDiseaseQuery] = useState("");
  const [group, setGroup] = useState("");
  const [error, setError] = useState("");

  useEffect(() => {
    Promise.all([
      getJson<Overview>(`/api/ancestry/overview?layer=${layer}`),
      getJson<HealthyEffect[]>(`/api/ancestry/healthy-effects?layer=${layer}&limit=20`),
      getJson<InteractionSummary[]>(`/api/ancestry/interaction-summaries?layer=${layer}&limit=12`),
    ]).then(([summary, effects, interactionSummary]) => { setOverview(summary); setHealthy(effects); setInteractionSummaries(interactionSummary); }).catch(() => setError("The ancestry-analysis results could not be loaded."));
  }, [layer]);

  useEffect(() => {
    const params = new URLSearchParams({ layer, limit: "30" });
    if (disease) params.set("disease", disease);
    if (group) params.set("group", group);
    getJson<MatchedResult[]>(`/api/ancestry/matched-results?${params}`).then(setMatched).catch(() => setError("The matched results could not be loaded."));
  }, [layer, disease, group]);

  useEffect(() => {
    const params = new URLSearchParams({ layer, result_set: interactionSet, limit: "30" });
    if (disease) params.set("disease", disease);
    if (group) params.set("group", group);
    getJson<Interaction[]>(`/api/ancestry/interactions?${params}`).then(setInteractions).catch(() => setError("Analysis C interactions could not be loaded."));
  }, [layer, interactionSet, disease, group]);

  const harmonization = useMemo(() => overview?.harmonization.filter((row) => row.layer.startsWith(layer === "broad" ? "layer1" : "layer2")) ?? [], [overview, layer]);
  const changeLayer = (nextLayer: Layer) => {
    setLayer(nextLayer);
    setDisease("");
    setDiseaseQuery("");
    setGroup("");
    setInteractionSet("high_priority");
    setError("");
  };

  return <div className="app-shell">
    <aside className="sidebar">
      <div className="brand"><span className="brand-mark"><Dna size={22} weight="bold" /></span><span>PRISM</span></div>
      <nav aria-label="Primary navigation"><Link className="nav-item" href="/"><ArrowLeft size={19} />Phenotype atlas</Link><Link className="nav-item active" href="/ancestry"><UsersThree size={19} />Ancestry analysis</Link><Link className="nav-item" href="/distributions"><ChartBar size={19} />Data distributions</Link><a className="nav-item" href="http://localhost:8000/docs"><Database size={19} />API documentation</a></nav>
      <div className="sidebar-note"><Flask size={20} weight="duotone" /><div><strong>Matched analysis</strong><span>Separates ancestry-associated variation from disease-associated effects.</span></div></div>
    </aside>
    <main><header className="topbar"><div><span className="status-dot" />Local analysis workspace</div><span className="version-tag">Ancestry-aware evidence</span></header>
      <div className="content ancestry-content">
        <section className="ancestry-hero"><p className="kicker">Ancestry analysis</p><h1>Separating ancestry-associated variation from disease phenotypes.</h1><p>Explore healthy-control differences across ancestry groups and disease comparisons matched to healthy controls from the same group.</p><div className="layer-switch" role="group" aria-label="Ancestry grouping level"><button className={layer === "broad" ? "active" : ""} onClick={() => changeLayer("broad")}>Broad ancestry groups</button><button className={layer === "subcategory" ? "active" : ""} onClick={() => changeLayer("subcategory")}>Detailed subcategories</button></div></section>
        {error && <div className="alert">{error}</div>}
        <section className="ancestry-metrics"><article><small>Healthy groups</small><b>{overview?.group_counts.length ?? 0}</b></article><article><small>Healthy features with FDR significance</small><b>{overview?.healthy_significant_features ?? 0}</b></article><article><small>Matched high-priority findings</small><b>{(overview?.matched_high_priority_findings ?? 0).toLocaleString()}</b></article><article><small>Syndrome–group comparisons</small><b>{(overview?.syndrome_group_comparisons ?? 0).toLocaleString()}</b></article></section>

        <section className="ancestry-grid"><article className="panel"><div className="panel-heading"><div><p className="eyebrow">Reference representation</p><h2>Healthy images by ancestry group</h2><p className="panel-description">The healthy-control sample available for ancestry-matched comparisons.</p></div></div><ResponsiveContainer width="100%" height={285}><BarChart data={overview?.group_counts ?? []} margin={{ top: 10, right: 10, bottom: 35, left: -10 }}><CartesianGrid stroke="#e6ebe8" vertical={false} strokeDasharray="3 5" /><XAxis dataKey="group" angle={-25} textAnchor="end" interval={0} tick={{ fontSize: 9, fill: "#71807c" }} /><YAxis tick={{ fontSize: 9, fill: "#71807c" }} /><Tooltip /><Bar dataKey="healthy_images" name="Healthy images" fill="#0d8f78" radius={[5,5,0,0]} /></BarChart></ResponsiveContainer></article>
          <article className="panel"><div className="panel-heading"><div><p className="eyebrow">Data harmonization</p><h2>Rows retained after mapping</h2><p className="panel-description">Coverage after source ethnicity labels were harmonized to the selected ancestry layer.</p></div></div><div className="harmonization-list">{harmonization.map((row) => <div key={row.dataset}><span><b>{row.dataset}</b><small>{row.kept_rows.toLocaleString()} of {row.total_rows.toLocaleString()} rows</small></span><strong>{number(row.percent_kept, 1)}%</strong><i><b style={{ width: `${row.percent_kept}%` }} /></i></div>)}</div></article></section>

        <section className="panel ancestry-results"><div className="panel-heading"><div><p className="eyebrow">Ancestry-matched disease evidence</p><h2>Top disease vs. matched-healthy effects</h2><p className="panel-description">Disease cohorts are compared only with healthy controls assigned to the same ancestry group.</p></div><span className="result-count">{matched.length} results</span></div>
          <div className="ancestry-filters"><label className="search-field"><span>Disease</span><MagnifyingGlass size={18} /><input list="ancestry-diseases" value={diseaseQuery} placeholder="Search diseases…" onChange={(event) => { const value = event.target.value; setDiseaseQuery(value); setDisease(overview?.diseases.find((item) => displayDisease(item).toLowerCase() === value.trim().toLowerCase()) ?? ""); }} /><datalist id="ancestry-diseases">{overview?.diseases.map((item) => <option key={item} value={displayDisease(item)} />)}</datalist></label><label><span>Ancestry group</span><select value={group} onChange={(event) => setGroup(event.target.value)}><option value="">All groups</option>{overview?.groups.map((item) => <option key={item}>{item}</option>)}</select></label></div>
          <div className="table-wrap"><table><thead><tr><th>Disease</th><th>Ancestry group</th><th>Feature</th><th>n disease / healthy</th><th>Means</th><th>Hedges’ g</th><th>95% CI</th><th>Global FDR q</th><th>Evidence</th></tr></thead><tbody>{matched.map((row, index) => <tr key={`${row.syndrome_name}-${row.group}-${row.feature}-${index}`}><td>{displayDisease(row.syndrome_name)}</td><td>{row.group}</td><td><strong>{displayFeature(row.feature)}</strong></td><td>{row.n_disease} / {row.n_healthy_matched}</td><td>{number(row.disease_mean, 3)} / {number(row.healthy_mean, 3)}</td><td className={row.hedges_g < 0 ? "negative" : "positive"}>{number(row.hedges_g)}</td><td>{number(row.hedges_g_ci_lower)} to {number(row.hedges_g_ci_upper)}</td><td>{qValue(row.global_fdr_q_value)}</td><td>{row.high_priority ? <span className="priority"><CheckCircle size={14} weight="fill" />High priority</span> : row.evidence_level}</td></tr>)}</tbody></table></div>
        </section>

        <section className="analysis-c-section">
          <div className="analysis-c-heading"><div><p className="eyebrow">Analysis C · disease × ancestry interaction</p><h2>Does a disease-associated facial effect differ across ancestry groups? <span className="info-tip" tabIndex={0}><Info size={16} /><span>An interaction tests whether the disease-versus-healthy difference changes between two ancestry groups. It is not simply a difference between ancestry groups.</span></span></h2><p>Analysis C compares each feature’s standardized disease effect between ancestry groups. Large absolute Δg values indicate stronger effect modification; direction reversals indicate that the disease effect changes sign.</p></div></div>
          <div className="interaction-metrics">
            <article><small>All tested</small><b>{(overview?.interaction_counts.all ?? 0).toLocaleString()}</b></article>
            <article><small>Within-syndrome FDR</small><b>{(overview?.interaction_counts.within_syndrome_fdr ?? 0).toLocaleString()}</b></article>
            <article><small>Global FDR</small><b>{(overview?.interaction_counts.global_fdr ?? 0).toLocaleString()}</b></article>
            <article><small>High priority</small><b>{(overview?.interaction_counts.high_priority ?? 0).toLocaleString()}</b></article>
            <article><small>High confidence</small><b>{(overview?.interaction_counts.high_confidence ?? 0).toLocaleString()}</b></article>
            <article><small>Direction reversals</small><b>{(overview?.interaction_counts.direction_reversals ?? 0).toLocaleString()}</b></article>
          </div>
          <section className="panel ancestry-results"><div className="panel-heading"><div><p className="eyebrow">Interaction evidence browser</p><h2>Ancestry-dependent disease effects</h2><p className="panel-description">The disease and ancestry filters above also apply here.</p></div><span className="result-count">{interactions.length} results</span></div>
            <div className="interaction-filter"><label><span>Analysis C result set</span><select value={interactionSet} onChange={(event) => setInteractionSet(event.target.value as InteractionSet)}><option value="high_priority">High-priority interactions</option><option value="high_confidence">High-confidence interactions</option><option value="global_fdr">Global FDR significant</option><option value="within_syndrome_fdr">Within-syndrome FDR significant</option><option value="direction_reversals">Direction reversals</option><option value="all">All tested interactions</option></select></label><p>Showing the largest absolute difference in Hedges’ g first.</p></div>
            <div className="table-wrap"><table><thead><tr><th>Disease</th><th>Feature</th><th>Groups</th><th>Disease n</th><th>g by group</th><th>Δ Hedges’ g</th><th>Interaction β</th><th>β 95% CI</th><th>Global FDR q</th><th>Flags</th></tr></thead><tbody>{interactions.map((row, index) => <tr key={`${row.syndrome_name}-${row.feature}-${row.group_1}-${row.group_2}-${index}`}><td>{displayDisease(row.syndrome_name)}</td><td><strong>{displayFeature(row.feature)}</strong></td><td>{row.group_1} → {row.group_2}</td><td>{row.n_disease_group_1} / {row.n_disease_group_2}</td><td>{number(row.hedges_g_group_1)} / {number(row.hedges_g_group_2)}</td><td className={row.delta_hedges_g_group2_minus_group1 < 0 ? "negative" : "positive"}>{number(row.delta_hedges_g_group2_minus_group1)}</td><td>{number(row.interaction_beta_group2_minus_group1, 3)}</td><td>{number(row.interaction_ci_lower, 3)} to {number(row.interaction_ci_upper, 3)}</td><td>{qValue(row.global_fdr_q_value)}</td><td><span className="interaction-flags">{row.direction_reversal && <b className="reversal">Reversal</b>}{row.high_priority_interaction && <b>Priority</b>}{row.high_confidence_interaction && <b>High confidence</b>}</span></td></tr>)}</tbody></table></div>
          </section>
          <section className="panel ancestry-results"><div className="panel-heading"><div><p className="eyebrow">Syndrome summary</p><h2>Diseases with the most high-priority interactions</h2><p className="panel-description">Aggregated directly from the Analysis C syndrome summary.</p></div></div><div className="table-wrap"><table><thead><tr><th>Disease</th><th>Tested</th><th>Within-syndrome significant</th><th>Global significant</th><th>High priority</th><th>High confidence</th><th>Direction reversals</th><th>Max |Δg|</th><th>Median |Δg|</th></tr></thead><tbody>{interactionSummaries.map((row) => <tr key={row.syndrome_name}><td><strong>{displayDisease(row.syndrome_name)}</strong></td><td>{row.comparisons_tested.toLocaleString()}</td><td>{row.significant_interactions}</td><td>{row.global_significant_interactions}</td><td>{row.high_priority_interactions}</td><td>{row.high_confidence_interactions}</td><td>{row.direction_reversals}</td><td>{number(row.max_abs_delta_g)}</td><td>{number(row.median_abs_delta_g)}</td></tr>)}</tbody></table></div></section>
        </section>

        <section className="panel ancestry-results"><div className="panel-heading"><div><p className="eyebrow">Healthy-control baseline</p><h2>Largest ancestry-associated feature differences</h2><p className="panel-description">These differences occur among healthy controls and may represent ancestry-associated facial variation rather than disease effects.</p></div><span className="result-count">Top {healthy.length}</span></div><div className="table-wrap"><table><thead><tr><th>Groups</th><th>Feature</th><th>n</th><th>Group means</th><th>Hedges’ g</th><th>Effect</th><th>FDR q</th></tr></thead><tbody>{healthy.map((row, index) => <tr key={`${row.group_1}-${row.group_2}-${row.feature}-${index}`}><td>{row.group_1} vs. {row.group_2}</td><td><strong>{displayFeature(row.feature)}</strong></td><td>{row.n_group_1} / {row.n_group_2}</td><td>{number(row.mean_group_1, 3)} / {number(row.mean_group_2, 3)}</td><td className={row.hedges_g < 0 ? "negative" : "positive"}>{number(row.hedges_g)}</td><td>{row.effect_category}</td><td>{qValue(row.within_pair_fdr_q_value)}</td></tr>)}</tbody></table></div></section>
      </div>
    </main>
  </div>;
}
