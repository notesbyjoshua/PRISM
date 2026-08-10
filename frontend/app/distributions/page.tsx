"use client";

import { useEffect, useMemo, useState } from "react";
import Link from "next/link";
import { ArrowLeft, ChartBar, Dna, Flask, Images, UsersThree, WarningCircle } from "@phosphor-icons/react";
import { API_URL, getJson } from "@/lib/api";

type DistributionPlot = {
  filename: string;
  title: string;
  description: string;
  category: "Dataset overview" | "Ethnicity comparisons";
};

export default function DistributionsPage() {
  const [plots, setPlots] = useState<DistributionPlot[]>([]);
  const [category, setCategory] = useState("All charts");
  const [error, setError] = useState("");

  useEffect(() => {
    getJson<DistributionPlot[]>("/api/distribution-plots")
      .then(setPlots)
      .catch(() => setError("The distribution charts could not be loaded. Make sure the API server is running."));
  }, []);

  const visible = useMemo(() => category === "All charts" ? plots : plots.filter((plot) => plot.category === category), [plots, category]);

  return <div className="app-shell">
    <aside className="sidebar">
      <div className="brand"><span className="brand-mark"><Dna size={22} weight="bold" /></span><span>PRISM</span></div>
      <nav aria-label="Primary navigation">
        <Link className="nav-item" href="/"><ArrowLeft size={19} />Phenotype atlas</Link>
        <Link className="nav-item" href="/ancestry"><UsersThree size={19} />Ancestry analysis</Link>
        <Link className="nav-item active" href="/distributions"><ChartBar size={19} />Data distributions</Link>
        <a className="nav-item" href="http://localhost:8000/docs"><Images size={19} />API documentation</a>
      </nav>
      <div className="sidebar-note"><Flask size={20} weight="duotone" /><div><strong>Dataset overview</strong><span>Charts are generated from the configured input CSV.</span></div></div>
    </aside>

    <main>
      <header className="topbar"><div><span className="status-dot" />Local analysis workspace</div><span className="version-tag">Dataset distributions</span></header>
      <div className="content distributions-content">
        <section className="distribution-hero">
          <div><p className="kicker">Cohort composition</p><h1>Explore the dataset distributions.</h1><p>Review disease cohort sizes, ethnicity representation, and pairwise ethnicity counts across disorders.</p></div>
          <div className="distribution-total"><b>{plots.length}</b><span>generated charts</span></div>
        </section>

        <div className="chart-tabs" role="group" aria-label="Filter distribution charts">
          {["All charts", "Dataset overview", "Ethnicity comparisons"].map((item) => <button key={item} className={category === item ? "active" : ""} onClick={() => setCategory(item)}>{item}</button>)}
        </div>

        {error && <div className="alert"><WarningCircle size={17} />{error}</div>}
        {!error && plots.length === 0 && <div className="empty-distributions"><ChartBar size={28} /><h2>No distribution charts found</h2><p>Run <code>data_distribution.py</code>, then refresh this page.</p></div>}

        <section className="distribution-cards" aria-live="polite">
          {visible.map((plot) => <article className={plot.category === "Dataset overview" ? "distribution-card overview-chart" : "distribution-card"} key={plot.filename}>
            <div className="distribution-card-copy"><span>{plot.category}</span><h2>{plot.title}</h2><p>{plot.description}</p></div>
            {/* Generated analysis output is served directly by the local FastAPI server. */}
            {/* eslint-disable-next-line @next/next/no-img-element */}
            <img src={`${API_URL}/api/distribution-plots/${encodeURIComponent(plot.filename)}`} alt={`${plot.title} data visualization`} loading="lazy" />
          </article>)}
        </section>
      </div>
    </main>
  </div>;
}
