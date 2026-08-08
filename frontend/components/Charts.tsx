"use client";

import {
  Bar, BarChart, CartesianGrid, Cell, ReferenceLine, ResponsiveContainer,
  Scatter, ScatterChart, Tooltip, XAxis, YAxis,
} from "recharts";
import { Result, Specificity } from "@/lib/api";
import { displayDisease, displayFeature } from "@/lib/labels";

const tooltipStyle = {
  background: "#17201f", border: "1px solid #344340", borderRadius: 12,
  color: "#fff", fontSize: 12,
};

const number = (value: number | null, digits = 2) => value == null ? "—" : value.toFixed(digits);
const qValue = (value: number | null) => value == null ? "—" : value < .001 ? "<0.001" : value.toPrecision(2);

function EffectTooltip({ active, payload }: { active?: boolean; payload?: Array<{ payload: Result & { effect: number } }> }) {
  if (!active || !payload?.length) return null;
  const row = payload[0].payload;
  return <div className="science-tooltip">
    <strong>{displayFeature(row.feature)}</strong>
    <span>{displayDisease(row.disease)}</span>
    <dl>
      <dt>Hedges’ g</dt><dd>{(row.hedges_g ?? 0) > 0 ? "+" : ""}{number(row.hedges_g)}</dd>
      <dt>95% CI</dt><dd>{number(row.hedges_g_ci_low)}–{number(row.hedges_g_ci_high)}</dd>
      <dt>AUROC</dt><dd>{number(row.roc_auc)}</dd>
      <dt>FDR q</dt><dd>{qValue(row.q_value)}</dd>
      <dt>PDS</dt><dd>{row.pds_score == null ? "—" : Math.round(row.pds_score)}</dd>
      <dt>Stability</dt><dd>{row.rank_stability == null ? "—" : `${Math.round(row.rank_stability * 100)}%`}</dd>
    </dl>
    <b>{(row.hedges_g ?? 0) >= 0 ? "↑ Larger" : "↓ Smaller"} in {displayDisease(row.disease)}</b>
    <small>Click to inspect the distributions</small>
  </div>;
}

export function EffectScatter({ results, onSelect }: { results: Result[]; onSelect: (row: Result) => void }) {
  const data = results.filter((row) => row.roc_auc !== null && row.hedges_g !== null)
    .map((row) => ({ ...row, effect: Math.abs(row.hedges_g ?? 0) }));
  return <ResponsiveContainer width="100%" height={280}>
    <ScatterChart margin={{ top: 12, right: 12, bottom: 10, left: -12 }}>
      <CartesianGrid stroke="#e6ebe8" strokeDasharray="3 5" />
      <XAxis type="number" dataKey="effect" name="|Hedges' g|" domain={[0, "auto"]} tick={{ fill: "#71807c", fontSize: 11 }} label={{ value: "Absolute Hedges’ g", position: "insideBottom", offset: -5, fill: "#71807c", fontSize: 10 }} />
      <YAxis type="number" dataKey="roc_auc" name="AUROC" domain={[0.5, 1]} tick={{ fill: "#71807c", fontSize: 11 }} />
      <ReferenceLine x={0.8} stroke="#c87d20" strokeDasharray="5 4" label={{ value: "large effect", fill: "#9a641e", fontSize: 10, position: "insideTopRight" }} />
      <Tooltip cursor={{ strokeDasharray: "3 3" }} content={<EffectTooltip />} />
      <Scatter data={data} onClick={(point) => onSelect((point as unknown as { payload: Result }).payload)} cursor="pointer">
        {data.map((row) => <Cell key={row.id} fill={row.high_priority ? "#0d8f78" : "#a8cfc6"} />)}
      </Scatter>
    </ScatterChart>
  </ResponsiveContainer>;
}

export function SpecificityBars({ data }: { data: Specificity[] }) {
  const compact = data.slice(0, 6).map((row) => ({ ...row, short: displayFeature(row.feature) }));
  return <ResponsiveContainer width="100%" height={280}>
    <BarChart data={compact} layout="vertical" margin={{ top: 8, right: 12, bottom: 2, left: 12 }}>
      <CartesianGrid stroke="#e6ebe8" horizontal={false} strokeDasharray="3 5" />
      <XAxis type="number" domain={[0, 1]} tick={{ fill: "#71807c", fontSize: 11 }} />
      <YAxis type="category" dataKey="short" width={132} tick={{ fill: "#53635f", fontSize: 10 }} />
      <Tooltip contentStyle={tooltipStyle} formatter={(value) => [Number(value).toFixed(2), "Specificity"]} />
      <Bar dataKey="specificity_score" name="Specificity" fill="#6d5bd0" radius={[0, 6, 6, 0]} />
    </BarChart>
  </ResponsiveContainer>;
}
