"use client";

import {
  Bar,
  BarChart,
  CartesianGrid,
  Cell,
  ResponsiveContainer,
  Scatter,
  ScatterChart,
  Tooltip,
  XAxis,
  YAxis,
} from "recharts";
import { Result, Specificity } from "@/lib/api";

const tooltipStyle = {
  background: "#17201f",
  border: "1px solid #344340",
  borderRadius: 12,
  color: "#fff",
  fontSize: 12,
};

export function EffectScatter({ results }: { results: Result[] }) {
  const data = results
    .filter((row) => row.roc_auc !== null && row.hedges_g !== null)
    .map((row) => ({ ...row, effect: Math.abs(row.hedges_g ?? 0), name: row.feature }));

  return (
    <ResponsiveContainer width="100%" height={260}>
      <ScatterChart margin={{ top: 12, right: 12, bottom: 4, left: -16 }}>
        <CartesianGrid stroke="#e6ebe8" strokeDasharray="3 5" />
        <XAxis type="number" dataKey="effect" name="|Hedges' g|" domain={[0, "auto"]} tick={{ fill: "#71807c", fontSize: 11 }} />
        <YAxis type="number" dataKey="roc_auc" name="AUROC" domain={[0.5, 1]} tick={{ fill: "#71807c", fontSize: 11 }} />
        <Tooltip cursor={{ strokeDasharray: "3 3" }} contentStyle={tooltipStyle} />
        <Scatter data={data} fill="#169b82">
          {data.map((row) => <Cell key={row.id} fill={row.high_priority ? "#0d8f78" : "#a8cfc6"} />)}
        </Scatter>
      </ScatterChart>
    </ResponsiveContainer>
  );
}

export function SpecificityBars({ data }: { data: Specificity[] }) {
  const compact = data.slice(0, 6).map((row) => ({
    ...row,
    short: row.feature.replaceAll("_", " "),
  }));
  return (
    <ResponsiveContainer width="100%" height={260}>
      <BarChart data={compact} layout="vertical" margin={{ top: 8, right: 12, bottom: 2, left: 12 }}>
        <CartesianGrid stroke="#e6ebe8" horizontal={false} strokeDasharray="3 5" />
        <XAxis type="number" domain={[0, 1]} tick={{ fill: "#71807c", fontSize: 11 }} />
        <YAxis type="category" dataKey="short" width={112} tick={{ fill: "#53635f", fontSize: 10 }} />
        <Tooltip contentStyle={tooltipStyle} />
        <Bar dataKey="specificity_score" name="Specificity" fill="#6d5bd0" radius={[0, 6, 6, 0]} />
      </BarChart>
    </ResponsiveContainer>
  );
}

