export const API_URL = process.env.NEXT_PUBLIC_API_URL ?? "http://localhost:8000";

export type Summary = {
  result_count: number;
  feature_count: number;
  disease_count: number;
  high_priority_count: number;
  average_auc: number | null;
  is_demo: boolean;
  data_source: string;
};

export type Result = {
  id: number;
  comparison: string;
  disease: string | null;
  feature: string;
  n_disease: number | null;
  n_healthy: number | null;
  hedges_g: number | null;
  hedges_g_ci_low: number | null;
  hedges_g_ci_high: number | null;
  cliffs_delta: number | null;
  roc_auc: number | null;
  q_value: number | null;
  pds_score: number | null;
  high_priority: boolean;
  rank_stability: number | null;
  robust_hedges_g: number | null;
};

export type Specificity = {
  feature: string;
  diseases_tested: number;
  diseases_with_effect: number;
  specificity_score: number | null;
  weighted_effect_score: number | null;
  median_abs_hedges_g: number | null;
};

export async function getJson<T>(path: string): Promise<T> {
  const response = await fetch(`${API_URL}${path}`);
  if (!response.ok) throw new Error(`API request failed: ${response.status}`);
  return response.json();
}
