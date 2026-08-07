import { Icon } from "@phosphor-icons/react";

type Props = {
  label: string;
  value: string;
  detail: string;
  icon: Icon;
  tone?: "mint" | "violet" | "amber";
};

export default function MetricCard({ label, value, detail, icon: MetricIcon, tone = "mint" }: Props) {
  return (
    <article className="metric-card">
      <div className={`metric-icon ${tone}`}><MetricIcon size={20} weight="duotone" /></div>
      <div>
        <p className="eyebrow">{label}</p>
        <strong className="metric-value">{value}</strong>
        <p className="metric-detail">{detail}</p>
      </div>
    </article>
  );
}

