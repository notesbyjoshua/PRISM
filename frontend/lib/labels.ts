const WORDS: Record<string, string> = {
  eb: "eyebrow",
  inter: "inter",
  canthal: "canthal",
  pupillary: "pupillary",
  r: "right",
  l: "left",
  asym: "asymmetry",
  ala: "alar",
  c2: "curvature",
  x: "horizontal",
  y: "vertical",
};

export const FEATURE_DISPLAY_NAMES: Record<string, string> = {
  eb_inter_distance: "Inter-eyebrow distance",
  inter_pupillary_distance: "Interpupillary distance",
  inter_canthal_distance: "Intercanthal distance",
  outer_canthal_distance: "Outer canthal distance",
  forehead_width_mid: "Mid-forehead width",
  forehead_width_upper: "Upper-forehead width",
  jaw_to_forehead_ratio: "Jaw-to-forehead ratio",
  nose_to_face_area_ratio: "Nose-to-face area ratio",
  canthal_to_pupillary_ratio: "Canthal-to-pupillary ratio",
};

export function displayDisease(value: string | null) {
  if (!value) return "—";
  const clean = value.replace(/__[a-z0-9]+$/i, "").replaceAll("_", " ").trim();
  return clean.toLowerCase().replace(/(^|[\s/-])\p{L}/gu, (letter) => letter.toUpperCase());
}

export function displayFeature(value: string) {
  if (FEATURE_DISPLAY_NAMES[value]) return FEATURE_DISPLAY_NAMES[value];
  const words = value.split("_").map((word) => WORDS[word] ?? word);
  const label = words.join(" ").replace("mean", "mean");
  return label.charAt(0).toUpperCase() + label.slice(1);
}

export const REGIONS = [
  "All regions",
  "Eyes & orbital region",
  "Eyebrows",
  "Nose",
  "Mouth & lips",
  "Midface",
  "Jaw & chin",
  "Forehead",
  "Global facial proportions",
] as const;

export type Region = (typeof REGIONS)[number];

export function featureRegion(feature: string): Region {
  if (feature.startsWith("eb_")) return "Eyebrows";
  if (/^(eye_|inter_canthal|inter_pupillary|outer_canthal|canthal_|iris_|gaze_)/.test(feature)) return "Eyes & orbital region";
  if (/^(nose_|nostril_|nasolabial|ala_|columella)/.test(feature)) return "Nose";
  if (/^(philtrum_|mouth_|upper_vermilion|lower_vermilion|vermilion_|cupid_|upper_lip|lower_lip|lip_)/.test(feature)) return "Mouth & lips";
  if (/^(midface_|malar_|cheek_)/.test(feature)) return "Midface";
  if (/^(jaw_|chin_)/.test(feature)) return "Jaw & chin";
  if (/^(forehead_|hairline_)/.test(feature)) return "Forehead";
  return "Global facial proportions";
}
