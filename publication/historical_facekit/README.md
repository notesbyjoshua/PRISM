# Historical FaceKit publication source

`extractor.py` is preserved byte-for-byte from author-confirmed FaceKit revision `e474628f3cb17a83edad59ac6e1d9e5438753032`, path `src/facekit/core/geometric/extractor.py`, originally hosted at https://github.com/HongzhuoChen/FaceKit . Its SHA-256 is `4dbdd486cb942f93017c442dec7325936933a7d296f3ad99c66550fc54eaf377`. The accompanying MIT license is preserved from the same revision.

The upstream server no longer provides this historical commit or the application's recorded submodule commit. This snapshot makes the existing publication formula/dependency audit reproducible without changing the application submodule, upgrading the historical 125-feature schema, or running image extraction. The generator checks the exact source hash before its synthetic formula/dependency probe. The workflow independently checks `SHA256SUMS`.

The snapshot is for source auditing of frozen measurements. It does not establish per-image execution, detector-model checksums, or authorization to reprocess patient images. Do not edit this source to update current FaceKit behavior; the application's implementation remains in its separate submodule.
