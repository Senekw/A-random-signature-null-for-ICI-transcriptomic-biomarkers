"""A random-signature null for ICI transcriptomic biomarkers.

The re-runnable null-calibration harness accompanying Lee, "Most published
immunotherapy-response gene signatures do not outperform random gene sets,
and collapse onto a single tumor-inflammation axis".

The test the package exists to make cheap:

    >>> from icinull import load_cohort, random_null_test
    >>> c = load_cohort("Mariathasan")
    >>> random_null_test(c.expr, c.labels(), my_genes, universe=c.universe)

See ``docs/using_the_harness.md`` for evaluating a new signature.
"""

from .io import (  # noqa: F401
    Cohort,
    expressed_universe,
    list_cohorts,
    load_cohort,
    load_config,
    load_signatures,
    repo_root,
)
from .nullmodel import NullResult, random_null_test, stouffer  # noqa: F401
from .scoring import (  # noqa: F401
    auroc,
    auroc_se_hanley_mcneil,
    logit,
    logit_se,
    mean_z_score,
    score_signature,
    ssgsea_score,
)

__version__ = "1.0.0"

__all__ = [
    "Cohort",
    "NullResult",
    "auroc",
    "auroc_se_hanley_mcneil",
    "expressed_universe",
    "list_cohorts",
    "load_cohort",
    "load_config",
    "load_signatures",
    "logit",
    "logit_se",
    "mean_z_score",
    "random_null_test",
    "repo_root",
    "score_signature",
    "ssgsea_score",
    "stouffer",
]
