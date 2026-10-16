---
title: API reference
---

# API reference

The package lives at `batchdx.cyto`. The public surface is organised by stage:
input and manifests (`io`), preprocessing (`preprocess`), the design audit
(`design`), diagnostics (`diagnostics`), correction backends (`adapters`), the
frozen sampling layer (`sampling`), metrics (`metrics`), the decision layer
(`decisions`), the report (`reporting`), and the schemas every stage exchanges
(`schemas`).

Every stage boundary is a contract: the schemas module defines what one stage may
hand to another, and the contracts are enforced, not conventional. Start with
[schemas](schemas.md) if you want to understand the data model, or with
[pipeline](pipeline.md) if you want to see the stages in order.