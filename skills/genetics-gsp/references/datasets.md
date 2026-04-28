---
datasets:
  - name: IBD
    s3_pickle: s3://tec-dev-usvga-11158-agenticboost-data-results-01/gsp/IBD_20260407/output/GSP.pkl.gz
    s3_xlsx: s3://tec-dev-usvga-11158-agenticboost-data-results-01/gsp/IBD_20260407/output/GSP.xlsx
    aws_profile: cmp-dev
    local_dir: ~/tmp/IBD_20260407
    indications:
      - efo_id: EFO_0003767
        name: Inflammatory bowel disease
        aliases: IBD
      - efo_id: EFO_0000729
        name: Ulcerative colitis
        aliases: UC, colitis
      - efo_id: EFO_0000384
        name: Crohn's disease
        aliases: CD, Crohns

  - name: SSc
    s3_pickle: s3://tec-dev-usvga-11158-agenticboost-data-results-01/gsp/ILD_IPF_SSc_20260325/output/GSP.pkl.gz
    s3_xlsx: s3://tec-dev-usvga-11158-agenticboost-data-results-01/gsp/ILD_IPF_SSc_20260325/output/GSP.xlsx
    aws_profile: cmp-dev
    local_dir: ~/tmp/ILD_IPF_SSc_20260325
    indications:
      - efo_id: EFO_0000768
        name: Idiopathic pulmonary fibrosis
        aliases: IPF
      - efo_id: EFO_0000717
        name: Systemic sclerosis
        aliases: SSc, scleroderma
      - efo_id: EFO_0004244
        name: Interstitial lung disease
        aliases: ILD

  - name: HS
    s3_pickle: s3://tec-dev-usvga-11158-agenticboost-data-results-01/gsp/EFO_1000710_MONDO_0006559_20260326_194012/output/GSP.pkl.gz
    s3_xlsx: s3://tec-dev-usvga-11158-agenticboost-data-results-01/gsp/EFO_1000710_MONDO_0006559_20260326_194012/output/GSP.xlsx
    aws_profile: cmp-dev
    local_dir: ~/tmp/EFO_1000710_MONDO_0006559_20260326_194012
    indications:
      - efo_id: EFO_1000710
        name: Hidradenitis suppurativa

  - name: AtD
    s3_pickle: s3://tec-dev-usvga-11158-agenticboost-data-results-01/gsp/AtD_20260410/output/GSP.pkl.gz
    s3_xlsx: s3://tec-dev-usvga-11158-agenticboost-data-results-01/gsp/AtD_20260410/output/GSP.xlsx
    aws_profile: cmp-dev
    local_dir: ~/tmp/AtD_20260410
    indications:
      - efo_id: EFO_0000274
        name: Atopic dermatitis
        aliases: AtD, atopic eczema, eczema
      - efo_id: EFO_1000651
        name: Recalcitrant atopic dermatitis
---

# GSP Dataset Registry

This file maps dataset names to their S3 locations and available indications. The YAML frontmatter is machine-readable; the table below is for human reference.

| Dataset | Indications | EFO IDs |
|---------|-------------|---------|
| **IBD** | Inflammatory bowel disease, Ulcerative colitis, Crohn's disease | EFO_0003767, EFO_0000729, EFO_0000384 |
| **SSc** | Idiopathic pulmonary fibrosis, Systemic sclerosis, Interstitial lung disease | EFO_0000768, EFO_0000717, EFO_0004244 |
| **HS** | Hidradenitis suppurativa | EFO_1000710 |
| **AtD** | Atopic dermatitis, Recalcitrant atopic dermatitis | EFO_0000274, EFO_1000651 |
