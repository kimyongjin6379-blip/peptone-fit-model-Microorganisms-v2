# PeptoMatch

**Genome + GEM-driven Peptone Recommendation System**

NCBI 유전체 정보, gapseq으로 자동 생성한 GEM(genome-scale metabolic model), 그리고 Growth Curve App에서 누적되는 실측 성장 데이터를 결합해 미생물 균주에 최적의 펩톤을 추천하고 검증하는 통합 플랫폼.

---

## Overview

PeptoMatch는 세 가지 정보 소스를 결합합니다:

1. **유전체 기반 demand 예측** — NCBI에서 가져온 균주 genome을 KO annotation해서 AA/비타민/뉴클레오티드 생합성 결손을 점수화 (균주가 무엇을 외부에서 받아야 하는가)
2. **펩톤 조성 기반 supply 점수** — 자사 분석 데이터(53종 × 89 features) 기반 영양 공급 점수 (균주에게 무엇을 줄 수 있는가)
3. **GEM 기반 FBA 시뮬레이션** — gapseq으로 만든 22종 SBML 모델로 배지+펩톤 조합에서의 성장률(μ) 예측 + shadow price 기반 병목 분석

여기에 **Growth Curve App에서 자동 적재되는 실측 OD/μ 데이터**가 누적되면서 추천 정확도를 점진적으로 보정합니다.

### 통합 파이프라인

```
NCBI Genome ──┬── KO Annotation ───→ Demand Score ──┐
              │                                       │
              └── gapseq → SBML GEM ──→ FBA μ ───────┤
                                                      ├──→ Match Score → Ranking
Peptone Composition ──→ Supply Score ─────────────────┤
                                                      │
Growth Curve App ──→ POST /api/ingest ──→ Growth DB ──┘
                                              ↓
                                        ML Calibration
                                    (Ridge → GP → XGBoost)
```

---

## Features

### 추천 엔진 (Production)
- 펩톤 조성(53종) × 균주 demand 매칭 점수
- Sempio 자사 소재 16종 우선 추천 (전체 비교 옵션)
- 펩톤 블렌딩 최적화 (`blend_optimizer.py`)
- 한/영 추천 이유 자동 생성 + PDF 리포트

### GEM 기반 FBA 시뮬레이션 (Production)
- **22종 GEM 보유** (gapseq로 자동 생성, `outputs/gem_cache/`)
- 균주별 배지+펩톤 조성에서 성장률 예측 (`fba_simulator.py`)
- Shadow price 기반 병목 영양소 분석
- 배지 정의: MRS / BHI / LB / TSB (rich · strict · cysteine 변형, `data/gapseq_media/`)

### 성장 곡선 데이터 적재 (Production)
- Growth Curve App → `POST /api/ingest` 자동 수신
- SQLite 기반 `experiments` / `growth_curves` / `growth_metrics` 스키마 (`growth_db.py`)
- 팀 약어 자동 해석 (LR → *L. rhamnosus*, N → SOY-N+ 등)
- 블렌딩 표기 파싱 (`1_75%+P_25%`)
- MRS 컨트롤 row 식별 (`is_control`)

### ML 보정 (대기 — 데이터 누적 후 활성화)
- Ridge baseline → Gaussian Process → XGBoost 단계별 학습 (`ml_train.py`, `ml_predict.py`)
- Feature: composition × genome demand × FBA μ
- Active learning 지원 (GP 불확실성 기반 다음 실험 제안)
- FBA 예측 vs 실측 Spearman ρ 검증 (`calibration.py`)

### Web UI (FastAPI + Tailwind)
- `gateway.py` — 다크 테마, Growth Curve App / 3P Tool과 통일된 디자인
- 추천 결과 + 실측 OD 동시 표시 (데이터가 있을 때)
- Streamlit UI는 legacy로 `app/streamlit_app.py`에 유지

---

## Current GEM Inventory (22종)

| Tier | 약어 | 균주 | 배지 |
|------|------|------|------|
| Core | LP, LR, LA, LC, LPC, LS | *Lactobacillus* (acidophilus / casei / paracasei / plantarum / rhamnosus / salivarius) | MRS |
| Core | LB, LF, LG, LH, LRE | *L. delbrueckii bulgaricus / fermentum / gasseri / helveticus* / *Lm. reuteri* | MRS |
| Core | STT, EF | *S. thermophilus* / *E. faecalis* | MRS / BHI |
| Core | BS, BC | *Bacillus subtilis / coagulans* | LB / TSB |
| Core | EC | *E. coli* | LB |
| Tier 2 | LL, EFM | *Lactococcus lactis* / *E. faecium* | MRS / BHI |
| Tier 2 | BFB, BFR, BFL, BFA | *Bifidobacterium bifidum / breve / longum / animalis subsp. lactis* | MRS |

확장 시 `scripts/gapseq_batch.sh` 또는 균주 단위 `gapseq doall` 실행 (WSL/Linux 필요).

---

## Project Structure

```
peptomatch/
├── README.md
├── pyproject.toml / requirements.txt
├── nixpacks.toml / Procfile / start.sh / runtime.txt   # Railway 배포
├── gateway.py                       # FastAPI 메인 앱 (UI + API)
│
├── config/config.yaml
│
├── data/
│   ├── composition_template.xlsx    # 펩톤 53종 × 89 features
│   ├── basal_media_composition.xlsx # MRS / BHI / LB / TSB 조성 + KEGG ID
│   ├── 신사업1팀 균주 리스트 (2024 ver.).xlsx
│   └── gapseq_media/                # gapseq 입력용 배지 csv
│       ├── MRS_rich.csv, MRS_rich_cys.csv, MRS_strict.csv
│       ├── BHI_rich.csv, BHI_strict.csv
│       └── LB_rich.csv, TSB_rich.csv
│
├── outputs/
│   └── gem_cache/                   # SBML GEM (22종 .xml)
│
├── templates/ + static/             # FastAPI Jinja2 템플릿 + 정적 자산
├── app/                             # Streamlit (legacy)
├── models/                          # 학습된 ML 모델 저장소 (보정 모델용)
│
├── scripts/
│   ├── gapseq_batch.sh              # GEM 일괄 생성 (WSL)
│   ├── download_genomes.sh          # NCBI genome 다운로드
│   ├── verify_gems.py               # GEM 로드 + μ 검증
│   ├── setup_kofamscan.sh           # KofamScan 환경 셋업
│   └── README_GAPSEQ_BATCH.md
│
└── src/peptomatch/
    │
    ├── # Recommendation core
    ├── scoring.py / explain.py / blend_optimizer.py
    ├── composition_features.py / genome_prior.py
    ├── kegg_pathway.py / kegg_client.py / kegg_viz.py
    ├── ko_annotator.py / taxonomy_priors.py
    │
    ├── # Strain & data sources
    ├── strain_db.py / ncbi_client.py / io_loaders.py
    │
    ├── # GEM / FBA
    ├── gem_manager.py / fba_simulator.py / media_config.py
    │
    ├── # Growth data ingestion
    ├── growth_db.py
    │
    ├── # ML calibration (대기 — 데이터 누적 후 활성)
    ├── ml_features.py / ml_train.py / ml_predict.py
    ├── calibration.py
    │
    ├── # Reports & utilities
    ├── report_pdf.py / compare.py / utils.py / cli.py
```

---

## Installation

### Requirements
- Python ≥ 3.11
- (Optional) [KofamScan](https://www.genome.jp/tools/kofamkoala/) — 정밀 KO annotation
- (Optional, 새 GEM 생성 시) WSL2 / Linux + [gapseq](https://github.com/jotech/gapseq)

### Setup

```bash
# 의존성 설치
pip install -e .
# 또는
pip install -r requirements.txt
```

`cobra` (FBA용)는 Linux/Mac에서 가장 안정적. Windows에서도 동작은 하나 일부 솔버 제약 있음.

---

## Usage

### Web UI (FastAPI)

```bash
uvicorn gateway:app --reload --host 0.0.0.0 --port 8000
```

- `/` — 대시보드
- `/recommend` — 펩톤 추천 폼 + 결과
- `/growth` — 적재된 실험 데이터 브라우저

### CLI

```bash
# 균주 목록 확인
peptomatch list-strains

# 펩톤 추천 (Sempio 소재만, 기본)
peptomatch recommend --strain "plantarum" --topk 5

# 전체 53종 비교
peptomatch recommend --strain "plantarum" --topk 10 --all-peptones

# 영문 출력
peptomatch recommend --strain "Bifidobacterium" --topk 5 --language en

# CSV 내보내기
peptomatch recommend --strain "subtilis" --format csv

# NCBI에서 균주 검색 및 추가
peptomatch search-ncbi --taxon "Lactobacillaceae" --limit 10 --add

# 유전체 prior 빌드
peptomatch build-priors
```

### Streamlit UI (legacy)

```bash
streamlit run app/streamlit_app.py
```

---

## API Endpoints

### Growth Data API (Growth Curve App에서 호출)

| Method | Path | 용도 |
|--------|------|------|
| `POST` | `/api/ingest` | 처리된 실험 데이터(JSON) 적재 |
| `GET` | `/api/growth/summary` | 적재 요약 |
| `GET` | `/api/growth/experiments` | 실험 목록 |
| `GET` | `/api/growth/curves` | 개별 성장 곡선 |
| `GET` | `/api/growth/ml-data` | ML 학습용 정형 데이터 |
| `GET` | `/api/growth/fba-data` | FBA 검증용 (균주, 펩톤, μ_max) |

### Recommendation API

| Method | Path | 용도 |
|--------|------|------|
| `POST` | `/api/recommend` | JSON in/out 추천 |

### Health

| Method | Path |
|--------|------|
| `GET` | `/healthz` / `/api/health` |

---

## Adding New GEMs (WSL / Linux)

```bash
# 1. genome 다운로드 (NCBI datasets CLI)
mkdir -p gapseq_work/<CODE> && cd gapseq_work/<CODE>
datasets download genome accession <GCF_ACCESSION> --include genome
unzip -o ncbi_dataset.zip
cp $(find ncbi_dataset -name "*.fna" | head -1) <CODE>.fna

# 2. gapseq 실행 (배지 = 실험 프로토콜과 일치)
gapseq doall \
  -m /mnt/d/folder1/peptomatch/data/gapseq_media/MRS_rich.csv \
  -v 1 \
  <CODE>.fna 2>&1 | tee <CODE>.log

# 3. peptomatch로 복사
cp <CODE>.xml /mnt/d/folder1/peptomatch/outputs/gem_cache/

# 4. 검증
python scripts/verify_gems.py
```

균주 그룹별 권장 배지:

| 균주 그룹 | 권장 gapseq media |
|-----------|-------------------|
| Lactobacillus / Lactococcus / Streptococcus | `MRS_rich.csv` |
| Bifidobacterium | `MRS_rich.csv` (실험 프로토콜이 MRS면) |
| Enterococcus | `BHI_rich.csv` |
| Bacillus | `LB_rich.csv` 또는 `TSB_rich.csv` |
| E. coli, Pseudomonas | `LB_rich.csv` |

---

## Configuration

`config/config.yaml` 주요 설정:

| 항목 | 설명 |
|------|------|
| `data.composition_file` | 펩톤 조성 Excel 경로 |
| `data.strain_file` | 균주 리스트 Excel 경로 |
| `peptone_filter` | 추천 후보 펩톤 목록 (Sempio 소재 16종) |
| `weights.*` | 매칭 스코어 가중치 |
| `ncbi.email` | NCBI API 사용자 이메일 |
| `annotation.strategy` | KO annotation 전략 (auto / kofamscan / gff3 / taxonomy) |
| `annotation.kofamscan_path` | KofamScan 실행파일 경로 |

### Sempio 소재 (기본 필터)

| 카테고리 | 제품 |
|----------|------|
| SOY 시리즈 | SOY-1, SOY-N+, SOY-L, SOY-P, SOY-BIO, SOY-BIO N50 |
| WHEAT / PEA / RICE | WHEAT-1, WHEAT-BIO, PEA-1, PEA-BIO, RICE-1, RICE-BIO |
| 동물성 | PPR Type2, PPR Type3, PPR Type4, Fish Collagen |

---

## Team Aliases (DB 자동 해석)

### 균주 약어 (`STRAIN_ALIASES` in `growth_db.py`)

| 약어 | 균주 | 약어 | 균주 |
|------|------|------|------|
| LP | *L. plantarum* | LRE | *Lm. reuteri* |
| LR | *L. rhamnosus* | LB | *L. delbrueckii* subsp. *bulgaricus* |
| LA | *L. acidophilus* | LF | *Lm. fermentum* |
| LC | *L. casei* | LH | *L. helveticus* |
| LPC | *L. paracasei* | LG | *L. gasseri* |
| LS | *L. salivarius* | EF | *E. faecalis* |
| STT | *S. thermophilus* | BS | *B. subtilis* |
| BC | *B. coagulans* | EC | *E. coli* |

(추가 균주는 `growth_db.py`의 `STRAIN_ALIASES` 딕셔너리에 등록)

### 펩톤 약어 (`PEPTONE_ALIASES`)

| 약어 | 정식명 | 약어 | 정식명 |
|------|--------|------|--------|
| 1 | SOY-1 | W | WHEAT-1 |
| N | SOY-N+ | R | RICE-1 |
| L | SOY-L | P | PEA-1 |
| B | SOY-B | PP | PPR Type4 |
| SP | SOY-P | | |

### 블렌딩 표기

```
형식: {약어}_{비율}%+{약어}_{비율}%
예시: 1_75%+P_25%   →  SOY-1 75% + PEA-1 25%
```

---

## Scoring Algorithm

매칭 스코어는 Supply(펩톤 조성) × Demand(균주 영양 요구)의 가중합으로 계산:

1. **FAA / TAA abundance** — 유리/총 아미노산 함량
2. **MW distribution** — 분자량 분포 × 펩타이드 수송체 활성
3. **AA-specific matching** — 필수 아미노산별 supply-demand 매칭
4. **Vitamin / Nucleotide** — 비타민·뉴클레오티드 공급-수요 매칭
5. **Transporter bonus** — 펩타이드 수송체 보유 시 추가 점수

향후 확장 (데이터 누적 후):
- **FBA-predicted μ** 가중치 추가
- **ML-corrected score** 레이어 (Ridge → GP → XGBoost)

---

## Data Sources

- **펩톤 조성**: 자사 분석 데이터 (89 features — FAA, TAA, MW, 비타민, 미네랄, 뉴클레오티드)
- **균주 유전체**: [NCBI Assembly Database](https://www.ncbi.nlm.nih.gov/assembly/)
- **GEM 생성**: [gapseq](https://github.com/jotech/gapseq) (NCBI genome → SBML)
- **FBA 엔진**: [COBRApy](https://opencobra.github.io/cobrapy/)
- **경로 정의**: [KEGG Pathway](https://www.kegg.jp/kegg/pathway.html), [ModelSEED Biochem](https://modelseed.org/)
- **KO annotation**: [KofamScan](https://www.genome.jp/tools/kofamkoala/) / GFF3 product matching

---

## Deployment (Railway)

```
nixpacks.toml  →  Python 3.11 + hmmer + ruby + gcc + sqlite + wget + gzip
start.sh       →  uvicorn gateway:app
Volume         →  /app/data (DB 영속화)
```

배포 환경 변수:
- `KOFAM_DIR` (선택, KofamScan 사용 시)
- `KOFAMSCAN_PATH` / `KOFAMSCAN_PROFILES`

---

## Roadmap

| 단계 | 상태 |
|------|------|
| Phase 1: Growth DB + ingestion API + 약어 해석 | ✅ 완료 |
| Phase 2: FastAPI UI + 성장 데이터 브라우저 | ✅ 완료 |
| Phase 3a: gapseq GEM 생성 (22종) | ✅ 완료 |
| Phase 3b: FBA 시뮬레이터 + 배지 변환 | ✅ 완료 |
| Phase 4a: FBA vs 실측 Spearman ρ 검증 | ⏳ 데이터 누적 후 |
| Phase 4b: ML 보정 모델 학습 | ⏳ 데이터 누적 후 |
| Phase 5: 배지 최적화 + 통합 추천 | 📋 장기 |

목표: **R² ≥ 0.7**, 실측 OD **15%↑** (vs 상업적 기본 제품)

---

## License

Internal use only — Sempio Foods Company
