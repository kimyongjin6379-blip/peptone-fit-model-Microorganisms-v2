# PeptoMatch

**Genome-driven Peptone Recommendation System**

NCBI 유전체 데이터 기반으로 미생물 균주에 최적의 펩톤을 추천하는 시스템입니다.

## Overview

PeptoMatch는 미생물의 유전체 정보(아미노산 생합성 경로, 비타민 합성, 수송체 등)를 분석하여 해당 균주의 영양 요구량(demand)을 예측하고, 펩톤의 조성 데이터(supply)와 매칭하여 최적의 펩톤을 추천합니다.

### 핵심 파이프라인

```
Genome (NCBI) → KO Annotation → Pathway Completeness → Demand Score
                                                            ↓
Peptone Composition → Feature Extraction → Supply Score → Match Score → Ranking
```

### 주요 기능

- **유전체 기반 추천**: NCBI Assembly DB에서 유전체를 가져와 KO annotation 수행
- **다중 annotation 전략**: KofamScan (primary) → GFF3 product matching (fallback) → Taxonomy prior (fallback)
- **Sempio 소재 필터**: 자사 펩톤 16종만 선별하여 추천 (전체 53종 비교도 가능)
- **균주 DB 확장**: NCBI에서 신규 균주 검색 및 추가
- **Web UI**: Streamlit 기반 대화형 인터페이스

## Project Structure

```
peptomatch/
├── README.md
├── pyproject.toml
├── requirements.txt
├── .gitignore
├── config/
│   └── config.yaml              # 설정 파일 (경로, 가중치, NCBI, annotation)
├── data/
│   ├── composition_template.xlsx # 펩톤 조성 데이터 (53종 × 89 features)
│   └── 신사업1팀 균주 리스트 (2024 ver.).xlsx  # 보유 균주 목록 (54종)
├── app/
│   └── streamlit_app.py         # Streamlit Web UI
└── src/
    └── peptomatch/
        ├── __init__.py
        ├── utils.py             # 설정 로드, 경로 처리, 데이터 클리닝
        ├── io_loaders.py        # Excel 데이터 로딩 (복합 헤더 파싱)
        ├── composition_features.py  # 펩톤 조성 → supply score 변환
        ├── kegg_pathway.py      # KEGG KO 경로 정의 (AA, 비타민, 뉴클레오타이드)
        ├── taxonomy_priors.py   # 30+ genus에 대한 문헌 기반 기본값
        ├── ncbi_client.py       # NCBI Datasets API v2 + Entrez client
        ├── ko_annotator.py      # KofamScan + GFF3 기반 KO annotation
        ├── genome_prior.py      # 유전체 → demand score 변환
        ├── strain_db.py         # SQLite 균주 DB (NCBI 확장 지원)
        ├── scoring.py           # Supply-Demand 매칭 및 추천 엔진
        ├── explain.py           # 추천 이유 생성 (한/영)
        └── cli.py               # CLI 인터페이스
```

## Installation

### Requirements

- Python >= 3.11
- (Optional) [KofamScan](https://www.genome.jp/tools/kofamkoala/) - 정밀 KO annotation용

### Setup

```bash
# 의존성 설치
pip install -e .

# 또는
pip install -r requirements.txt
```

## Usage

### CLI

```bash
# 균주 목록 확인
peptomatch list-strains

# 펩톤 추천 (Sempio 소재만, 기본)
peptomatch recommend --strain "plantarum" --topk 5

# 펩톤 추천 (전체 53종)
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

### Streamlit Web UI

```bash
streamlit run app/streamlit_app.py
```

3개 탭으로 구성:
1. **펩톤 추천**: 균주 선택 → 분석 → 추천 결과 (차트 + 테이블 + CSV 다운로드)
2. **균주 브라우저**: 보유 균주 목록 조회 + NCBI 검색/추가
3. **펩톤 탐색기**: 펩톤 간 성분 비교 (FAA 차트 + 요약 테이블)

## Configuration

`config/config.yaml` 주요 설정:

| 항목 | 설명 |
|------|------|
| `data.composition_file` | 펩톤 조성 Excel 경로 |
| `data.strain_file` | 균주 리스트 Excel 경로 |
| `peptone_filter` | 추천 후보 펩톤 목록 (Sempio 소재 16종) |
| `weights.*` | 매칭 스코어 가중치 |
| `ncbi.email` | NCBI API 사용자 이메일 |
| `annotation.strategy` | KO annotation 전략 (auto/kofamscan/gff3/taxonomy) |
| `annotation.kofamscan_path` | KofamScan 실행파일 경로 (설치 시) |

### Sempio 소재 (기본 필터)

| 카테고리 | 제품 |
|----------|------|
| SOY 시리즈 | SOY-1, SOY-N+, SOY-L, SOY-P, SOY-BIO, SOY-BIO N50 |
| WHEAT/PEA/RICE | WHEAT-1, WHEAT-BIO, PEA-1, PEA-BIO, RICE-1, RICE-BIO |
| 동물성 | PPR Type2, PPR Type3, PPR Type4, Fish Collagen |

## Scoring Algorithm

매칭 스코어는 Supply(펩톤 조성) × Demand(균주 영양 요구)의 가중합으로 계산됩니다:

1. **FAA/TAA abundance**: 유리/총 아미노산 함량
2. **MW distribution**: 분자량 분포 (저분자 선호도 × 수송체 활성)
3. **AA-specific matching**: 필수 아미노산별 supply-demand 매칭
4. **Vitamin/Nucleotide**: 비타민·뉴클레오타이드 공급-수요 매칭
5. **Transporter bonus**: 펩타이드 수송체 보유 시 추가 점수

## Data Sources

- **펩톤 조성**: 자사 분석 데이터 (89 features: FAA, TAA, MW, 비타민, 미네랄, 뉴클레오타이드)
- **균주 유전체**: [NCBI Assembly Database](https://www.ncbi.nlm.nih.gov/assembly/)
- **경로 정의**: [KEGG Pathway](https://www.kegg.jp/kegg/pathway.html) (참조 지식으로만 사용)
- **KO annotation**: [KofamScan](https://www.genome.jp/tools/kofamkoala/) / GFF3 product matching

## KofamScan Setup (Optional)

KofamScan을 설치하면 보다 정확한 KO annotation이 가능합니다:

```bash
# KofamScan 설치 (Linux/Mac)
# https://www.genome.jp/tools/kofamkoala/ 에서 다운로드

# HMM 프로파일 다운로드 (~4GB)
wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
tar xzf profiles.tar.gz

# config.yaml에 경로 설정
# annotation:
#   kofamscan_path: "/path/to/kofamscan"
#   kofamscan_profiles: "/path/to/profiles"
```

KofamScan 미설치 시 GFF3 product matching → Taxonomy prior 순으로 자동 fallback됩니다.

## License

Internal use only - Sempio Foods Company
