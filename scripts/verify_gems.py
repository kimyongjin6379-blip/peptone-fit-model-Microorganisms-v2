"""
10종 GEM 검증 스크립트

outputs/gem_cache/*.xml 파일을 모두 로드하고 기본 μ를 출력.
ALLmed/LBmed 조건(gap-fill 때 쓴 배지)의 이론적 최대 성장률.
"""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GEM_DIR = ROOT / "outputs" / "gem_cache"

# peptomatch 패키지 import 가능하도록
sys.path.insert(0, str(ROOT / "src"))

try:
    from peptomatch.fba_simulator import FBASimulator
except ImportError as e:
    print(f"[ERROR] peptomatch.fba_simulator import 실패: {e}")
    print("COBRApy 직접 사용으로 폴백합니다...\n")
    FBASimulator = None

if FBASimulator is None:
    import cobra
    cobra.Configuration().solver = "glpk"

    def predict(xml_path):
        model = cobra.io.read_sbml_model(str(xml_path))
        sol = model.optimize()
        return sol.objective_value, len(model.reactions), len(model.metabolites)

    print(f"{'Strain':<8}{'mu (h^-1)':<14}{'Reactions':<12}{'Metabolites':<12}")
    print("-" * 46)
    for xml in sorted(GEM_DIR.glob("*.xml")):
        try:
            mu, n_rxn, n_met = predict(xml)
            mu_str = f"{mu:.3f}" if mu is not None else "FAIL"
            print(f"{xml.stem:<8}{mu_str:<14}{n_rxn:<12}{n_met:<12}")
        except Exception as e:
            print(f"{xml.stem:<8}ERROR: {e}")
else:
    print(f"{'Strain':<8}{'mu (h^-1)':<14}")
    print("-" * 22)
    for xml in sorted(GEM_DIR.glob("*.xml")):
        try:
            sim = FBASimulator(str(xml))
            mu = sim.predict_growth()
            print(f"{xml.stem:<8}{mu:.3f}")
        except Exception as e:
            print(f"{xml.stem:<8}ERROR: {e}")
