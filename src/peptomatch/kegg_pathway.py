"""KEGG pathway analysis for biosynthesis completeness calculation."""

import logging
from typing import Any, Optional

logger = logging.getLogger("peptone_recommender")


# KEGG KO (KEGG Orthology) IDs for amino acid biosynthesis pathways
# Reference: https://www.genome.jp/kegg/pathway.html

AA_BIOSYNTHESIS_KOS = {
    # Essential amino acids
    "His": {
        "name": "Histidine biosynthesis",
        "pathway": "M00026",
        "kos": [
            "K00765",  # hisG; ATP phosphoribosyltransferase
            "K02501",  # hisI; phosphoribosyl-AMP cyclohydrolase
            "K01814",  # hisA; phosphoribosylformimino-5-aminoimidazole carboxamide ribotide isomerase
            "K02500",  # hisH; imidazole glycerol-phosphate synthase
            "K01693",  # hisB; imidazoleglycerol-phosphate dehydratase
            "K00013",  # hisD; histidinol dehydrogenase
            "K01089",  # hisC; histidinol-phosphate aminotransferase
        ],
        "essential_kos": ["K00765", "K00013"],  # Key enzymes
    },
    "Ile": {
        "name": "Isoleucine biosynthesis",
        "pathway": "M00019",
        "kos": [
            "K01754",  # ilvA; threonine dehydratase
            "K01652",  # ilvB; acetolactate synthase I/II/III large subunit
            "K01653",  # ilvH; acetolactate synthase small subunit
            "K00053",  # ilvC; ketol-acid reductoisomerase
            "K01687",  # ilvD; dihydroxy-acid dehydratase
            "K00826",  # ilvE; branched-chain amino acid aminotransferase
        ],
        "essential_kos": ["K01754", "K01652", "K00826"],
    },
    "Leu": {
        "name": "Leucine biosynthesis",
        "pathway": "M00019",
        "kos": [
            "K01652",  # ilvB; acetolactate synthase
            "K00053",  # ilvC; ketol-acid reductoisomerase
            "K01687",  # ilvD; dihydroxy-acid dehydratase
            "K01649",  # leuA; 2-isopropylmalate synthase
            "K01703",  # leuC; 3-isopropylmalate dehydratase large subunit
            "K01704",  # leuD; 3-isopropylmalate dehydratase small subunit
            "K00052",  # leuB; 3-isopropylmalate dehydrogenase
            "K00826",  # ilvE; branched-chain amino acid aminotransferase
        ],
        "essential_kos": ["K01649", "K00052", "K00826"],
    },
    "Lys": {
        "name": "Lysine biosynthesis (DAP pathway)",
        "pathway": "M00016",
        "kos": [
            "K00928",  # lysC; aspartate kinase
            "K00133",  # asd; aspartate-semialdehyde dehydrogenase
            "K01714",  # dapA; dihydrodipicolinate synthase
            "K00215",  # dapB; dihydrodipicolinate reductase
            "K01439",  # dapE; succinyl-diaminopimelate desuccinylase
            "K01778",  # dapF; diaminopimelate epimerase
            "K01586",  # lysA; diaminopimelate decarboxylase
        ],
        "essential_kos": ["K01714", "K01586"],
    },
    "Met": {
        "name": "Methionine biosynthesis",
        "pathway": "M00017",
        "kos": [
            "K00928",  # lysC; aspartate kinase
            "K00133",  # asd; aspartate-semialdehyde dehydrogenase
            "K00003",  # hom; homoserine dehydrogenase
            "K00651",  # metA; homoserine O-succinyltransferase
            "K01739",  # metB; cystathionine gamma-synthase
            "K01760",  # metC; cystathionine beta-lyase
            "K00548",  # metH; methionine synthase
            "K00549",  # metE; 5-methyltetrahydropteroyltriglutamate--homocysteine methyltransferase
        ],
        "essential_kos": ["K00651", "K00548", "K00549"],
    },
    "Phe": {
        "name": "Phenylalanine biosynthesis",
        "pathway": "M00024",
        "kos": [
            "K01609",  # trpC; indole-3-glycerol phosphate synthase (shared)
            "K00800",  # aroA; 3-phosphoshikimate 1-carboxyvinyltransferase
            "K01736",  # aroC; chorismate synthase
            "K04518",  # pheA2; prephenate dehydratase
            "K00832",  # tyrB; aromatic-amino-acid transaminase
            "K01850",  # pheA; chorismate mutase
        ],
        "essential_kos": ["K01736", "K04518", "K00832"],
    },
    "Thr": {
        "name": "Threonine biosynthesis",
        "pathway": "M00018",
        "kos": [
            "K00928",  # lysC; aspartate kinase
            "K00133",  # asd; aspartate-semialdehyde dehydrogenase
            "K00003",  # hom; homoserine dehydrogenase
            "K00872",  # thrB; homoserine kinase
            "K01733",  # thrC; threonine synthase
        ],
        "essential_kos": ["K00872", "K01733"],
    },
    "Trp": {
        "name": "Tryptophan biosynthesis",
        "pathway": "M00023",
        "kos": [
            "K01657",  # trpE; anthranilate synthase component I
            "K01658",  # trpG; anthranilate synthase component II
            "K00766",  # trpD; anthranilate phosphoribosyltransferase
            "K01817",  # trpF; phosphoribosylanthranilate isomerase
            "K01609",  # trpC; indole-3-glycerol phosphate synthase
            "K01695",  # trpA; tryptophan synthase alpha chain
            "K01696",  # trpB; tryptophan synthase beta chain
        ],
        "essential_kos": ["K01657", "K01695", "K01696"],
    },
    "Val": {
        "name": "Valine biosynthesis",
        "pathway": "M00019",
        "kos": [
            "K01652",  # ilvB; acetolactate synthase I/II/III large subunit
            "K01653",  # ilvH; acetolactate synthase small subunit
            "K00053",  # ilvC; ketol-acid reductoisomerase
            "K01687",  # ilvD; dihydroxy-acid dehydratase
            "K00826",  # ilvE; branched-chain amino acid aminotransferase
        ],
        "essential_kos": ["K01652", "K00053", "K00826"],
    },
    # Non-essential but important
    "Ala": {
        "name": "Alanine biosynthesis",
        "pathway": "M00030",
        "kos": [
            "K00814",  # GPT; alanine transaminase
            "K00259",  # ald; alanine dehydrogenase
        ],
        "essential_kos": ["K00814"],
    },
    "Arg": {
        "name": "Arginine biosynthesis",
        "pathway": "M00028",
        "kos": [
            "K00611",  # argF; ornithine carbamoyltransferase
            "K01940",  # argG; argininosuccinate synthase
            "K01755",  # argH; argininosuccinate lyase
            "K00930",  # argB; acetylglutamate kinase
            "K00145",  # argC; N-acetyl-gamma-glutamyl-phosphate reductase
            "K00821",  # argD; acetylornithine aminotransferase
        ],
        "essential_kos": ["K00611", "K01940", "K01755"],
    },
    "Asn": {
        "name": "Asparagine biosynthesis",
        "pathway": "M00029",
        "kos": [
            "K01953",  # asnB; asparagine synthase
            "K01914",  # asnA; aspartate--ammonia ligase
        ],
        "essential_kos": ["K01953"],
    },
    "Asp": {
        "name": "Aspartate biosynthesis",
        "pathway": "M00029",
        "kos": [
            "K00813",  # aspC; aspartate aminotransferase
            "K00812",  # aspB; aspartate aminotransferase
        ],
        "essential_kos": ["K00813"],
    },
    "Cys": {
        "name": "Cysteine biosynthesis",
        "pathway": "M00021",
        "kos": [
            "K00640",  # cysE; serine O-acetyltransferase
            "K01738",  # cysK; cysteine synthase A
            "K12339",  # cysM; cysteine synthase B
        ],
        "essential_kos": ["K00640", "K01738"],
    },
    "Glu": {
        "name": "Glutamate biosynthesis",
        "pathway": "M00027",
        "kos": [
            "K00262",  # gdhA; glutamate dehydrogenase (NADP+)
            "K00265",  # gltB; glutamate synthase large subunit
            "K00266",  # gltD; glutamate synthase small subunit
            "K01915",  # glnA; glutamine synthetase
        ],
        "essential_kos": ["K00262", "K01915"],
    },
    "Gln": {
        "name": "Glutamine biosynthesis",
        "pathway": "M00027",
        "kos": [
            "K01915",  # glnA; glutamine synthetase
        ],
        "essential_kos": ["K01915"],
    },
    "Gly": {
        "name": "Glycine biosynthesis",
        "pathway": "M00020",
        "kos": [
            "K00600",  # glyA; serine hydroxymethyltransferase
            "K00281",  # gcvP; glycine dehydrogenase
        ],
        "essential_kos": ["K00600"],
    },
    "Pro": {
        "name": "Proline biosynthesis",
        "pathway": "M00015",
        "kos": [
            "K00931",  # proB; glutamate 5-kinase
            "K00147",  # proA; glutamate-5-semialdehyde dehydrogenase
            "K00286",  # proC; pyrroline-5-carboxylate reductase
        ],
        "essential_kos": ["K00931", "K00286"],
    },
    "Ser": {
        "name": "Serine biosynthesis",
        "pathway": "M00020",
        "kos": [
            "K00058",  # serA; D-3-phosphoglycerate dehydrogenase
            "K00831",  # serC; phosphoserine aminotransferase
            "K01079",  # serB; phosphoserine phosphatase
        ],
        "essential_kos": ["K00058", "K01079"],
    },
    "Tyr": {
        "name": "Tyrosine biosynthesis",
        "pathway": "M00025",
        "kos": [
            "K00800",  # aroA; 3-phosphoshikimate 1-carboxyvinyltransferase
            "K01736",  # aroC; chorismate synthase
            "K00220",  # tyrA; prephenate dehydrogenase
            "K00832",  # tyrB; aromatic-amino-acid transaminase
        ],
        "essential_kos": ["K00220", "K00832"],
    },
}

# Vitamin B biosynthesis KOs
VITAMIN_BIOSYNTHESIS_KOS = {
    "B1": {  # Thiamine
        "name": "Thiamine biosynthesis",
        "pathway": "M00127",
        "kos": [
            "K00941",  # thiD; phosphomethylpyrimidine kinase
            "K00788",  # thiE; thiamine-phosphate pyrophosphorylase
            "K00946",  # thiL; thiamine-monophosphate kinase
            "K03147",  # thiC; phosphomethylpyrimidine synthase
        ],
        "essential_kos": ["K00788", "K00946"],
    },
    "B2": {  # Riboflavin
        "name": "Riboflavin biosynthesis",
        "pathway": "M00125",
        "kos": [
            "K00794",  # ribH; riboflavin synthase
            "K00793",  # ribE; riboflavin synthase alpha chain
            "K02858",  # ribD; diaminohydroxyphosphoribosylaminopyrimidine deaminase
            "K00082",  # ribA; GTP cyclohydrolase II
        ],
        "essential_kos": ["K00794", "K00082"],
    },
    "B3": {  # Niacin/NAD
        "name": "NAD biosynthesis",
        "pathway": "M00115",
        "kos": [
            "K00767",  # nadC; nicotinate-nucleotide pyrophosphorylase
            "K00969",  # nadD; nicotinate-nucleotide adenylyltransferase
            "K01916",  # nadE; NAD+ synthase
            "K03517",  # nadB; L-aspartate oxidase
        ],
        "essential_kos": ["K00767", "K01916"],
    },
    "B5": {  # Pantothenate
        "name": "Pantothenate biosynthesis",
        "pathway": "M00119",
        "kos": [
            "K00077",  # panE; 2-dehydropantoate 2-reductase
            "K00867",  # panK; pantothenate kinase
            "K01918",  # panC; pantoate--beta-alanine ligase
        ],
        "essential_kos": ["K00867", "K01918"],
    },
    "B6": {  # Pyridoxine
        "name": "Pyridoxine biosynthesis",
        "pathway": "M00124",
        "kos": [
            "K00275",  # pdxH; pyridoxamine 5'-phosphate oxidase
            "K03472",  # pdxJ; pyridoxine 5-phosphate synthase
            "K03473",  # pdxA; 4-hydroxythreonine-4-phosphate dehydrogenase
        ],
        "essential_kos": ["K00275", "K03472"],
    },
    "B7": {  # Biotin
        "name": "Biotin biosynthesis",
        "pathway": "M00123",
        "kos": [
            "K00652",  # bioF; 8-amino-7-oxononanoate synthase
            "K00833",  # bioA; adenosylmethionine-8-amino-7-oxononanoate transaminase
            "K01012",  # bioB; biotin synthase
            "K00796",  # bioD; dethiobiotin synthetase
        ],
        "essential_kos": ["K00652", "K01012"],
    },
    "B9": {  # Folate
        "name": "Folate biosynthesis",
        "pathway": "M00126",
        "kos": [
            "K00796",  # folK; 2-amino-4-hydroxy-6-hydroxymethyldihydropteridine pyrophosphokinase
            "K00950",  # folC; dihydrofolate synthase
            "K01930",  # folE; GTP cyclohydrolase I
            "K00287",  # folA; dihydrofolate reductase
        ],
        "essential_kos": ["K00796", "K00287"],
    },
    "B12": {  # Cobalamin
        "name": "Cobalamin biosynthesis",
        "pathway": "M00122",
        "kos": [
            "K02232",  # cobA; uroporphyrinogen III methyltransferase
            "K02224",  # cobI; precorrin-2 C20-methyltransferase
            "K02229",  # cobO; cob(I)alamin adenosyltransferase
        ],
        "essential_kos": ["K02232", "K02229"],
    },
}

# Nucleotide biosynthesis KOs
NUCLEOTIDE_BIOSYNTHESIS_KOS = {
    "purine": {
        "name": "Purine biosynthesis",
        "pathway": "M00048",
        "kos": [
            "K00764",  # purF; amidophosphoribosyltransferase
            "K01945",  # purD; phosphoribosylamine--glycine ligase
            "K01587",  # purN; phosphoribosylglycinamide formyltransferase
            "K01952",  # purL; phosphoribosylformylglycinamidine synthase
            "K01933",  # purM; phosphoribosylformylglycinamidine cyclo-ligase
            "K00602",  # purH; phosphoribosylaminoimidazolecarboxamide formyltransferase
            "K01756",  # purB; adenylosuccinate lyase
        ],
        "essential_kos": ["K00764", "K01952"],
    },
    "pyrimidine": {
        "name": "Pyrimidine biosynthesis",
        "pathway": "M00051",
        "kos": [
            "K01955",  # carB; carbamoyl-phosphate synthase large subunit
            "K00609",  # pyrB; aspartate carbamoyltransferase
            "K01465",  # pyrC; dihydroorotase
            "K00254",  # pyrD; dihydroorotate dehydrogenase
            "K01591",  # pyrF; orotidine-5'-phosphate decarboxylase
        ],
        "essential_kos": ["K01955", "K00609"],
    },
}

# Peptide/amino acid transporter KOs
TRANSPORTER_KOS = {
    "oligopeptide": {
        "name": "Oligopeptide transport system",
        "kos": [
            "K15580",  # oppA; oligopeptide transport system substrate-binding protein
            "K15581",  # oppB; oligopeptide transport system permease protein
            "K15582",  # oppC; oligopeptide transport system permease protein
            "K15583",  # oppD; oligopeptide transport system ATP-binding protein
            "K15584",  # oppF; oligopeptide transport system ATP-binding protein
        ],
    },
    "dipeptide": {
        "name": "Dipeptide transport system",
        "kos": [
            "K02035",  # dppA; dipeptide transport system substrate-binding protein
            "K02036",  # dppB; dipeptide transport system permease protein
            "K02037",  # dppC; dipeptide transport system permease protein
            "K02038",  # dppD; dipeptide transport system ATP-binding protein
            "K02039",  # dppF; dipeptide transport system ATP-binding protein
        ],
    },
    "amino_acid_general": {
        "name": "General amino acid transport",
        "kos": [
            "K09969",  # aapJ; general L-amino acid transport system substrate-binding protein
            "K09970",  # aapQ; general L-amino acid transport system permease protein
            "K09971",  # aapM; general L-amino acid transport system permease protein
            "K09972",  # aapP; general L-amino acid transport system ATP-binding protein
        ],
    },
    "branched_chain": {
        "name": "Branched-chain amino acid transport",
        "kos": [
            "K01999",  # livK; branched-chain amino acid transport system substrate-binding protein
            "K01997",  # livH; branched-chain amino acid transport system permease protein
            "K01998",  # livM; branched-chain amino acid transport system permease protein
            "K01995",  # livG; branched-chain amino acid transport system ATP-binding protein
            "K01996",  # livF; branched-chain amino acid transport system ATP-binding protein
        ],
    },
}


def calculate_pathway_completeness(
    ko_list: list[str],
    pathway_kos: dict[str, Any],
    use_essential_only: bool = False
) -> float:
    """Calculate completeness of a biosynthesis pathway.

    Supports two modes:
    - essential_ko_groups (OR groups): each group is an alternative way to
      fulfill one function. Score = fraction of groups with at least one hit.
      Example: glucose transport has multiple alternatives (ptsG, galP, glcU).
    - essential_kos (legacy AND mode): all listed KOs must be present.

    Args:
        ko_list: List of KO IDs found in the genome
        pathway_kos: Dictionary with 'kos' and optionally
            'essential_ko_groups' or 'essential_kos' keys
        use_essential_only: If True, check essential KOs/groups

    Returns:
        Completeness score (0.0 to 1.0)
    """
    # OR-group mode: each group represents one function with alternatives
    if use_essential_only and "essential_ko_groups" in pathway_kos:
        groups = pathway_kos["essential_ko_groups"]
        if not groups:
            return 1.0
        satisfied = sum(
            1 for group in groups
            if any(ko in ko_list for ko in group)
        )
        return satisfied / len(groups)

    # Legacy AND mode
    if use_essential_only and "essential_kos" in pathway_kos:
        target_kos = pathway_kos["essential_kos"]
    else:
        target_kos = pathway_kos["kos"]

    if not target_kos:
        return 1.0

    found = sum(1 for ko in target_kos if ko in ko_list)
    return found / len(target_kos)


def calculate_aa_biosynthesis(ko_list: list[str]) -> dict[str, float]:
    """Calculate amino acid biosynthesis completeness.

    Args:
        ko_list: List of KO IDs found in the genome

    Returns:
        Dictionary of AA -> completeness (0.0 to 1.0)
    """
    ko_set = set(ko_list)

    result = {}
    for aa, pathway_info in AA_BIOSYNTHESIS_KOS.items():
        completeness = calculate_pathway_completeness(
            ko_set, pathway_info, use_essential_only=True
        )
        result[aa] = completeness

    return result


def calculate_vitamin_biosynthesis(ko_list: list[str]) -> dict[str, float]:
    """Calculate vitamin biosynthesis completeness.

    Args:
        ko_list: List of KO IDs found in the genome

    Returns:
        Dictionary of vitamin -> completeness (0.0 to 1.0)
    """
    ko_set = set(ko_list)

    result = {}
    for vitamin, pathway_info in VITAMIN_BIOSYNTHESIS_KOS.items():
        completeness = calculate_pathway_completeness(
            ko_set, pathway_info, use_essential_only=True
        )
        result[vitamin] = completeness

    return result


def calculate_nucleotide_biosynthesis(ko_list: list[str]) -> float:
    """Calculate nucleotide biosynthesis completeness.

    Args:
        ko_list: List of KO IDs found in the genome

    Returns:
        Average completeness (0.0 to 1.0)
    """
    ko_set = set(ko_list)

    purine = calculate_pathway_completeness(
        ko_set, NUCLEOTIDE_BIOSYNTHESIS_KOS["purine"], use_essential_only=True
    )
    pyrimidine = calculate_pathway_completeness(
        ko_set, NUCLEOTIDE_BIOSYNTHESIS_KOS["pyrimidine"], use_essential_only=True
    )

    return (purine + pyrimidine) / 2


def calculate_transporter_score(ko_list: list[str]) -> float:
    """Calculate peptide/amino acid transporter score.

    Args:
        ko_list: List of KO IDs found in the genome

    Returns:
        Normalized transporter score (0.0 to 1.0)
    """
    ko_set = set(ko_list)

    total_found = 0
    total_possible = 0

    for transporter_info in TRANSPORTER_KOS.values():
        kos = transporter_info["kos"]
        total_possible += len(kos)
        total_found += sum(1 for ko in kos if ko in ko_set)

    if total_possible == 0:
        return 0.5  # Default

    return total_found / total_possible


def analyze_ko_annotations(ko_list: list[str]) -> dict[str, Any]:
    """Perform complete analysis of KO annotations.

    Args:
        ko_list: List of KO IDs found in the genome

    Returns:
        Complete analysis result
    """
    return {
        "aa_biosynthesis": calculate_aa_biosynthesis(ko_list),
        "vitamin_biosynthesis": calculate_vitamin_biosynthesis(ko_list),
        "nucleotide_biosynthesis": calculate_nucleotide_biosynthesis(ko_list),
        "transporter_score": calculate_transporter_score(ko_list),
        "sugar_metabolism": calculate_sugar_metabolism(ko_list),
        "organic_acid_metabolism": calculate_organic_acid_metabolism(ko_list),
        "mineral_transport": calculate_mineral_transport(ko_list),
        "total_kos": len(set(ko_list)),
        "source": "gcf_annotation"
    }


# Sugar metabolism / transport KOs
SUGAR_METABOLISM_KOS = {
    "glucose": {
        "name": "Glucose utilization",
        "kos": [
            "K02778",  # ptsG; glucose PTS system EIICB
            "K00844",  # glk; glucokinase
            "K00134",  # gapA; glyceraldehyde-3-phosphate dehydrogenase
            "K02777",  # crr; glucose PTS system EIIA
            "K00845",  # HK; hexokinase (eukaryotic)
            "K07025",  # galP; galactose:H+ symporter (alt glucose transporter)
            "K20118",  # glcU; glucose uptake protein
            "K02750",  # PTS-Glc-EIIA; glucose PTS system EIIA (alt)
        ],
        # OR groups: any ONE group satisfied → can use glucose
        "essential_ko_groups": [
            ["K02778"],           # ptsG (main PTS)
            ["K02777", "K02750"], # crr or alt PTS EIIA
            ["K07025"],           # galP (alternative transporter)
            ["K20118"],           # glcU (alternative transporter)
            ["K00844"],           # glucokinase (phosphorylation)
            ["K00845"],           # hexokinase (eukaryotic)
        ],
    },
    "fructose": {
        "name": "Fructose utilization",
        "kos": [
            "K02770",  # fruA; fructose PTS system EIIBC
            "K00882",  # fruK; 1-phosphofructokinase
            "K02768",  # fruB; fructose PTS system EIIA
            "K00847",  # scrK; fructokinase
        ],
        "essential_ko_groups": [
            ["K02770"],  # fruA (PTS)
            ["K00882"],  # fruK
            ["K00847"],  # fructokinase (alt)
        ],
    },
    "sucrose": {
        "name": "Sucrose utilization",
        "kos": [
            "K01193",  # sacA/INV; beta-fructofuranosidase (invertase)
            "K02810",  # scrA; sucrose PTS system EIIBCA
            "K00690",  # sacB; levansucrase
            "K01187",  # malZ; alpha-glucosidase
        ],
        "essential_ko_groups": [
            ["K01193"],  # invertase
            ["K02810"],  # sucrose PTS
            ["K00690"],  # levansucrase (alt)
        ],
    },
    "lactose": {
        "name": "Lactose utilization",
        "kos": [
            "K01220",  # lacZ; beta-galactosidase (class 1)
            "K02786",  # lacY; lactose permease (MFS transporter)
            "K01784",  # galE; UDP-glucose 4-epimerase
            "K00849",  # galK; galactokinase
            "K12111",  # lacA; beta-galactosidase (class 2)
            "K12308",  # bgaB; beta-galactosidase (class 3)
            "K02793",  # PTS-Lac-EIIA; lactose PTS
            "K02794",  # PTS-Lac-EIIB; lactose PTS
            "K02795",  # PTS-Lac-EIIC; lactose PTS
        ],
        "essential_ko_groups": [
            ["K01220"],          # lacZ
            ["K12111"],          # lacA (alt beta-gal)
            ["K12308"],          # bgaB (alt beta-gal)
            ["K02793"],          # lactose PTS EIIA
            ["K02786"],          # lacY permease
        ],
    },
    "maltose": {
        "name": "Maltose utilization",
        "kos": [
            "K01182",  # malZ; alpha-glucosidase
            "K10111",  # malE; maltose/maltodextrin transport protein
            "K10112",  # malF; maltose transport system permease
            "K10113",  # malG; maltose transport system permease
            "K01187",  # malL; oligo-1,6-glucosidase
            "K01176",  # amyA; alpha-amylase
        ],
        "essential_ko_groups": [
            ["K01182"],  # malZ
            ["K10111"],  # malE
            ["K01187"],  # malL (alt)
            ["K01176"],  # amylase (alt)
        ],
    },
}

# Organic acid metabolism KOs
ORGANIC_ACID_METABOLISM_KOS = {
    "lactate": {
        "name": "Lactate utilization/tolerance",
        "kos": [
            "K00016",  # ldh; L-lactate dehydrogenase
            "K07246",  # lctP; lactate permease (lldP)
            "K00101",  # lldD; L-lactate dehydrogenase (cytochrome)
            "K03778",  # dld; D-lactate dehydrogenase
        ],
        "essential_ko_groups": [
            ["K00016"],  # ldh (any lactate dehydrogenase)
            ["K00101"],  # lldD (alt)
            ["K03778"],  # dld (alt)
            ["K07246"],  # lctP permease
        ],
    },
    "citrate": {
        "name": "Citrate utilization",
        "kos": [
            "K01647",  # gltA; citrate synthase
            "K01643",  # citE; citrate lyase beta subunit
            "K01644",  # citF; citrate lyase alpha subunit
            "K16353",  # citP; citrate transporter (citM)
        ],
        "essential_ko_groups": [
            ["K01643"],  # citE
            ["K01644"],  # citF
            ["K16353"],  # citP transporter
        ],
    },
    "acetate": {
        "name": "Acetate utilization",
        "kos": [
            "K00925",  # ackA; acetate kinase
            "K01895",  # acs; acetyl-CoA synthetase
            "K00625",  # pta; phosphotransacetylase
        ],
        "essential_ko_groups": [
            ["K00925"],  # ackA
            ["K01895"],  # acs (alt)
            ["K00625"],  # pta
        ],
    },
    "succinate": {
        "name": "Succinate utilization",
        "kos": [
            "K00239",  # sdhA; succinate dehydrogenase flavoprotein subunit
            "K00240",  # sdhB; succinate dehydrogenase iron-sulfur subunit
            "K13924",  # dctA; C4-dicarboxylate transporter
        ],
        "essential_ko_groups": [
            ["K00239"],  # sdhA
            ["K13924"],  # dctA transporter
        ],
    },
    "malate": {
        "name": "Malate utilization",
        "kos": [
            "K00024",  # mdh; malate dehydrogenase
            "K00029",  # maeB; malate dehydrogenase (NADP+)
            "K01571",  # mleA; malolactic enzyme
        ],
        "essential_ko_groups": [
            ["K00024"],  # mdh
            ["K00029"],  # maeB (alt)
            ["K01571"],  # mleA (malolactic, LAB-specific)
        ],
    },
}

# Mineral transporter KOs
MINERAL_TRANSPORT_KOS = {
    "Mg": {
        "name": "Magnesium transport",
        "kos": [
            "K06213",  # mgtE; Mg2+ transporter
            "K01531",  # mgtA; Mg2+-importing ATPase
            "K06215",  # corA; Mg2+ transporter CorA
        ],
        "essential_ko_groups": [
            ["K06213"],  # mgtE
            ["K01531"],  # mgtA (alt)
            ["K06215"],  # corA (alt)
        ],
    },
    "Mn": {
        "name": "Manganese transport",
        "kos": [
            "K11707",  # mntH; Mn2+ transporter (NRAMP family)
            "K11709",  # mntA; Mn2+/Zn2+ ABC transporter substrate-binding
            "K11710",  # mntB; Mn2+/Zn2+ ABC transporter permease
            "K11604",  # sitA; Mn2+/Fe2+ ABC transporter
        ],
        "essential_ko_groups": [
            ["K11707"],  # mntH
            ["K11709"],  # mntA (alt)
            ["K11604"],  # sitA (alt, also Fe)
        ],
    },
    "Fe": {
        "name": "Iron transport",
        "kos": [
            "K02010",  # feoB; Fe2+ transporter FeoB
            "K02012",  # feoA; Fe2+ transporter FeoA
            "K01992",  # afuA; Fe3+ ABC transporter substrate-binding
            "K11604",  # sitA; Mn2+/Fe2+ ABC transporter
        ],
        "essential_ko_groups": [
            ["K02010"],  # feoB
            ["K01992"],  # afuA (alt)
            ["K11604"],  # sitA (alt)
        ],
    },
    "K": {
        "name": "Potassium transport",
        "kos": [
            "K03455",  # kdpA; K+-transporting ATPase subunit A
            "K03456",  # kdpB; K+-transporting ATPase subunit B
            "K03549",  # trkA; Trk K+ transport system
            "K06217",  # kup; K+ uptake protein (low affinity)
        ],
        "essential_ko_groups": [
            ["K03455"],  # kdpA
            ["K03549"],  # trkA (alt)
            ["K06217"],  # kup (alt)
        ],
    },
    "Ca": {
        "name": "Calcium transport",
        "kos": [
            "K01529",  # ATP2B; Ca2+-transporting ATPase
            "K07300",  # chaA; Ca2+/H+ antiporter
        ],
        "essential_ko_groups": [
            ["K01529"],  # ATP2B
            ["K07300"],  # chaA (alt)
        ],
    },
    "Na": {
        "name": "Sodium transport/homeostasis",
        "kos": [
            "K03313",  # nhaA; Na+/H+ antiporter
            "K03315",  # nhaB; Na+/H+ antiporter NhaB
        ],
        "essential_ko_groups": [
            ["K03313"],  # nhaA
            ["K03315"],  # nhaB (alt)
        ],
    },
}


def calculate_sugar_metabolism(ko_list: list[str]) -> dict[str, float]:
    """Calculate sugar metabolism capability for each sugar type.

    Args:
        ko_list: List of KO IDs found in the genome

    Returns:
        Dictionary of sugar_name -> utilization completeness (0.0 to 1.0)
    """
    ko_set = set(ko_list)
    result = {}
    for sugar, pathway_info in SUGAR_METABOLISM_KOS.items():
        completeness = calculate_pathway_completeness(
            ko_set, pathway_info, use_essential_only=True
        )
        result[sugar] = completeness
    return result


def calculate_organic_acid_metabolism(ko_list: list[str]) -> dict[str, float]:
    """Calculate organic acid metabolism capability.

    Args:
        ko_list: List of KO IDs found in the genome

    Returns:
        Dictionary of acid_name -> utilization completeness (0.0 to 1.0)
    """
    ko_set = set(ko_list)
    result = {}
    for acid, pathway_info in ORGANIC_ACID_METABOLISM_KOS.items():
        completeness = calculate_pathway_completeness(
            ko_set, pathway_info, use_essential_only=True
        )
        result[acid] = completeness
    return result


def calculate_mineral_transport(ko_list: list[str]) -> dict[str, float]:
    """Calculate mineral transport capability.

    Args:
        ko_list: List of KO IDs found in the genome

    Returns:
        Dictionary of mineral_name -> transport completeness (0.0 to 1.0)
    """
    ko_set = set(ko_list)
    result = {}
    for mineral, pathway_info in MINERAL_TRANSPORT_KOS.items():
        completeness = calculate_pathway_completeness(
            ko_set, pathway_info, use_essential_only=True
        )
        result[mineral] = completeness
    return result


def get_all_pathway_kos() -> set[str]:
    """Get all KO IDs used in pathway analysis.

    Returns:
        Set of all KO IDs
    """
    all_kos = set()

    for pathway_info in AA_BIOSYNTHESIS_KOS.values():
        all_kos.update(pathway_info["kos"])

    for pathway_info in VITAMIN_BIOSYNTHESIS_KOS.values():
        all_kos.update(pathway_info["kos"])

    for pathway_info in NUCLEOTIDE_BIOSYNTHESIS_KOS.values():
        all_kos.update(pathway_info["kos"])

    for transporter_info in TRANSPORTER_KOS.values():
        all_kos.update(transporter_info["kos"])

    for pathway_info in SUGAR_METABOLISM_KOS.values():
        all_kos.update(pathway_info["kos"])

    for pathway_info in ORGANIC_ACID_METABOLISM_KOS.values():
        all_kos.update(pathway_info["kos"])

    for pathway_info in MINERAL_TRANSPORT_KOS.values():
        all_kos.update(pathway_info["kos"])

    return all_kos
