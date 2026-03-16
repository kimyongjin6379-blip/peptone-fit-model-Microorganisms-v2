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

    Args:
        ko_list: List of KO IDs found in the genome
        pathway_kos: Dictionary with 'kos' and 'essential_kos' keys
        use_essential_only: If True, only check essential KOs

    Returns:
        Completeness score (0.0 to 1.0)
    """
    if use_essential_only and "essential_kos" in pathway_kos:
        target_kos = pathway_kos["essential_kos"]
    else:
        target_kos = pathway_kos["kos"]

    if not target_kos:
        return 1.0  # No KOs to check, assume complete

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
        "total_kos": len(set(ko_list)),
        "source": "gcf_annotation"
    }


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

    return all_kos
