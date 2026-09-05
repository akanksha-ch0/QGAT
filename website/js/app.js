/* ============================================================
   QGAT — Application Logic
   QTL and Gene Annotation Tool
   ============================================================ */

// ---- Mock Data ----
const MOCK_GENES = [
  {
    symbol: "DGAT1",
    id: "ENSBTAG00000026356",
    transcript_id: "ENSBTAT00000035128.4",
    chr: "Chr14",
    start: 1795425,
    end: 1804838,
    strand: "+",
    biotype: "protein_coding",
    overlap_bp: 9414,
    overlap_pct: 100.0,
    distance_to_tss: 0,
    qtl: "Milk Fat QTL",
    trait: "Milk fat percentage",
    score: 96,
    cds_length: 1467,
    protein_length: 489,
    exon_count: 17,
    key_variant: "p.Ala232Lys (K232A) in Exon 8 — Causal variant for milk fat content & yield in dairy cattle",
    pathways: ["Glycerolipid metabolism", "Fat digestion and absorption", "Triacylglycerol biosynthesis"],
    function: "Diacylglycerol O-acyltransferase 1, catalyzes the final step in triglyceride synthesis. Overlapping causal mutation on bovine chromosome 14.",
    exons: [
      { num: 1, start: 1795425, end: 1795640, len: 216, type: "5_prime_utr", phase: -1 },
      { num: 2, start: 1796100, end: 1796245, len: 146, type: "cds", phase: 0 },
      { num: 3, start: 1796850, end: 1797010, len: 161, type: "cds", phase: 2 },
      { num: 4, start: 1797520, end: 1797680, len: 161, type: "cds", phase: 1 },
      { num: 5, start: 1798300, end: 1798450, len: 151, type: "cds", phase: 0 },
      { num: 6, start: 1798900, end: 1799040, len: 141, type: "cds", phase: 2 },
      { num: 7, start: 1799500, end: 1799650, len: 151, type: "cds", phase: 1 },
      { num: 8, start: 1800200, end: 1800360, len: 161, type: "cds", phase: 0, is_causal: true, variant_note: "K232A Causal Mutation" },
      { num: 9, start: 1800900, end: 1801040, len: 141, type: "cds", phase: 2 },
      { num: 10, start: 1801550, end: 1801690, len: 141, type: "cds", phase: 1 },
      { num: 11, start: 1802200, end: 1802340, len: 141, type: "cds", phase: 0 },
      { num: 12, start: 1802800, end: 1802930, len: 131, type: "cds", phase: 2 },
      { num: 13, start: 1803350, end: 1803470, len: 121, type: "cds", phase: 1 },
      { num: 14, start: 1803850, end: 1803970, len: 121, type: "cds", phase: 0 },
      { num: 15, start: 1804150, end: 1804260, len: 111, type: "cds", phase: 2 },
      { num: 16, start: 1804420, end: 1804520, len: 101, type: "cds", phase: 1 },
      { num: 17, start: 1804680, end: 1804838, len: 159, type: "3_prime_utr", phase: -1 }
    ]
  },
  {
    symbol: "SCD",
    id: "ENSBTAG00000006757",
    transcript_id: "ENSBTAT00000008892.3",
    chr: "Chr26",
    start: 21144971,
    end: 21161782,
    strand: "+",
    biotype: "protein_coding",
    overlap_bp: 16812,
    overlap_pct: 100.0,
    distance_to_tss: 0,
    qtl: "Fat Composition QTL",
    trait: "Fatty acid composition",
    score: 92,
    cds_length: 1080,
    protein_length: 359,
    exon_count: 6,
    key_variant: "p.Ala293Val in Exon 5 — Significantly modifies monounsaturated fatty acid ratios in beef and milk",
    pathways: ["PPAR signaling pathway", "Fatty acid metabolism", "Biosynthesis of unsaturated fatty acids"],
    function: "Stearoyl-CoA desaturase, introduces cis-double bonds into stearoyl- and palmitoyl-CoA. Directly impacts beef marbling, tenderness, and healthy fatty acid profile.",
    exons: [
      { num: 1, start: 21144971, end: 21145220, len: 250, type: "5_prime_utr", phase: -1 },
      { num: 2, start: 21147800, end: 21148010, len: 211, type: "cds", phase: 0 },
      { num: 3, start: 21151200, end: 21151430, len: 231, type: "cds", phase: 1 },
      { num: 4, start: 21154500, end: 21154740, len: 241, type: "cds", phase: 2 },
      { num: 5, start: 21158100, end: 21158360, len: 261, type: "cds", phase: 0, is_causal: true, variant_note: "A293V Marbling Variant" },
      { num: 6, start: 21160900, end: 21161782, len: 883, type: "3_prime_utr", phase: -1 }
    ]
  },
  {
    symbol: "ABCG2",
    id: "ENSBTAG00000000578",
    transcript_id: "ENSBTAT00000000755.5",
    chr: "Chr6",
    start: 37956970,
    end: 38055839,
    strand: "+",
    biotype: "protein_coding",
    overlap_bp: 42100,
    overlap_pct: 42.6,
    distance_to_tss: 120,
    qtl: "Milk Yield QTL",
    trait: "Milk production",
    score: 89,
    cds_length: 1974,
    protein_length: 657,
    exon_count: 16,
    key_variant: "p.Tyr581Ser (Y581S) in Exon 14 — Direct effect on milk yield and somatic cell score in dairy cattle",
    pathways: ["ABC transporters", "Bile secretion", "Drug metabolism"],
    function: "ATP-binding cassette sub-family G member 2. Highly conserved xenobiotic and organic anion exporter in bovine mammary gland.",
    exons: [
      { num: 1, start: 37956970, end: 37957190, len: 221, type: "5_prime_utr", phase: -1 },
      { num: 2, start: 37965000, end: 37965150, len: 151, type: "cds", phase: 0 },
      { num: 3, start: 37972000, end: 37972160, len: 161, type: "cds", phase: 2 },
      { num: 4, start: 37980000, end: 37980140, len: 141, type: "cds", phase: 1 },
      { num: 5, start: 37988000, end: 37988150, len: 151, type: "cds", phase: 0 },
      { num: 6, start: 37996000, end: 37996160, len: 161, type: "cds", phase: 2 },
      { num: 7, start: 38004000, end: 38004140, len: 141, type: "cds", phase: 1 },
      { num: 8, start: 38012000, end: 38012150, len: 151, type: "cds", phase: 0 },
      { num: 9, start: 38020000, end: 38020150, len: 151, type: "cds", phase: 2 },
      { num: 10, start: 38028000, end: 38028140, len: 141, type: "cds", phase: 1 },
      { num: 11, start: 38035000, end: 38035150, len: 151, type: "cds", phase: 0 },
      { num: 12, start: 38042000, end: 38042140, len: 141, type: "cds", phase: 2 },
      { num: 13, start: 38047000, end: 38047150, len: 151, type: "cds", phase: 1 },
      { num: 14, start: 38051000, end: 38051180, len: 181, type: "cds", phase: 0, is_causal: true, variant_note: "Y581S Causal Variant" },
      { num: 15, start: 38053500, end: 38053640, len: 141, type: "cds", phase: 2 },
      { num: 16, start: 38055200, end: 38055839, len: 640, type: "3_prime_utr", phase: -1 }
    ]
  },
  {
    symbol: "MSTN",
    id: "ENSBTAG00000000563",
    transcript_id: "ENSBTAT00000000735.3",
    chr: "Chr2",
    start: 6444220,
    end: 6450160,
    strand: "+",
    biotype: "protein_coding",
    overlap_bp: 5941,
    overlap_pct: 100.0,
    distance_to_tss: 0,
    qtl: "Muscle Mass QTL",
    trait: "Muscle development",
    score: 81,
    cds_length: 1128,
    protein_length: 375,
    exon_count: 3,
    key_variant: "nt821(del11) in Exon 3 — 11-bp deletion causing premature stop codon & double-muscling phenotype in Belgian Blue / Piedmontese",
    pathways: ["TGF-beta signaling pathway", "MAPK signaling", "Skeletal muscle hypertrophy"],
    function: "Myostatin (GDF-8), essential negative regulator of skeletal muscle growth. Loss-of-function variants dramatically expand muscle fiber number and carcass yield.",
    exons: [
      { num: 1, start: 6444220, end: 6444602, len: 383, type: "5_prime_utr", phase: -1 },
      { num: 2, start: 6446450, end: 6446824, len: 375, type: "cds", phase: 0 },
      { num: 3, start: 6449100, end: 6450160, len: 1061, type: "cds", phase: 1, is_causal: true, variant_note: "nt821(del11) Double Muscling" }
    ]
  },
  {
    symbol: "LEP",
    id: "ENSBTAG00000005765",
    transcript_id: "ENSBTAT00000007584.4",
    chr: "Chr4",
    start: 93266599,
    end: 93283605,
    strand: "-",
    biotype: "protein_coding",
    overlap_bp: 17007,
    overlap_pct: 100.0,
    distance_to_tss: 0,
    qtl: "Body Weight QTL",
    trait: "Growth rate",
    score: 87,
    cds_length: 504,
    protein_length: 167,
    exon_count: 3,
    key_variant: "c.73C>T (Arg25Cys) in Exon 2 — Associated with subcutaneous backfat depth and voluntary feed intake in beef steers",
    pathways: ["Adipocytokine signaling", "Neuroactive ligand-receptor", "AMPK signaling pathway"],
    function: "Leptin, adipocyte-secreted hormone signaling energy reserves to the hypothalamus. Regulates satiety, metabolic rate, and adipose deposition.",
    exons: [
      { num: 1, start: 93266599, end: 93266780, len: 182, type: "5_prime_utr", phase: -1 },
      { num: 2, start: 93275200, end: 93275372, len: 173, type: "cds", phase: 0, is_causal: true, variant_note: "Arg25Cys Feed Intake Variant" },
      { num: 3, start: 93282800, end: 93283605, len: 806, type: "3_prime_utr", phase: 1 }
    ]
  },
  {
    symbol: "GHR",
    id: "ENSBTAG00000001335",
    transcript_id: "ENSBTAT00000001740.4",
    chr: "Chr20",
    start: 31890736,
    end: 32162460,
    strand: "+",
    biotype: "protein_coding",
    overlap_bp: 271725,
    overlap_pct: 90.5,
    distance_to_tss: 0,
    qtl: "Growth QTL",
    trait: "Body size",
    score: 85,
    cds_length: 1917,
    protein_length: 638,
    exon_count: 9,
    key_variant: "p.Phe279Tyr in Exon 8 — Modulates IGF-1 release and adult skeletal stature across cattle breeds",
    pathways: ["JAK-STAT signaling pathway", "Cytokine-cytokine receptor", "PI3K-Akt signaling"],
    function: "Growth hormone receptor. Transduces somatotropin signaling to drive hepatic IGF1 production, regulating post-weaning gain and frame size.",
    exons: [
      { num: 1, start: 31890736, end: 31890950, len: 215, type: "5_prime_utr", phase: -1 },
      { num: 2, start: 31920000, end: 31920180, len: 181, type: "cds", phase: 0 },
      { num: 3, start: 31950000, end: 31950190, len: 191, type: "cds", phase: 1 },
      { num: 4, start: 31985000, end: 31985180, len: 181, type: "cds", phase: 0 },
      { num: 5, start: 32025000, end: 32025170, len: 171, type: "cds", phase: 2 },
      { num: 6, start: 32065000, end: 32065210, len: 211, type: "cds", phase: 1 },
      { num: 7, start: 32105000, end: 32105190, len: 191, type: "cds", phase: 0 },
      { num: 8, start: 32135000, end: 32135220, len: 221, type: "cds", phase: 2, is_causal: true, variant_note: "F279Y Stature Variant" },
      { num: 9, start: 32161200, end: 32162460, len: 1261, type: "3_prime_utr", phase: -1 }
    ]
  },
  {
    symbol: "CAST",
    id: "ENSBTAG00000010312",
    transcript_id: "ENSBTAT00000013589.4",
    chr: "Chr7",
    start: 98350000,
    end: 98470000,
    strand: "+",
    biotype: "protein_coding",
    overlap_bp: 120000,
    overlap_pct: 100.0,
    distance_to_tss: 0,
    qtl: "Meat Tenderness QTL",
    trait: "Meat quality",
    score: 83,
    cds_length: 2145,
    protein_length: 714,
    exon_count: 8,
    key_variant: "CAST:c.2870A>G in 3' UTR — Primary genetic marker for Warner-Bratzler shear force and beef tenderness",
    pathways: ["Calpain-calpastatin system", "Post-mortem tenderization", "Apoptosis"],
    function: "Calpastatin, endogenous specific inhibitor of calpains (calcium-activated neutral proteases). Central regulator of muscle protein turnover and beef tenderness.",
    exons: [
      { num: 1, start: 98350000, end: 98350250, len: 251, type: "5_prime_utr", phase: -1 },
      { num: 2, start: 98370000, end: 98370280, len: 281, type: "cds", phase: 0 },
      { num: 3, start: 98390000, end: 98390310, len: 311, type: "cds", phase: 2 },
      { num: 4, start: 98410000, end: 98410300, len: 301, type: "cds", phase: 1 },
      { num: 5, start: 98430000, end: 98430290, len: 291, type: "cds", phase: 0 },
      { num: 6, start: 98445000, end: 98445320, len: 321, type: "cds", phase: 2 },
      { num: 7, start: 98458000, end: 98458310, len: 311, type: "cds", phase: 1 },
      { num: 8, start: 98468500, end: 98470000, len: 1501, type: "3_prime_utr", phase: -1, is_causal: true, variant_note: "Tenderness Marker Marker" }
    ]
  },
  {
    symbol: "PLAG1",
    id: "ENSBTAG00000019711",
    transcript_id: "ENSBTAT00000025981.3",
    chr: "Chr14",
    start: 25006221,
    end: 25059820,
    strand: "+",
    biotype: "protein_coding",
    overlap_bp: 53600,
    overlap_pct: 100.0,
    distance_to_tss: 0,
    qtl: "Stature QTL",
    trait: "Body height",
    score: 79,
    cds_length: 1503,
    protein_length: 500,
    exon_count: 5,
    key_variant: "rs109231213 upstream regulatory variant — Major determinant of skeletal height & birth weight across Bos taurus and Bos indicus",
    pathways: ["Transcriptional regulation", "Wnt signaling", "Skeletal development"],
    function: "Pleiomorphic adenoma gene 1, zinc-finger transcription factor governing fetal growth, withers height, and weight gain across ruminants.",
    exons: [
      { num: 1, start: 25006221, end: 25006450, len: 230, type: "5_prime_utr", phase: -1 },
      { num: 2, start: 25018000, end: 25018320, len: 321, type: "cds", phase: 0 },
      { num: 3, start: 25032000, end: 25032410, len: 411, type: "cds", phase: 1 },
      { num: 4, start: 25045000, end: 25045420, len: 421, type: "cds", phase: 2 },
      { num: 5, start: 25058500, end: 25059820, len: 1321, type: "3_prime_utr", phase: 0 }
    ]
  },
  {
    symbol: "CSN1S1",
    id: "ENSBTAG00000003063",
    transcript_id: "ENSBTAT00000004018.4",
    chr: "Chr6",
    start: 87141556,
    end: 87159808,
    strand: "+",
    biotype: "protein_coding",
    overlap_bp: 18253,
    overlap_pct: 100.0,
    distance_to_tss: 0,
    qtl: "Milk Protein QTL",
    trait: "Milk protein content",
    score: 77,
    cds_length: 645,
    protein_length: 214,
    exon_count: 19,
    key_variant: "CSN1S1*B vs *C allele — Direct impact on cheese yield, curd firmness, and rennet coagulation time",
    pathways: ["Mammary gland lactation", "Calcium signaling", "Protein processing in ER"],
    function: "Alpha-S1-casein, dominant phosphoprotein of livestock milk. Key nutritional determinant of micelle structure and cheese-making traits.",
    exons: [
      { num: 1, start: 87141556, end: 87141680, len: 125, type: "5_prime_utr", phase: -1 },
      { num: 2, start: 87143000, end: 87143120, len: 121, type: "cds", phase: 0 },
      { num: 3, start: 87145000, end: 87145100, len: 101, type: "cds", phase: 1 },
      { num: 4, start: 87147500, end: 87147580, len: 81, type: "cds", phase: 2 },
      { num: 5, start: 87150000, end: 87150090, len: 91, type: "cds", phase: 0 },
      { num: 6, start: 87153000, end: 87153100, len: 101, type: "cds", phase: 1 },
      { num: 7, start: 87156000, end: 87156100, len: 101, type: "cds", phase: 2 },
      { num: 8, start: 87158500, end: 87159808, len: 1309, type: "3_prime_utr", phase: -1 }
    ]
  },
  {
    symbol: "SLC11A1",
    id: "ENSBTAG00000001158",
    transcript_id: "ENSBTAT00000001512.4",
    chr: "Chr2",
    start: 113500000,
    end: 113516000,
    strand: "-",
    biotype: "protein_coding",
    overlap_bp: 16000,
    overlap_pct: 100.0,
    distance_to_tss: 0,
    qtl: "Disease Resistance QTL",
    trait: "Disease resistance",
    score: 74,
    cds_length: 1647,
    protein_length: 548,
    exon_count: 15,
    key_variant: "3' UTR GT-microsatellite polymorphism — Conveys innate macrophage resistance against Brucella abortus and tuberculosis",
    pathways: ["Phagosome maturation", "Mineral absorption", "Innate antimicrobial immunity"],
    function: "Solute carrier family 11 member 1 (NRAMP1). Divalent transition metal transporter regulating intracellular bacterial pathogen survival in macrophages.",
    exons: [
      { num: 1, start: 113500000, end: 113500220, len: 221, type: "5_prime_utr", phase: -1 },
      { num: 2, start: 113502500, end: 113502680, len: 181, type: "cds", phase: 0 },
      { num: 3, start: 113505000, end: 113505160, len: 161, type: "cds", phase: 2 },
      { num: 4, start: 113507500, end: 113507670, len: 171, type: "cds", phase: 1 },
      { num: 5, start: 113510000, end: 113510160, len: 161, type: "cds", phase: 0 },
      { num: 6, start: 113512500, end: 113512680, len: 181, type: "cds", phase: 2 },
      { num: 7, start: 113515000, end: 113516000, len: 1001, type: "3_prime_utr", phase: -1, is_causal: true, variant_note: "Resistance Microsatellite" }
    ]
  },
  {
    symbol: "CAPN1",
    id: "ENSBTAG00000009768",
    transcript_id: "ENSBTAT00000012879.4",
    chr: "Chr29",
    start: 44068000,
    end: 44093000,
    strand: "+",
    biotype: "protein_coding",
    overlap_bp: 25000,
    overlap_pct: 100.0,
    distance_to_tss: 0,
    qtl: "Meat Quality QTL",
    trait: "Meat tenderness",
    score: 68,
    cds_length: 2148,
    protein_length: 715,
    exon_count: 22,
    key_variant: "CAPN1_316 (p.Ala316Gly) and CAPN1_4751 markers — Major commercial DNA markers for tender beef assurance",
    pathways: ["Calcium signaling", "Calpain proteolytic cascade", "Skeletal muscle proteolysis"],
    function: "Calpain 1 (mu-calpain), calcium-activated intracellular cysteine protease. Primary enzyme responsible for post-mortem myofibrillar degradation and meat tenderization.",
    exons: [
      { num: 1, start: 44068000, end: 44068180, len: 181, type: "5_prime_utr", phase: -1 },
      { num: 2, start: 44071500, end: 44071680, len: 181, type: "cds", phase: 0 },
      { num: 3, start: 44075000, end: 44075190, len: 191, type: "cds", phase: 2 },
      { num: 4, start: 44079000, end: 44079210, len: 211, type: "cds", phase: 1 },
      { num: 5, start: 44083000, end: 44083200, len: 201, type: "cds", phase: 0 },
      { num: 6, start: 44087000, end: 44087220, len: 221, type: "cds", phase: 2, is_causal: true, variant_note: "CAPN1_316 Marker" },
      { num: 7, start: 44091500, end: 44093000, len: 1501, type: "3_prime_utr", phase: -1 }
    ]
  }
];

const MOCK_QTLS = [
  { id: "QTL:12450", species: "cattle", speciesName: "Cattle (Bos taurus)", chr: "Chr14", start: 1750000, end: 1850000, trait: "Milk fat percentage", category: "Milk Production", pubmed: "PMID:11752453", overlap_pct: 98.4, candidate: "DGAT1", desc: "Major QTL affecting milk fat percentage, fat yield, and triglyceride synthesis in dairy cattle breeds worldwide." },
  { id: "QTL:38902", species: "cattle", speciesName: "Cattle (Bos taurus)", chr: "Chr6", start: 37900000, end: 38200000, trait: "Milk yield & protein content", category: "Milk Production", pubmed: "PMID:15654098", overlap_pct: 94.2, candidate: "ABCG2", desc: "Significant QTL associated with lactation volume, somatic cell score, and milk protein yield." },
  { id: "QTL:08941", species: "cattle", speciesName: "Cattle (Bos taurus)", chr: "Chr2", start: 6400000, end: 6500000, trait: "Muscle hypertrophy & double muscling", category: "Meat / Muscle", pubmed: "PMID:9367989", overlap_pct: 100.0, candidate: "MSTN", desc: "Regulates skeletal muscle growth and carcass weight. Inactivation results in double-muscling phenotype." },
  { id: "QTL:45210", species: "cattle", speciesName: "Cattle (Bos taurus)", chr: "Chr20", start: 31800000, end: 32200000, trait: "Growth rate & mature body size", category: "Growth & Size", pubmed: "PMID:17991950", overlap_pct: 88.5, candidate: "GHR", desc: "Growth hormone receptor locus associated with post-weaning average daily gain and frame size." },
  { id: "QTL:19042", species: "cattle", speciesName: "Cattle (Bos taurus)", chr: "Chr14", start: 24900000, end: 25150000, trait: "Stature & skeletal height", category: "Conformation", pubmed: "PMID:24608383", overlap_pct: 91.0, candidate: "PLAG1", desc: "Associated with hip height, body length, and early skeletal development across Bos taurus and Bos indicus." },
  { id: "QTL:67219", species: "cattle", speciesName: "Cattle (Bos taurus)", chr: "Chr26", start: 21100000, end: 21200000, trait: "Fatty acid composition & marbling", category: "Meat Quality", pubmed: "PMID:16327092", overlap_pct: 86.7, candidate: "SCD", desc: "Affects monounsaturated fatty acid ratios and beef tenderness." },
  { id: "QTL:89104", species: "goat", speciesName: "Goat (Capra hircus)", chr: "Chr6", start: 86000000, end: 87500000, trait: "Milk casein synthesis & coagulation", category: "Milk Production", pubmed: "PMID:11244199", overlap_pct: 95.0, candidate: "CSN1S1", desc: "Key casein locus in goats determining curd firming and cheese production efficiency." },
  { id: "QTL:78214", species: "sheep", speciesName: "Sheep (Ovis aries)", chr: "ChrX", start: 48000000, end: 49500000, trait: "Ovulation rate & litter size", category: "Reproduction", pubmed: "PMID:11912440", overlap_pct: 92.3, candidate: "BMP15", desc: "FecX locus modulating ovulation rate and twin/triplet lambing frequency in prolific sheep." },
  { id: "QTL:54120", species: "pig", speciesName: "Pig (Sus scrofa)", chr: "Chr2", start: 1400000, end: 1600000, trait: "Muscle mass & lean meat percentage", category: "Meat / Muscle", pubmed: "PMID:14574411", overlap_pct: 97.1, candidate: "IGF2", desc: "Major regulatory mutation in intron 3 of IGF2 influencing muscle accretion and fat deposition." },
  { id: "QTL:33190", species: "chicken", speciesName: "Chicken (Gallus gallus)", chr: "Chr4", start: 72000000, end: 74000000, trait: "Body weight & feed conversion ratio", category: "Growth & Size", pubmed: "PMID:15582963", overlap_pct: 89.0, candidate: "CCKAR", desc: "QTL governing feed intake regulation and rapid growth in broiler chickens." },
  { id: "QTL:91240", species: "horse", speciesName: "Horse (Equus caballus)", chr: "Chr18", start: 66000000, end: 67000000, trait: "Sprint racing performance & speed", category: "Performance", pubmed: "PMID:20147631", overlap_pct: 96.5, candidate: "MSTN", desc: "The equine speed gene determining optimal racing distance (sprint vs staying) in Thoroughbred horses." }
];

const MOCK_QTL_REGIONS = [
  { id: "QTL:12450", chr: "Chr14", start: 1500000, end: 3000000, trait: "Milk fat percentage", color: "rgba(139,92,246,0.35)" },
  { id: "QTL:38902", chr: "Chr6", start: 37000000, end: 42000000, trait: "Milk production", color: "rgba(139,92,246,0.25)" },
  { id: "QTL:19042", chr: "Chr14", start: 24000000, end: 26000000, trait: "Body height", color: "rgba(139,92,246,0.3)" },
  { id: "QTL:08941", chr: "Chr2", start: 6000000, end: 7000000, trait: "Muscle development", color: "rgba(139,92,246,0.25)" },
  { id: "QTL:45210", chr: "Chr20", start: 31500000, end: 32500000, trait: "Growth rate", color: "rgba(139,92,246,0.3)" },
];

const TRAIT_CATEGORIES = [
  { name: "Milk Production", icon: "🥛", count: 4521, color: "#0d9488" },
  { name: "Growth & Body Size", icon: "📈", count: 3287, color: "#8b5cf6" },
  { name: "Reproduction", icon: "🧬", count: 2456, color: "#f43f5e" },
  { name: "Disease Resistance", icon: "🛡️", count: 1893, color: "#10b981" },
  { name: "Milk Composition", icon: "🧈", count: 1654, color: "#f59e0b" },
  { name: "Meat Quality", icon: "🥩", count: 1432, color: "#06b6d4" },
  { name: "Heat Tolerance", icon: "🌡️", count: 876, color: "#ef4444" },
  { name: "Feed Efficiency", icon: "🌾", count: 743, color: "#84cc16" },
  { name: "Coat & Morphology", icon: "🎨", count: 654, color: "#a855f7" },
  { name: "Carcass Traits", icon: "📐", count: 1123, color: "#ec4899" },
  { name: "Immune Response", icon: "🔬", count: 567, color: "#14b8a6" },
  { name: "Fertility & Calving", icon: "🐄", count: 1245, color: "#6366f1" },
];

const ASSEMBLIES = {
  human: ["GRCh38 (hg38)", "GRCh37 (hg19)"],
  mouse: ["GRCm39 (mm39)", "GRCm38 (mm10)"],
  cattle: ["ARS-UCD1.2", "ARS-UCD1.3", "UMD3.1"],
  buffalo: ["UOA_WB_1"],
  goat: ["ARS1 (CHIR_1.0)", "ARS1.2"],
  sheep: ["Oar_v3.1", "Oar_rambouillet_v1.0"],
  pig: ["Sscrofa11.1", "Sscrofa10.2"],
  chicken: ["GRCg6a", "GRCg7b"],
  horse: ["EquCab3.0", "EquCab2.0"],
};

// ---- App State ----
const state = {
  currentPage: "home",
  annotationType: "gene",
  selectedAnimal: "cattle",
  activeResultsTab: "genes",
  geneFormat: "range",
  qtlFormat: "range",
  combinedFormat: "range",
  resultsLoaded: false,
  genomeZoom: 1,
  genomePan: 0,
  hoveredGene: null,
  ontologyGenes: [],
  ontologySpecies: "btaurus"
};

// ============================================================
//   ROUTING
// ============================================================
function navigateTo(page) {
  state.currentPage = page;

  // Hide all sections
  document.querySelectorAll(".page-section").forEach(s => {
    s.classList.remove("active");
  });

  // Show target
  const target = document.getElementById(`section-${page}`);
  if (target) {
    // Small delay so CSS transition works
    requestAnimationFrame(() => {
      target.classList.add("active");
    });
  }

  // Update nav links
  document.querySelectorAll(".nav__link").forEach(l => {
    l.classList.toggle("active", l.dataset.nav === page);
  });

  // Update nav style
  const nav = document.getElementById("mainNav");
  if (page === "home") {
    nav.className = "nav nav--transparent";
    window.scrollTo(0, 0);
  } else {
    nav.className = "nav nav--solid";
    window.scrollTo(0, 0);
  }

  // Trigger results animations if navigating to results
  if (page === "results" && !state.resultsLoaded) {
    state.resultsLoaded = true;
    setTimeout(animateResults, 300);
  }

  if (page === "results") {
    setTimeout(() => drawGenomeBrowser(), 100);
  }

  if (page === "ontology") {
    updateOntologyView();
  }
}

// ============================================================
//   HERO CANVAS — Animated Genomic Visualization
// ============================================================
class HeroAnimation {
  constructor(canvas) {
    this.canvas = canvas;
    this.ctx = canvas.getContext("2d");
    this.particles = [];
    this.helixPoints = [];
    this.geneBlocks = [];
    this.qtlRegions = [];
    this.tracks = [];
    this.frame = 0;
    this.resize();
    this.init();
    window.addEventListener("resize", () => this.resize());
  }

  resize() {
    const dpr = window.devicePixelRatio || 1;
    this.canvas.width = this.canvas.offsetWidth * dpr;
    this.canvas.height = this.canvas.offsetHeight * dpr;
    this.ctx.scale(dpr, dpr);
    this.w = this.canvas.offsetWidth;
    this.h = this.canvas.offsetHeight;
  }

  init() {
    // Particles
    for (let i = 0; i < 80; i++) {
      this.particles.push({
        x: Math.random() * this.w,
        y: Math.random() * this.h,
        size: Math.random() * 2 + 0.5,
        speedX: (Math.random() - 0.5) * 0.3,
        speedY: (Math.random() - 0.5) * 0.3,
        opacity: Math.random() * 0.4 + 0.1,
      });
    }

    // Gene blocks along bottom tracks
    const trackY = this.h * 0.68;
    for (let i = 0; i < 12; i++) {
      this.geneBlocks.push({
        x: 80 + i * (this.w - 160) / 12 + Math.random() * 30,
        y: trackY + (i % 3) * 25,
        width: 30 + Math.random() * 60,
        height: 8,
        color: i % 3 === 0 ? "rgba(13,148,136," : i % 3 === 1 ? "rgba(6,182,212," : "rgba(16,185,129,",
        opacity: 0.3 + Math.random() * 0.4,
      });
    }

    // QTL regions
    for (let i = 0; i < 4; i++) {
      this.qtlRegions.push({
        x: 120 + i * (this.w - 240) / 4,
        width: 80 + Math.random() * 120,
        opacity: 0.08 + Math.random() * 0.08,
      });
    }
  }

  drawHelix(time) {
    const ctx = this.ctx;
    const centerX = this.w * 0.75;
    const startY = this.h * 0.15;
    const helixH = this.h * 0.55;
    const amplitude = 40;
    const frequency = 0.025;

    for (let strand = 0; strand < 2; strand++) {
      const phase = strand * Math.PI;
      ctx.beginPath();
      ctx.strokeStyle = strand === 0 ? "rgba(13,148,136,0.2)" : "rgba(139,92,246,0.2)";
      ctx.lineWidth = 2;

      for (let y = 0; y < helixH; y += 2) {
        const x = centerX + Math.sin(y * frequency + time * 0.001 + phase) * amplitude;
        const absY = startY + y;
        if (y === 0) ctx.moveTo(x, absY);
        else ctx.lineTo(x, absY);
      }
      ctx.stroke();
    }

    // Rungs
    for (let y = 0; y < helixH; y += 20) {
      const x1 = centerX + Math.sin(y * frequency + time * 0.001) * amplitude;
      const x2 = centerX + Math.sin(y * frequency + time * 0.001 + Math.PI) * amplitude;
      const absY = startY + y;
      const gradient = ctx.createLinearGradient(x1, absY, x2, absY);
      gradient.addColorStop(0, "rgba(13,148,136,0.15)");
      gradient.addColorStop(0.5, "rgba(6,182,212,0.08)");
      gradient.addColorStop(1, "rgba(139,92,246,0.15)");
      ctx.beginPath();
      ctx.strokeStyle = gradient;
      ctx.lineWidth = 1;
      ctx.moveTo(x1, absY);
      ctx.lineTo(x2, absY);
      ctx.stroke();
    }
  }

  drawTracks(time) {
    const ctx = this.ctx;
    const trackY = this.h * 0.65;

    // Chromosome line
    ctx.beginPath();
    ctx.strokeStyle = "rgba(255,255,255,0.08)";
    ctx.lineWidth = 2;
    ctx.moveTo(60, trackY);
    ctx.lineTo(this.w - 60, trackY);
    ctx.stroke();

    // Coordinate ticks
    for (let i = 0; i < 20; i++) {
      const x = 60 + i * (this.w - 120) / 20;
      ctx.beginPath();
      ctx.strokeStyle = "rgba(255,255,255,0.05)";
      ctx.moveTo(x, trackY - 5);
      ctx.lineTo(x, trackY + 5);
      ctx.stroke();
    }

    // QTL highlight regions
    this.qtlRegions.forEach(q => {
      const grad = ctx.createLinearGradient(q.x, trackY - 50, q.x, trackY + 60);
      grad.addColorStop(0, `rgba(139,92,246,0)`);
      grad.addColorStop(0.3, `rgba(139,92,246,${q.opacity})`);
      grad.addColorStop(0.7, `rgba(139,92,246,${q.opacity})`);
      grad.addColorStop(1, `rgba(139,92,246,0)`);
      ctx.fillStyle = grad;
      ctx.fillRect(q.x, trackY - 50, q.width, 110);
    });

    // Gene blocks
    this.geneBlocks.forEach((g, i) => {
      const pulse = Math.sin(time * 0.002 + i) * 0.1;
      ctx.fillStyle = g.color + (g.opacity + pulse) + ")";
      ctx.beginPath();
      ctx.roundRect(g.x, g.y, g.width, g.height, 3);
      ctx.fill();
    });

    // Marker dots
    for (let i = 0; i < 25; i++) {
      const x = 80 + i * (this.w - 160) / 25;
      const y = trackY + 55;
      const pulse = Math.sin(time * 0.003 + i * 0.5) * 0.15;
      ctx.fillStyle = `rgba(245,158,11,${0.25 + pulse})`;
      ctx.beginPath();
      ctx.arc(x, y, 2.5, 0, Math.PI * 2);
      ctx.fill();
    }

    // Labels
    ctx.fillStyle = "rgba(255,255,255,0.15)";
    ctx.font = "11px 'JetBrains Mono', monospace";
    ctx.fillText("Chromosome 14", 60, trackY - 15);
    ctx.fillText("Genes", 60, this.h * 0.68 - 8);
    ctx.fillText("QTL", 60, trackY - 40);
    ctx.fillText("Markers", 60, trackY + 62);
  }

  drawParticles(time) {
    const ctx = this.ctx;
    this.particles.forEach(p => {
      p.x += p.speedX;
      p.y += p.speedY;
      if (p.x < 0) p.x = this.w;
      if (p.x > this.w) p.x = 0;
      if (p.y < 0) p.y = this.h;
      if (p.y > this.h) p.y = 0;

      const pulse = Math.sin(time * 0.002 + p.x * 0.01) * 0.1;
      ctx.fillStyle = `rgba(13,148,136,${p.opacity + pulse})`;
      ctx.beginPath();
      ctx.arc(p.x, p.y, p.size, 0, Math.PI * 2);
      ctx.fill();
    });

    // Draw connections between nearby particles
    for (let i = 0; i < this.particles.length; i++) {
      for (let j = i + 1; j < this.particles.length; j++) {
        const dx = this.particles[i].x - this.particles[j].x;
        const dy = this.particles[i].y - this.particles[j].y;
        const dist = Math.sqrt(dx * dx + dy * dy);
        if (dist < 100) {
          ctx.beginPath();
          ctx.strokeStyle = `rgba(13,148,136,${0.04 * (1 - dist / 100)})`;
          ctx.lineWidth = 0.5;
          ctx.moveTo(this.particles[i].x, this.particles[i].y);
          ctx.lineTo(this.particles[j].x, this.particles[j].y);
          ctx.stroke();
        }
      }
    }
  }

  animate(time) {
    this.ctx.clearRect(0, 0, this.w, this.h);
    this.drawParticles(time);
    this.drawHelix(time);
    this.drawTracks(time);
    this.frame = requestAnimationFrame(t => this.animate(t));
  }

  start() {
    this.frame = requestAnimationFrame(t => this.animate(t));
  }

  stop() {
    cancelAnimationFrame(this.frame);
  }
}

// ============================================================
//   GENOME BROWSER — Interactive Canvas
// ============================================================
function drawGenomeBrowser() {
  const canvas = document.getElementById("genomeCanvas");
  if (!canvas) return;
  const ctx = canvas.getContext("2d");
  const dpr = window.devicePixelRatio || 1;
  canvas.width = canvas.offsetWidth * dpr;
  canvas.height = canvas.offsetHeight * dpr;
  ctx.scale(dpr, dpr);
  const W = canvas.offsetWidth;
  const H = canvas.offsetHeight;

  // Viewport
  const zoom = state.genomeZoom;
  const pan = state.genomePan;
  const viewStart = 0 + pan;
  const viewEnd = 50000000 / zoom + pan;
  const bpToX = (bp) => ((bp - viewStart) / (viewEnd - viewStart)) * W;

  ctx.clearRect(0, 0, W, H);

  // Background
  ctx.fillStyle = "#0a1628";
  ctx.fillRect(0, 0, W, H);

  // Grid lines
  const gridInterval = Math.pow(10, Math.floor(Math.log10((viewEnd - viewStart) / 10)));
  ctx.strokeStyle = "rgba(255,255,255,0.04)";
  ctx.lineWidth = 1;
  for (let bp = Math.ceil(viewStart / gridInterval) * gridInterval; bp < viewEnd; bp += gridInterval) {
    const x = bpToX(bp);
    ctx.beginPath();
    ctx.moveTo(x, 0);
    ctx.lineTo(x, H);
    ctx.stroke();
  }

  // Coordinate axis
  const axisY = 30;
  ctx.strokeStyle = "rgba(255,255,255,0.15)";
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(0, axisY);
  ctx.lineTo(W, axisY);
  ctx.stroke();

  // Coordinate labels
  ctx.fillStyle = "rgba(255,255,255,0.4)";
  ctx.font = "10px 'JetBrains Mono', monospace";
  const labelInterval = gridInterval;
  for (let bp = Math.ceil(viewStart / labelInterval) * labelInterval; bp < viewEnd; bp += labelInterval) {
    const x = bpToX(bp);
    ctx.fillText(formatBp(bp), x + 3, axisY - 6);
    ctx.beginPath();
    ctx.strokeStyle = "rgba(255,255,255,0.15)";
    ctx.moveTo(x, axisY - 3);
    ctx.lineTo(x, axisY + 3);
    ctx.stroke();
  }

  // Chromosome label
  ctx.fillStyle = "rgba(255,255,255,0.6)";
  ctx.font = "bold 12px 'Plus Jakarta Sans', sans-serif";
  ctx.fillText("Chromosome 14", 10, 18);

  // -- QTL track --
  const qtlTrackY = 50;
  ctx.fillStyle = "rgba(255,255,255,0.25)";
  ctx.font = "10px 'Plus Jakarta Sans', sans-serif";
  ctx.fillText("QTL Regions", 10, qtlTrackY + 4);

  const chr14Qtls = MOCK_QTL_REGIONS.filter(q => q.chr === "Chr14");
  chr14Qtls.forEach(q => {
    const x1 = bpToX(q.start);
    const x2 = bpToX(q.end);
    const w = Math.max(x2 - x1, 4);
    // Background glow
    const grad = ctx.createLinearGradient(x1, qtlTrackY + 10, x1, qtlTrackY + 30);
    grad.addColorStop(0, "rgba(139,92,246,0.15)");
    grad.addColorStop(1, "rgba(139,92,246,0.05)");
    ctx.fillStyle = grad;
    ctx.fillRect(x1, qtlTrackY + 10, w, 20);
    // Border
    ctx.fillStyle = "rgba(139,92,246,0.5)";
    ctx.fillRect(x1, qtlTrackY + 10, w, 2);
    ctx.fillRect(x1, qtlTrackY + 28, w, 2);
    // Label
    if (w > 40) {
      ctx.fillStyle = "rgba(139,92,246,0.7)";
      ctx.font = "9px 'JetBrains Mono', monospace";
      ctx.fillText(q.id, x1 + 4, qtlTrackY + 23);
    }
  });

  // -- Gene track --
  const geneTrackY = 90;
  ctx.fillStyle = "rgba(255,255,255,0.25)";
  ctx.font = "10px 'Plus Jakarta Sans', sans-serif";
  ctx.fillText("Genes", 10, geneTrackY + 4);

  // Create some genes on Chr14
  const chr14Genes = [
    { symbol: "DGAT1", start: 1795425, end: 1804838, candidate: true },
    { symbol: "PLAG1", start: 25006221, end: 25059820, candidate: true },
    { symbol: "FOXH1", start: 1850000, end: 1870000, candidate: false },
    { symbol: "CPSF1", start: 2100000, end: 2180000, candidate: false },
    { symbol: "CYHR1", start: 2350000, end: 2380000, candidate: false },
    { symbol: "ZC3H3", start: 3500000, end: 3650000, candidate: false },
    { symbol: "LY6K", start: 5200000, end: 5250000, candidate: false },
    { symbol: "PENK", start: 8600000, end: 8620000, candidate: false },
    { symbol: "TG", start: 9400000, end: 9690000, candidate: false },
    { symbol: "MOS", start: 12000000, end: 12010000, candidate: false },
    { symbol: "RGS20", start: 15800000, end: 16050000, candidate: false },
    { symbol: "NCOA2", start: 19000000, end: 19350000, candidate: false },
    { symbol: "KHDRBS3", start: 24200000, end: 24500000, candidate: false },
    { symbol: "INTS1", start: 26500000, end: 26800000, candidate: false },
    { symbol: "DEPTOR", start: 32000000, end: 32200000, candidate: false },
    { symbol: "XKR4", start: 400000, end: 500000, candidate: false },
  ];

  const geneItems = [];
  chr14Genes.forEach((g, i) => {
    const x1 = bpToX(g.start);
    const x2 = bpToX(g.end);
    const w = Math.max(x2 - x1, 6);
    const row = i % 4;
    const y = geneTrackY + 12 + row * 26;

    // Store for hover
    geneItems.push({ ...g, pixelX: x1, pixelY: y, pixelW: w, pixelH: 14 });

    // Gene body
    const color = g.candidate ? "rgba(244,63,94," : "rgba(13,148,136,";
    ctx.fillStyle = color + "0.6)";
    ctx.beginPath();
    ctx.roundRect(x1, y, w, 14, 3);
    ctx.fill();

    // Gene label
    if (w > 25) {
      ctx.fillStyle = g.candidate ? "rgba(244,63,94,0.9)" : "rgba(13,148,136,0.9)";
      ctx.font = "bold 9px 'JetBrains Mono', monospace";
      ctx.fillText(g.symbol, x1 + 3, y + 11);
    }

    // Direction arrow
    ctx.fillStyle = color + "0.4)";
    ctx.beginPath();
    ctx.moveTo(x1 + w, y + 7);
    ctx.lineTo(x1 + w + 5, y + 7);
    ctx.lineTo(x1 + w, y + 2);
    ctx.fill();
  });

  // -- Marker track --
  const markerY = geneTrackY + 125;
  ctx.fillStyle = "rgba(255,255,255,0.25)";
  ctx.font = "10px 'Plus Jakarta Sans', sans-serif";
  ctx.fillText("Markers", 10, markerY + 4);

  for (let i = 0; i < 50; i++) {
    const bp = Math.random() * (viewEnd - viewStart) + viewStart;
    const x = bpToX(bp);
    ctx.fillStyle = "rgba(245,158,11,0.35)";
    ctx.beginPath();
    ctx.arc(x, markerY + 15, 3, 0, Math.PI * 2);
    ctx.fill();
  }

  // -- Coverage track (simulated) --
  const coverageY = markerY + 40;
  ctx.fillStyle = "rgba(255,255,255,0.25)";
  ctx.font = "10px 'Plus Jakarta Sans', sans-serif";
  ctx.fillText("Signal", 10, coverageY + 4);

  ctx.beginPath();
  ctx.moveTo(0, coverageY + 45);
  for (let x = 0; x < W; x += 3) {
    const bp = viewStart + (x / W) * (viewEnd - viewStart);
    let val = 0;
    // Simulate peaks near genes
    chr14Genes.forEach(g => {
      const mid = (g.start + g.end) / 2;
      const dist = Math.abs(bp - mid) / 500000;
      if (dist < 3) val += Math.exp(-dist * dist) * 30;
    });
    val += Math.random() * 3;
    ctx.lineTo(x, coverageY + 45 - val);
  }
  ctx.lineTo(W, coverageY + 45);
  ctx.closePath();
  const coverageGrad = ctx.createLinearGradient(0, coverageY + 10, 0, coverageY + 45);
  coverageGrad.addColorStop(0, "rgba(6,182,212,0.3)");
  coverageGrad.addColorStop(1, "rgba(6,182,212,0.02)");
  ctx.fillStyle = coverageGrad;
  ctx.fill();

  // Store gene items for hover detection
  canvas._geneItems = geneItems;
}

function formatBp(bp) {
  if (bp >= 1000000) return (bp / 1000000).toFixed(1) + " Mb";
  if (bp >= 1000) return (bp / 1000).toFixed(0) + " kb";
  return bp + " bp";
}

// ============================================================
//   PROCESSING ANIMATION
// ============================================================
const PROCESSING_STEPS = [
  { text: "Parsing coordinates & flanking window parameters...", duration: 900 },
  { text: "Querying Animal QTLdb (Cattle, Goat, Sheep, Pig, etc.)...", duration: 1100 },
  { text: "Mapping Ensembl/NCBI GTF gene models & biotypes...", duration: 1200 },
  { text: "Computing overlap_bp, overlap_pct, and distance_to_tss...", duration: 1000 },
  { text: "Generating QGAT interactive genomic tracks & candidate report...", duration: 800 },
];

function runProcessingAnimation() {
  return new Promise(resolve => {
    const overlay = document.getElementById("processingOverlay");
    const textEl = document.getElementById("processingText");
    const stepEl = document.getElementById("processingStep");
    const barEl = document.getElementById("processingBar");
    overlay.classList.add("active");

    let stepIndex = 0;
    function nextStep() {
      if (stepIndex >= PROCESSING_STEPS.length) {
        overlay.classList.remove("active");
        resolve();
        return;
      }
      const step = PROCESSING_STEPS[stepIndex];
      textEl.textContent = step.text;
      stepEl.textContent = `Step ${stepIndex + 1} of ${PROCESSING_STEPS.length}`;
      barEl.style.width = `${((stepIndex + 1) / PROCESSING_STEPS.length) * 100}%`;
      stepIndex++;
      setTimeout(nextStep, step.duration);
    }
    nextStep();
  });
}

// ============================================================
//   RESULTS ANIMATIONS
// ============================================================
function animateResults() {
  // Animate count-up for summary cards
  document.querySelectorAll(".summary-card__value[data-count]").forEach(el => {
    const target = parseInt(el.dataset.count);
    if (isNaN(target)) return;
    let current = 0;
    const increment = Math.max(1, Math.ceil(target / 30));
    const timer = setInterval(() => {
      current += increment;
      if (current >= target) {
        current = target;
        clearInterval(timer);
      }
      el.textContent = current.toLocaleString();
    }, 30);
    el.style.animation = "countUp 0.5s var(--ease) both";
  });

  // Populate both tables
  populateGenesTable();
  populateQTLsTable();

  // Populate candidate cards
  populateCandidateCards();

  // Activate the appropriate tab
  switchResultsTab(state.activeResultsTab || "genes");
}

function switchResultsTab(tabName) {
  state.activeResultsTab = tabName;
  document.querySelectorAll(".results-nav-btn").forEach(btn => {
    btn.classList.toggle("active", btn.dataset.resultsTab === tabName);
  });
  document.querySelectorAll(".results-tab-content").forEach(content => {
    content.classList.remove("active");
  });

  const tabContentMap = {
    genes: "resultsGenesView",
    qtls: "resultsQTLsView",
    browser: "resultsBrowserView"
  };

  const activeContent = document.getElementById(tabContentMap[tabName] || "resultsGenesView");
  if (activeContent) activeContent.classList.add("active");

  if (tabName === "browser") {
    setTimeout(drawGenomeBrowser, 50);
  }
}

// ---- Populate List of Genes Table ----
function populateGenesTable(genesToDisplay) {
  const genes = genesToDisplay || MOCK_GENES;
  const tbody = document.getElementById("genesTableBody");
  if (!tbody) return;
  tbody.innerHTML = "";

  genes.forEach(gene => {
    const tr = document.createElement("tr");
    tr.innerHTML = `
      <td><span class="table-gene-name">${gene.symbol}</span></td>
      <td><span class="table-mono">${gene.id}</span></td>
      <td>${gene.chr}</td>
      <td><span class="table-mono">${gene.start.toLocaleString()}–${gene.end.toLocaleString()}</span></td>
      <td><span class="table-badge table-badge--teal">${gene.biotype}</span></td>
      <td><span class="table-mono">${(gene.overlap_bp || (gene.end - gene.start + 1)).toLocaleString()} bp</span></td>
      <td><strong>${gene.overlap_pct || 100}%</strong></td>
      <td><span class="table-mono">${gene.distance_to_tss !== undefined ? gene.distance_to_tss + ' bp' : '0 bp'}</span></td>
      <td><span class="table-mono">${gene.strand || '+'}</span></td>
      <td><button class="btn btn--sm btn--primary" style="padding:4px 10px; font-size:0.75rem; white-space:nowrap;">🧬 Gene Structure</button></td>
    `;
    tr.addEventListener("click", () => openGeneDrawer(gene));
    tbody.appendChild(tr);
  });

  const countBadge = document.getElementById("genesCountBadge");
  if (countBadge) countBadge.textContent = genes.length;
  const totalEl = document.getElementById("totalGenesCount");
  if (totalEl) totalEl.textContent = genes.length;
  const showingEl = document.getElementById("showingGenesCount");
  if (showingEl) showingEl.textContent = `1–${genes.length}`;
}

// ---- Populate List of QTLs Table ----
function populateQTLsTable(qtlsToDisplay) {
  const qtls = qtlsToDisplay || (
    state.selectedAnimal
      ? MOCK_QTLS.filter(q => q.species === state.selectedAnimal).concat(MOCK_QTLS.filter(q => q.species !== state.selectedAnimal))
      : MOCK_QTLS
  );
  const tbody = document.getElementById("qtlsTableBody");
  if (!tbody) return;
  tbody.innerHTML = "";

  qtls.forEach(qtl => {
    const tr = document.createElement("tr");
    tr.innerHTML = `
      <td><span class="table-mono" style="color:var(--lavender); font-weight:700;">${qtl.id}</span></td>
      <td>
        <div style="display:flex; align-items:center; gap:8px;">
          <img src="assets/animals/${qtl.species}.jpg" alt="${qtl.species}" class="table-animal-avatar">
          <span>${qtl.speciesName || qtl.species}</span>
        </div>
      </td>
      <td>${qtl.chr}</td>
      <td><span class="table-mono">${qtl.start.toLocaleString()}–${qtl.end.toLocaleString()}</span></td>
      <td><strong>${qtl.trait}</strong></td>
      <td><span class="table-badge table-badge--purple">${qtl.category}</span></td>
      <td><strong>${qtl.overlap_pct || 95}%</strong></td>
      <td><span class="table-mono" style="color:var(--teal); font-weight:600;">${qtl.pubmed}</span></td>
      <td><button class="btn btn--sm btn--outline" style="padding:4px 10px; font-size:0.75rem;">Inspect</button></td>
    `;
    tr.addEventListener("click", () => openQTLDrawer(qtl));
    tbody.appendChild(tr);
  });

  const countBadge = document.getElementById("qtlsCountBadge");
  if (countBadge) countBadge.textContent = qtls.length;
  const totalEl = document.getElementById("totalQTLsCount");
  if (totalEl) totalEl.textContent = qtls.length;
  const showingEl = document.getElementById("showingQTLsCount");
  if (showingEl) showingEl.textContent = `1–${qtls.length}`;

  const spLabel = document.getElementById("statAnimalSpecies");
  if (spLabel) {
    const activeCard = document.querySelector(".animal-select-card.selected");
    spLabel.textContent = activeCard ? activeCard.querySelector(".animal-select-card__name").textContent : "Cattle";
  }
}

function populateCandidateCards() {
  const grid = document.getElementById("candidatesGrid");
  if (!grid) return;
  grid.innerHTML = "";
  const candidates = MOCK_GENES.filter(g => g.score >= 75).slice(0, 6);
  candidates.forEach(gene => {
    const circumference = 2 * Math.PI * 17;
    const offset = circumference - (gene.score / 100) * circumference;
    const strokeColor = gene.score >= 90 ? "var(--teal)" : gene.score >= 80 ? "var(--turquoise)" : "var(--lavender)";

    const evidence = [
      { label: "QTL overlap", has: gene.score > 60 },
      { label: "Functional relevance", has: gene.score > 70 },
      { label: "Expression support", has: gene.score > 75 },
      { label: "Trait association", has: gene.score > 80 },
      { label: "Conserved region", has: gene.score > 85 },
    ];

    const card = document.createElement("div");
    card.className = "candidate-card";
    card.innerHTML = `
      <div class="candidate-card__top">
        <div class="candidate-card__gene">${gene.symbol}</div>
        <div class="candidate-card__score">
          <div class="score-ring">
            <svg viewBox="0 0 40 40">
              <circle class="score-ring__bg" cx="20" cy="20" r="17"></circle>
              <circle class="score-ring__fill" cx="20" cy="20" r="17"
                stroke="${strokeColor}"
                stroke-dasharray="${circumference}"
                stroke-dashoffset="${offset}"></circle>
            </svg>
            <div class="score-ring__text">${gene.score}</div>
          </div>
        </div>
      </div>
      <div class="candidate-card__location">${gene.chr}: ${gene.start.toLocaleString()}–${gene.end.toLocaleString()}</div>
      <div class="candidate-card__evidence">
        ${evidence.map(e => `
          <div class="evidence-item">
            <div class="evidence-item__check evidence-item__check--${e.has ? 'yes' : 'no'}">
              ${e.has ? '✓' : '–'}
            </div>
            ${e.label}
          </div>
        `).join("")}
      </div>
    `;
    card.addEventListener("click", () => openGeneDrawer(gene));
    grid.appendChild(card);
  });
}

// ============================================================
//   REAL GENE STRUCTURE VISUALIZER (Exon-Intron Architecture)
// ============================================================
function generateGeneStructureSVG(gene) {
  if (!gene.exons || !gene.exons.length) {
    return `<div style="padding:1rem; color:var(--text-muted); font-size:0.85rem;">Single exon transcript model</div>`;
  }

  const geneSpan = Math.max(1, gene.end - gene.start + 1);
  const svgW = 700;
  const svgH = 110;
  const padX = 40;
  const usableW = svgW - (padX * 2);
  const midY = 48;

  function toX(bp) {
    return padX + ((bp - gene.start) / geneSpan) * usableW;
  }

  // Splicing chevron lines
  let intronPaths = "";
  for (let i = 0; i < gene.exons.length - 1; i++) {
    const e1 = gene.exons[i];
    const e2 = gene.exons[i + 1];
    const x1 = toX(e1.end);
    const x2 = toX(e2.start);
    const mx = (x1 + x2) / 2;
    const peakY = (i % 2 === 0) ? midY - 14 : midY + 14;
    intronPaths += `<path d="M ${x1.toFixed(1)} ${midY} L ${mx.toFixed(1)} ${peakY} L ${x2.toFixed(1)} ${midY}" fill="none" stroke="#94a3b8" stroke-width="1.8" />`;
  }

  // Exon rects and labels
  let exonRects = "";
  let exonLabels = "";
  let mutationPins = "";

  gene.exons.forEach(ex => {
    const x = toX(ex.start);
    const rawW = toX(ex.end) - x;
    const w = Math.max(8, rawW);
    const isUtr = ex.type === "5_prime_utr" || ex.type === "3_prime_utr";
    const h = isUtr ? 16 : 28;
    const y = midY - (h / 2);
    const fill = isUtr ? "#06b6d4" : "#0d9488";

    exonRects += `
      <g class="exon-element" style="cursor:pointer;">
        <rect x="${x.toFixed(1)}" y="${y.toFixed(1)}" width="${w.toFixed(1)}" height="${h}" rx="3" fill="${fill}" stroke="#0f172a" stroke-width="1" opacity="0.95">
          <title>Exon ${ex.num}: ${ex.start.toLocaleString()}–${ex.end.toLocaleString()} (${ex.len} bp, ${ex.type.toUpperCase()})</title>
        </rect>
      </g>
    `;

    if (w >= 16 && !isUtr) {
      exonLabels += `<text x="${(x + w / 2).toFixed(1)}" y="${(midY + 4).toFixed(1)}" text-anchor="middle" font-size="9" font-weight="700" fill="#ffffff" pointer-events="none">E${ex.num}</text>`;
    }

    if (ex.is_causal) {
      const pinX = x + w / 2;
      mutationPins += `
        <g>
          <line x1="${pinX.toFixed(1)}" y1="${(midY - 16).toFixed(1)}" x2="${pinX.toFixed(1)}" y2="${(midY - 30).toFixed(1)}" stroke="#f43f5e" stroke-width="2" />
          <circle cx="${pinX.toFixed(1)}" cy="${(midY - 32).toFixed(1)}" r="5" fill="#f43f5e" />
          <text x="${pinX.toFixed(1)}" y="${(midY - 40).toFixed(1)}" text-anchor="middle" font-size="9" font-weight="800" fill="#f43f5e">★ ${ex.variant_note || 'QTL'}</text>
        </g>
      `;
    }
  });

  const tick1 = gene.start;
  const tick2 = Math.round(gene.start + geneSpan / 2);
  const tick3 = gene.end;

  return `
    <svg viewBox="0 0 ${svgW} ${svgH}" width="100%" height="${svgH}" xmlns="http://www.w3.org/2000/svg">
      <line x1="${padX}" y1="${midY}" x2="${svgW - padX}" y2="${midY}" stroke="#cbd5e1" stroke-width="1.5" />
      <path d="M ${svgW - padX + 6} ${midY} L ${svgW - padX} ${midY - 4} L ${svgW - padX} ${midY + 4} Z" fill="#64748b" />
      <text x="${svgW - padX + 10}" y="${midY + 3}" font-size="9" font-weight="700" fill="#64748b">3'</text>
      <text x="${padX - 18}" y="${midY + 3}" font-size="9" font-weight="700" fill="#64748b">5'</text>
      ${intronPaths}
      ${exonRects}
      ${exonLabels}
      ${mutationPins}
      <line x1="${padX}" y1="${svgH - 12}" x2="${svgW - padX}" y2="${svgH - 12}" stroke="#94a3b8" stroke-width="1" />
      <text x="${padX}" y="${svgH - 2}" font-size="8.5" fill="#64748b" font-family="monospace">${tick1.toLocaleString()} bp</text>
      <text x="${svgW / 2}" y="${svgH - 2}" font-size="8.5" fill="#64748b" font-family="monospace" text-anchor="middle">${tick2.toLocaleString()} bp</text>
      <text x="${svgW - padX}" y="${svgH - 2}" font-size="8.5" fill="#64748b" font-family="monospace" text-anchor="end">${tick3.toLocaleString()} bp</text>
    </svg>
  `;
}

// ============================================================
//   GENE DETAIL DRAWER (With Real Gene Structure & No Ontology)
// ============================================================
function openGeneDrawer(gene) {
  const drawer = document.getElementById("geneDrawer");
  const backdrop = document.getElementById("drawerBackdrop");
  document.getElementById("drawerGeneName").textContent = `Gene: ${gene.symbol}`;
  document.getElementById("drawerGeneId").textContent = `${gene.id} • ${gene.transcript_id || 'Canonical'} • ${gene.biotype}`;

  const exonsTableHtml = (gene.exons && gene.exons.length) ? `
    <div style="max-height:180px; overflow-y:auto; border:1px solid var(--border); border-radius:var(--radius); margin-top:0.75rem;">
      <table class="exon-table">
        <thead>
          <tr>
            <th>Exon</th>
            <th>Genomic Coordinates</th>
            <th>Length</th>
            <th>Type</th>
            <th>Phase</th>
            <th>Splice Signal</th>
          </tr>
        </thead>
        <tbody>
          ${gene.exons.map(ex => `
            <tr>
              <td><strong>E${ex.num}</strong></td>
              <td>${ex.start.toLocaleString()} – ${ex.end.toLocaleString()}</td>
              <td>${ex.len} bp</td>
              <td><span class="table-badge table-badge--${ex.type === 'cds' ? 'teal' : 'sky'}">${ex.type === 'cds' ? 'Coding CDS' : ex.type.toUpperCase()}</span></td>
              <td>${ex.phase >= 0 ? ex.phase : '–'}</td>
              <td><span class="table-mono">GT-AG</span></td>
            </tr>
          `).join("")}
        </tbody>
      </table>
    </div>
  ` : '';

  const body = document.getElementById("drawerBody");
  body.innerHTML = `
    <!-- Real Gene Structure Schematic -->
    <div class="drawer-section">
      <div class="gene-structure-header">
        <div class="gene-structure-title">
          <span>🧬</span> Exon-Intron Gene Architecture
        </div>
        <div class="gene-structure-badges">
          <span class="table-badge table-badge--teal">${gene.exon_count || (gene.exons ? gene.exons.length : 1)} Exons</span>
          <span class="table-badge table-badge--purple">${gene.cds_length ? gene.cds_length + ' bp CDS' : 'Coding'}</span>
          <span class="table-badge table-badge--sky">${gene.protein_length ? gene.protein_length + ' aa' : 'Canonical'}</span>
        </div>
      </div>

      <div class="gene-structure-svg-wrap">
        ${generateGeneStructureSVG(gene)}
      </div>

      <div class="gene-structure-legend">
        <div class="legend-item"><div class="legend-color legend-color--utr"></div> 5'/3' UTR</div>
        <div class="legend-item"><div class="legend-color legend-color--cds"></div> Coding Exon (CDS)</div>
        <div class="legend-item"><div class="legend-color legend-color--intron"></div> Splicing Junction</div>
        <div class="legend-item"><div class="legend-color legend-color--mutation"></div> Functional / Causal QTL Variant</div>
      </div>

      ${exonsTableHtml}
    </div>

    <!-- Genomic Coordinates & Overlap -->
    <div class="drawer-section">
      <div class="drawer-section__title">Genomic Coordinates &amp; Overlap</div>
      <div class="drawer-info-row">
        <span class="drawer-info-row__label">Chromosome</span>
        <span class="drawer-info-row__value mono">${gene.chr}</span>
      </div>
      <div class="drawer-info-row">
        <span class="drawer-info-row__label">Gene Span</span>
        <span class="drawer-info-row__value mono">${gene.start.toLocaleString()} – ${gene.end.toLocaleString()} (${((gene.end - gene.start + 1)/1000).toFixed(1)} kb)</span>
      </div>
      <div class="drawer-info-row">
        <span class="drawer-info-row__label">Overlap Length</span>
        <span class="drawer-info-row__value mono">${(gene.overlap_bp || (gene.end - gene.start + 1)).toLocaleString()} bp (${gene.overlap_pct || 100}%)</span>
      </div>
      <div class="drawer-info-row">
        <span class="drawer-info-row__label">Distance to TSS</span>
        <span class="drawer-info-row__value mono">${gene.distance_to_tss || 0} bp</span>
      </div>
      <div class="drawer-info-row">
        <span class="drawer-info-row__label">Strand</span>
        <span class="drawer-info-row__value mono">${gene.strand || '+'}</span>
      </div>
      <div class="drawer-info-row">
        <span class="drawer-info-row__label">Ensembl Transcript</span>
        <span class="drawer-info-row__value mono">${gene.transcript_id || 'Canonical'}</span>
      </div>
    </div>

    <!-- Functional Description & Known Variants -->
    <div class="drawer-section">
      <div class="drawer-section__title">Functional Description &amp; Causal Variants</div>
      <p style="font-size:0.88rem; color:var(--text-secondary); line-height:1.7; margin-bottom:0.75rem;">
        ${gene.function}
      </p>
      ${gene.key_variant ? `
        <div style="background:rgba(244,63,94,0.06); border:1px solid rgba(244,63,94,0.25); border-radius:var(--radius); padding:10px 12px; margin-top:0.5rem; font-size:0.82rem;">
          <strong style="color:#e11d48;">★ Causal / Key Diagnostic Variant:</strong><br>
          <span style="color:var(--navy); font-weight:600;">${gene.key_variant}</span>
        </div>
      ` : ''}
    </div>

    <!-- Animal QTLdb Overlap -->
    <div class="drawer-section">
      <div class="drawer-section__title">Animal QTLdb Association</div>
      <div class="drawer-info-row">
        <span class="drawer-info-row__label">Associated QTL</span>
        <span class="drawer-info-row__value"><strong>${gene.qtl}</strong></span>
      </div>
      <div class="drawer-info-row">
        <span class="drawer-info-row__label">Agricultural Trait</span>
        <span class="drawer-info-row__value">${gene.trait}</span>
      </div>
      <div class="drawer-info-row">
        <span class="drawer-info-row__label">Candidate Confidence</span>
        <span class="drawer-info-row__value"><strong style="color:var(--teal);">${gene.score}/100</strong></span>
      </div>
    </div>

    <!-- Biological Pathways -->
    <div class="drawer-section">
      <div class="drawer-section__title">Biological Pathways</div>
      <div style="display:flex; flex-wrap:wrap; gap:6px;">
        ${(gene.pathways || []).map(p => `<span class="drawer-tag drawer-tag--purple">${p}</span>`).join("")}
      </div>
    </div>
  `;

  drawer.classList.add("open");
  backdrop.classList.add("visible");
}

// ============================================================
//   QTL DETAIL DRAWER (With Real Livestock Photo Hero)
// ============================================================
function openQTLDrawer(qtl) {
  const drawer = document.getElementById("geneDrawer");
  const backdrop = document.getElementById("drawerBackdrop");
  document.getElementById("drawerGeneName").textContent = `${qtl.id}: ${qtl.trait}`;
  document.getElementById("drawerGeneId").textContent = `${qtl.speciesName} • Animal QTLdb`;

  const body = document.getElementById("drawerBody");
  body.innerHTML = `
    <!-- Livestock Photo Hero Banner -->
    <div class="qtl-drawer-hero">
      <img src="assets/animals/${qtl.species}.jpg" alt="${qtl.speciesName}">
      <div class="qtl-drawer-hero-overlay">
        <span class="table-badge table-badge--teal" style="align-self:flex-start; margin-bottom:4px;">${qtl.speciesName}</span>
        <h3 style="margin:0; font-size:1.05rem; color:#fff; font-weight:700;">${qtl.id} • ${qtl.trait}</h3>
        <span style="font-size:0.75rem; opacity:0.85;">Pre-Packaged Livestock Animal QTLdb</span>
      </div>
    </div>

    <div class="drawer-section">
      <div class="drawer-section__title">Animal QTLdb Coordinates</div>
      <div class="drawer-info-row">
        <span class="drawer-info-row__label">Species</span>
        <span class="drawer-info-row__value">${qtl.speciesName}</span>
      </div>
      <div class="drawer-info-row">
        <span class="drawer-info-row__label">Chromosome</span>
        <span class="drawer-info-row__value mono">${qtl.chr}</span>
      </div>
      <div class="drawer-info-row">
        <span class="drawer-info-row__label">QTL Span</span>
        <span class="drawer-info-row__value mono">${qtl.start.toLocaleString()} – ${qtl.end.toLocaleString()} (${((qtl.end - qtl.start + 1) / 1000).toFixed(1)} kb)</span>
      </div>
      <div class="drawer-info-row">
        <span class="drawer-info-row__label">Region Overlap</span>
        <span class="drawer-info-row__value"><strong>${qtl.overlap_pct}%</strong></span>
      </div>
    </div>

    <div class="drawer-section">
      <div class="drawer-section__title">Trait Classification</div>
      <div class="drawer-info-row">
        <span class="drawer-info-row__label">Trait Name</span>
        <span class="drawer-info-row__value"><strong>${qtl.trait}</strong></span>
      </div>
      <div class="drawer-info-row">
        <span class="drawer-info-row__label">Trait Category</span>
        <span class="drawer-info-row__value"><span class="table-badge table-badge--purple">${qtl.category}</span></span>
      </div>
      <p style="font-size:0.88rem; color:var(--text-secondary); line-height:1.7; margin-top:0.75rem;">
        ${qtl.desc}
      </p>
    </div>

    <div class="drawer-section">
      <div class="drawer-section__title">Candidate Gene Mapping</div>
      <div class="drawer-info-row">
        <span class="drawer-info-row__label">Causal / Candidate Gene</span>
        <span class="drawer-info-row__value"><strong style="color:var(--teal); font-size:1.05rem;">${qtl.candidate}</strong></span>
      </div>
    </div>

    <div class="drawer-section">
      <div class="drawer-section__title">Literature &amp; External Database Links</div>
      <div class="drawer-info-row">
        <span class="drawer-info-row__label">PubMed Reference</span>
        <span class="drawer-info-row__value"><a href="https://pubmed.ncbi.nlm.nih.gov/${qtl.pubmed.replace('PMID:','')}/" target="_blank" style="color:var(--teal); text-decoration:underline; font-weight:600;">${qtl.pubmed} ↗</a></span>
      </div>
      <div class="drawer-info-row">
        <span class="drawer-info-row__label">Animal QTLdb Source</span>
        <span class="drawer-info-row__value"><a href="https://www.animalgenome.org/cgi-bin/QTLdb/index" target="_blank" style="color:var(--teal); text-decoration:underline;">animalgenome.org/QTLdb ↗</a></span>
      </div>
    </div>
  `;

  drawer.classList.add("open");
  backdrop.classList.add("visible");
}

function closeGeneDrawer() {
  document.getElementById("geneDrawer").classList.remove("open");
  document.getElementById("drawerBackdrop").classList.remove("visible");
}

// ============================================================
//   GENE ONTOLOGY & FUNCTIONAL ENRICHMENT TOOLS
// ============================================================
const ORGANISM_CONFIG = {
  btaurus: { name: "Bos taurus (Cattle)", gprofiler: "btaurus", david: "Bos taurus", assem: "ARS-UCD1.2" },
  chircus: { name: "Capra hircus (Goat)", gprofiler: "chircus", david: "Capra hircus", assem: "CHIR_1.0" },
  oaries: { name: "Ovis aries (Sheep)", gprofiler: "oaries", david: "Ovis aries", assem: "OAR3.1" },
  sscrofa: { name: "Sus scrofa (Pig)", gprofiler: "sscrofa", david: "Sus scrofa", assem: "Sscrofa11.1" },
  ggallus: { name: "Gallus gallus (Chicken)", gprofiler: "ggallus", david: "Gallus gallus", assem: "GRCg4" },
  ecaballus: { name: "Equus caballus (Horse)", gprofiler: "ecaballus", david: "Equus caballus", assem: "EquCab2" }
};

const CITATION_DATA = {
  gprofiler: "Raudvere U, Kolberg L, Kuzmin I, et al. (2019). g:Profiler: a web server for functional enrichment analysis and conversions of gene lists (2019 update). Nucleic Acids Res, 47(W1):W191–W198. DOI: 10.1093/nar/gkz369",
  david: "Sherman BT, Hao M, Qiu J, et al. (2022). DAVID: a web server for functional enrichment analysis and functional annotation of gene lists (2021 update). Nucleic Acids Res, 50(W1):W216–W221. DOI: 10.1093/nar/gkac194",
  panther: "Thomas PD, Ebert D, Muruganujan A, et al. (2022). PANTHER: Making genome-scale phylogenetics accessible to all. Protein Sci, 31(1):8–22. DOI: 10.1002/pro.4218",
  shinygo: "Ge SX, Jung D, Yao R. (2020). ShinyGO: a graphical gene-set enrichment tool for animals and plants. Bioinformatics, 36(8):2628–2629. DOI: 10.1093/bioinformatics/btz931"
};

function getOntologyGeneSymbols() {
  const textarea = document.getElementById("ontologyGenesTextarea");
  if (!textarea || !textarea.value.trim()) {
    return (state.ontologyGenes && state.ontologyGenes.length)
      ? state.ontologyGenes.map(g => g.symbol)
      : MOCK_GENES.map(g => g.symbol);
  }
  return textarea.value.trim().split(/[\s,;\n\t]+/).filter(Boolean);
}

function updateOntologyView() {
  const textarea = document.getElementById("ontologyGenesTextarea");
  const countBadge = document.getElementById("ontologyGeneCountBadge");

  if (textarea && (!textarea.value || textarea.value.trim() === "")) {
    const genes = (state.ontologyGenes && state.ontologyGenes.length)
      ? state.ontologyGenes
      : MOCK_GENES;
    textarea.value = genes.map(g => g.symbol).join(" ");
  }

  if (countBadge) {
    const list = getOntologyGeneSymbols();
    countBadge.textContent = `${list.length} gene${list.length === 1 ? '' : 's'} loaded`;
  }
}

function launchGProfiler() {
  const genes = getOntologyGeneSymbols();
  if (!genes.length) {
    alert("Please enter or load at least one gene symbol first.");
    return;
  }
  const spKey = (state.selectedAnimal || "cattle").toLowerCase();
  const matchedOrg = Object.keys(ORGANISM_CONFIG).find(k => k.startsWith(spKey.slice(0, 3))) || "btaurus";
  const orgInfo = ORGANISM_CONFIG[matchedOrg] || ORGANISM_CONFIG.btaurus;

  const url = `https://biit.cs.ut.ee/gprofiler/gost?organism=${orgInfo.gprofiler}&query=${encodeURIComponent(genes.join(" "))}`;
  window.open(url, "_blank");
}

function launchDAVID() {
  const genes = getOntologyGeneSymbols();
  if (!genes.length) {
    alert("Please enter or load at least one gene symbol first.");
    return;
  }
  navigator.clipboard.writeText(genes.join("\n")).then(() => {
    alert(`Copied ${genes.length} gene symbols to your clipboard!\n\nOpening DAVID Bioinformatics Resources...\n1. In DAVID, paste your genes into Box 1.\n2. Select "OFFICIAL_GENE_SYMBOL" in Box 2.\n3. Select "Gene List" in Box 3 and click Submit List.`);
    window.open("https://david.ncifcrf.gov/summary.jsp", "_blank");
  }).catch(() => {
    window.open("https://david.ncifcrf.gov/summary.jsp", "_blank");
  });
}

function launchPanther() {
  const genes = getOntologyGeneSymbols();
  if (!genes.length) {
    alert("Please enter or load at least one gene symbol first.");
    return;
  }
  navigator.clipboard.writeText(genes.join("\n")).then(() => {
    alert(`Copied ${genes.length} gene symbols to clipboard!\n\nOpening PANTHER Classification System...\nPaste your genes into the query box and select your organism.`);
    window.open("http://www.pantherdb.org/tools/compareToRefList.jsp", "_blank");
  }).catch(() => {
    window.open("http://www.pantherdb.org/tools/compareToRefList.jsp", "_blank");
  });
}

function launchShinyGO() {
  const genes = getOntologyGeneSymbols();
  if (!genes.length) {
    alert("Please enter or load at least one gene symbol first.");
    return;
  }
  navigator.clipboard.writeText(genes.join("\n")).then(() => {
    alert(`Copied ${genes.length} gene symbols to clipboard!\n\nOpening ShinyGO v0.80...\nPaste your gene list and choose species.`);
    window.open("http://bioinformatics.sdstate.edu/go/", "_blank");
  }).catch(() => {
    window.open("http://bioinformatics.sdstate.edu/go/", "_blank");
  });
}

function copyCitation(toolKey, btnEl) {
  const text = CITATION_DATA[toolKey];
  if (!text) return;
  navigator.clipboard.writeText(text).then(() => {
    if (btnEl) {
      const orig = btnEl.textContent;
      btnEl.textContent = "✓ Citation Copied!";
      setTimeout(() => { btnEl.textContent = orig; }, 2000);
    }
  });
}



// ============================================================
//   EVENT HANDLERS
// ============================================================
function init() {
  // ---- Navigation ----
  document.querySelectorAll("[data-nav]").forEach(el => {
    el.addEventListener("click", e => {
      e.preventDefault();
      navigateTo(el.dataset.nav);
    });
  });

  // ---- Scroll-based nav style ----
  window.addEventListener("scroll", () => {
    const nav = document.getElementById("mainNav");
    if (state.currentPage === "home") {
      if (window.scrollY > 80) {
        nav.className = "nav nav--solid";
      } else {
        nav.className = "nav nav--transparent";
      }
    }

    // Animate-on-scroll
    document.querySelectorAll(".animate-on-scroll").forEach(el => {
      const rect = el.getBoundingClientRect();
      if (rect.top < window.innerHeight - 60) {
        el.classList.add("visible");
      }
    });
  });

  // ---- Hero Canvas ----
  const heroCanvas = document.getElementById("heroCanvas");
  if (heroCanvas) {
    const heroAnim = new HeroAnimation(heroCanvas);
    heroAnim.start();
  }

  // ---- Type selector (Gene / QTL / Combined) ----
  document.querySelectorAll(".type-card").forEach(card => {
    card.addEventListener("click", () => {
      const type = card.dataset.type;
      state.annotationType = type;
      document.querySelectorAll(".type-card").forEach(c => c.classList.remove("selected"));
      card.classList.add("selected");
      document.querySelectorAll(".input-panel").forEach(p => p.classList.remove("active"));
      const panelMap = {
        gene: "panelGene",
        qtl: "panelQTL",
        combined: "panelCombined"
      };
      const target = document.getElementById(panelMap[type] || "panelGene");
      if (target) target.classList.add("active");
    });
  });

  // ---- Coordinate Format Toggles (Range vs Position + Window) ----
  function setupFormatToggle(toggleId, windowGroupId, textareaId, labelId, isQtl) {
    const toggle = document.getElementById(toggleId);
    if (!toggle) return;
    const windowGroup = document.getElementById(windowGroupId);
    const textarea = document.getElementById(textareaId);
    const label = document.getElementById(labelId);

    toggle.querySelectorAll(".format-toggle__btn").forEach(btn => {
      btn.addEventListener("click", () => {
        toggle.querySelectorAll(".format-toggle__btn").forEach(b => b.classList.remove("active"));
        btn.classList.add("active");
        const format = btn.dataset.format;

        if (isQtl) {
          state.qtlFormat = format;
        } else if (toggleId === "combinedFormatToggle") {
          state.combinedFormat = format;
        } else {
          state.geneFormat = format;
        }

        if (format === "position") {
          if (windowGroup) windowGroup.style.display = "block";
          if (label) {
            label.innerHTML = isQtl
              ? 'Genomic Coordinates <span class="form-sublabel">(Chromosome &lt;tab/space&gt; Position [± Window Applied])</span>'
              : 'Paste Coordinates <span class="form-sublabel">(Chromosome &lt;tab/space&gt; Position [± Window Applied])</span>';
          }
          if (textarea && !textarea.value.trim()) {
            textarea.placeholder = "14\t1801116\n6\t38024921\n2\t6445100\n20\t31920050\n26\t21150000";
          }
        } else {
          if (windowGroup) windowGroup.style.display = "none";
          if (label) {
            label.innerHTML = isQtl
              ? 'Genomic Coordinates <span class="form-sublabel">(Chromosome &lt;tab/space&gt; Start &lt;tab/space&gt; End)</span>'
              : 'Paste Coordinates <span class="form-sublabel">(Chromosome &lt;tab/space&gt; Start &lt;tab/space&gt; End)</span>';
          }
          if (textarea && !textarea.value.trim()) {
            textarea.placeholder = isQtl
              ? "14\t1500000\t3200000\n6\t37000000\t42000000\n2\t6000000\t7000000\n20\t31500000\t32500000\n26\t21100000\t21200000"
              : "14\t1790000\t1820000\n6\t37900000\t38100000\n2\t6440000\t6460000\n20\t31880000\t32180000\n26\t21140000\t21170000";
          }
        }
      });
    });
  }

  setupFormatToggle("geneFormatToggle", "geneWindowGroup", "geneCoordsInput", "geneInputLabel", false);
  setupFormatToggle("qtlFormatToggle", "qtlWindowGroup", "qtlCoordsInput", "qtlInputLabel", true);
  setupFormatToggle("combinedFormatToggle", "combinedWindowGroup", "combinedCoordsInput", null, false);

  // ---- Window Preset Chips ----
  function setupWindowChips(containerId, inputId) {
    const container = document.getElementById(containerId);
    const input = document.getElementById(inputId);
    if (!container || !input) return;
    container.querySelectorAll(".chip").forEach(chip => {
      chip.addEventListener("click", () => {
        container.querySelectorAll(".chip").forEach(c => c.classList.remove("active"));
        chip.classList.add("active");
        if (chip.dataset.window) {
          input.value = chip.dataset.window;
        }
      });
    });
  }
  setupWindowChips("geneWindowPresetChips", "geneCustomWindowInput");
  setupWindowChips("qtlWindowPresetChips", "qtlCustomWindowInput");

  // ---- Chips generic toggle ----
  document.querySelectorAll(".chip-group:not(#geneWindowPresetChips):not(#qtlWindowPresetChips) .chip").forEach(chip => {
    chip.addEventListener("click", () => chip.classList.toggle("active"));
  });

  // ---- Animal Selection Cards (Animal QTLdb Livestock) ----
  const ANIMAL_PRESETS = {
    cattle: {
      range: "14\t1750000\t1850000\n6\t37900000\t38200000\n2\t6400000\t6500000\n20\t31800000\t32200000\n26\t21100000\t21200000",
      pos: "14\t1801116\n6\t38024921\n2\t6445100\n20\t31920050\n26\t21150000"
    },
    goat: {
      range: "6\t86000000\t87500000\n3\t15000000\t16500000\n11\t42000000\t43500000",
      pos: "6\t87145000\n3\t15800000\n11\t42900000"
    },
    sheep: {
      range: "X\t48000000\t49500000\n6\t35000000\t36500000\n2\t118000000\t120000000",
      pos: "X\t48750000\n6\t35800000\n2\t119000000"
    },
    pig: {
      range: "2\t1400000\t1600000\n6\t49000000\t51000000\n4\t62000000\t64000000",
      pos: "2\t1500000\n6\t50200000\n4\t63100000"
    },
    chicken: {
      range: "4\t72000000\t74000000\n1\t54000000\t56000000\n3\t31000000\t33000000",
      pos: "4\t73000000\n1\t55100000\n3\t32000000"
    },
    horse: {
      range: "18\t66000000\t67000000\n3\t104000000\t106000000\n9\t45000000\t47000000",
      pos: "18\t66500000\n3\t105000000\n9\t46200000"
    }
  };

  document.querySelectorAll("#qtlAnimalGrid .animal-select-card").forEach(card => {
    card.addEventListener("click", () => {
      document.querySelectorAll("#qtlAnimalGrid .animal-select-card").forEach(c => c.classList.remove("selected"));
      card.classList.add("selected");
      state.selectedAnimal = card.dataset.species;

      // Update animal species indicator
      const spLabel = document.getElementById("statAnimalSpecies");
      if (spLabel) spLabel.textContent = card.querySelector(".animal-select-card__name").textContent;

      // If coordinates input has text or is empty, we can update placeholder
      const qtlInput = document.getElementById("qtlCoordsInput");
      if (qtlInput && ANIMAL_PRESETS[state.selectedAnimal]) {
        qtlInput.placeholder = (state.qtlFormat === "position")
          ? ANIMAL_PRESETS[state.selectedAnimal].pos
          : ANIMAL_PRESETS[state.selectedAnimal].range;
      }
    });
  });

  // ---- Preset & Example Fill Buttons ----
  const loadGeneBtn = document.getElementById("loadGeneRangePresetBtn");
  if (loadGeneBtn) {
    loadGeneBtn.addEventListener("click", () => {
      const input = document.getElementById("geneCoordsInput");
      if (input) {
        if (state.geneFormat === "position") {
          input.value = "14\t1800000\n6\t38000000\n2\t6445000\n20\t32000000\n26\t21150000\n4\t93270000";
        } else {
          input.value = "14\t1790000\t1820000\n6\t37900000\t38100000\n2\t6440000\t6460000\n20\t31880000\t32180000\n26\t21140000\t21170000\n4\t93260000\t93290000";
        }
      }
    });
  }

  const loadQtlBtn = document.getElementById("loadQTLRangePresetBtn");
  if (loadQtlBtn) {
    loadQtlBtn.addEventListener("click", () => {
      const input = document.getElementById("qtlCoordsInput");
      const sp = state.selectedAnimal || "cattle";
      const preset = ANIMAL_PRESETS[sp] || ANIMAL_PRESETS.cattle;
      if (input) {
        input.value = (state.qtlFormat === "position") ? preset.pos : preset.range;
      }
    });
  }

  const loadCombinedBtn = document.getElementById("loadCombinedPresetBtn");
  if (loadCombinedBtn) {
    loadCombinedBtn.addEventListener("click", () => {
      const input = document.getElementById("combinedCoordsInput");
      if (input) {
        if (state.combinedFormat === "position") {
          input.value = "14\t1800000\n6\t38000000\n2\t6445000\n20\t32000000";
        } else {
          input.value = "14\t1750000\t1850000\n6\t37900000\t38200000\n2\t6400000\t6500000\n20\t31800000\t32200000";
        }
      }
    });
  }

  function loadGlobalExample() {
    const geneInput = document.getElementById("geneCoordsInput");
    if (geneInput) {
      geneInput.value = "14\t1790000\t1820000\n6\t37900000\t38100000\n2\t6440000\t6460000\n20\t31880000\t32180000\n26\t21140000\t21170000";
    }
    const qtlInput = document.getElementById("qtlCoordsInput");
    if (qtlInput) {
      qtlInput.value = "14\t1750000\t1850000\n6\t37900000\t38200000\n2\t6400000\t6500000\n20\t31800000\t32200000\n26\t21100000\t21200000";
    }
    const combInput = document.getElementById("combinedCoordsInput");
    if (combInput) {
      combInput.value = "14\t1750000\t1850000\n6\t37900000\t38200000\n2\t6400000\t6500000\n20\t31800000\t32200000";
    }
    navigateTo("annotate");
  }
  const globalExBtn = document.getElementById("loadExampleBtn");
  if (globalExBtn) globalExBtn.addEventListener("click", loadGlobalExample);
  const heroExBtn = document.getElementById("heroExampleBtn");
  if (heroExBtn) heroExBtn.addEventListener("click", loadGlobalExample);

  // ---- File uploads setup ----
  function setupUpload(zoneId, inputId, infoId, nameId) {
    const zone = document.getElementById(zoneId);
    const input = document.getElementById(inputId);
    const info = document.getElementById(infoId);
    const nameEl = document.getElementById(nameId);
    if (!zone || !input || !info || !nameEl) return;

    zone.addEventListener("click", () => input.click());
    zone.addEventListener("dragover", e => { e.preventDefault(); zone.classList.add("drag-over"); });
    zone.addEventListener("dragleave", () => zone.classList.remove("drag-over"));
    zone.addEventListener("drop", e => {
      e.preventDefault();
      zone.classList.remove("drag-over");
      if (e.dataTransfer.files.length) handleFile(e.dataTransfer.files[0]);
    });
    input.addEventListener("change", () => {
      if (input.files.length) handleFile(input.files[0]);
    });

    function handleFile(file) {
      nameEl.textContent = `${file.name} (${(file.size / 1024).toFixed(1)} KB)`;
      info.classList.add("visible");
    }
  }
  setupUpload("geneCoordsUploadZone", "geneCoordsFileInput", "geneCoordsFileInfo", "geneCoordsFileName");
  setupUpload("qtlUploadZone", "qtlFileInput", "qtlFileInfo", "qtlFileName");

  // ---- Execute Gene Annotation (Produces List of Genes from GTF) ----
  const runGene = document.getElementById("runGeneAnnotation");
  if (runGene) {
    runGene.addEventListener("click", async () => {
      const input = document.getElementById("geneCoordsInput");
      if (input && !input.value.trim()) {
        if (state.geneFormat === "position") {
          input.value = "14\t1800000\n6\t38000000\n2\t6445000\n20\t32000000\n26\t21150000";
        } else {
          input.value = "14\t1790000\t1820000\n6\t37900000\t38100000\n2\t6440000\t6460000\n20\t31880000\t32180000\n26\t21140000\t21170000";
        }
      }
      state.activeResultsTab = "genes";
      state.resultsLoaded = false;
      await runProcessingAnimation();
      navigateTo("results");
      switchResultsTab("genes");
    });
  }

  // ---- Execute QTL Annotation (Produces List of QTLs from Animal QTLdb) ----
  const runQTL = document.getElementById("runQTLAnnotation");
  if (runQTL) {
    runQTL.addEventListener("click", async () => {
      const input = document.getElementById("qtlCoordsInput");
      const sp = state.selectedAnimal || "cattle";
      const preset = ANIMAL_PRESETS[sp] || ANIMAL_PRESETS.cattle;
      if (input && !input.value.trim()) {
        input.value = (state.qtlFormat === "position") ? preset.pos : preset.range;
      }
      state.activeResultsTab = "qtls";
      state.resultsLoaded = false;
      await runProcessingAnimation();
      navigateTo("results");
      switchResultsTab("qtls");
    });
  }

  // ---- Execute Combined Annotation ----
  const runCombined = document.getElementById("runCombinedAnnotation");
  if (runCombined) {
    runCombined.addEventListener("click", async () => {
      state.activeResultsTab = "genes";
      state.resultsLoaded = false;
      await runProcessingAnimation();
      navigateTo("results");
      switchResultsTab("genes");
    });
  }

  // ---- Results Sub-Navigation Tabs ----
  document.querySelectorAll(".results-nav-btn").forEach(btn => {
    btn.addEventListener("click", () => {
      switchResultsTab(btn.dataset.resultsTab);
    });
  });

  // ---- Table Search / Filtering ----
  const genesSearch = document.getElementById("genesSearchInput");
  if (genesSearch) {
    genesSearch.addEventListener("input", e => {
      const query = e.target.value.toLowerCase();
      const rows = document.querySelectorAll("#genesTableBody tr");
      rows.forEach(row => {
        row.style.display = row.textContent.toLowerCase().includes(query) ? "" : "none";
      });
    });
  }

  const qtlsSearch = document.getElementById("qtlsSearchInput");
  if (qtlsSearch) {
    qtlsSearch.addEventListener("input", e => {
      const query = e.target.value.toLowerCase();
      const rows = document.querySelectorAll("#qtlsTableBody tr");
      rows.forEach(row => {
        row.style.display = row.textContent.toLowerCase().includes(query) ? "" : "none";
      });
    });
  }

  // ---- Export Buttons ----
  const exportGenesCSVBtn = document.getElementById("exportGenesCSV");
  if (exportGenesCSVBtn) {
    exportGenesCSVBtn.addEventListener("click", () => {
      let csv = "Gene Symbol,Gene ID,Chromosome,Start,End,Biotype,Overlap bp,Overlap %,Distance to TSS,Strand\n";
      MOCK_GENES.forEach(g => {
        csv += `${g.symbol},${g.id},${g.chr},${g.start},${g.end},${g.biotype},${g.overlap_bp},${g.overlap_pct}%,${g.distance_to_tss},${g.strand}\n`;
      });
      downloadFile(csv, "qgat_genes_list.csv", "text/csv");
    });
  }

  const exportGenesTSVBtn = document.getElementById("exportGenesTSV");
  if (exportGenesTSVBtn) {
    exportGenesTSVBtn.addEventListener("click", () => {
      let tsv = "Gene Symbol\tGene ID\tChromosome\tStart\tEnd\tBiotype\tOverlap bp\tOverlap %\tDistance to TSS\tStrand\n";
      MOCK_GENES.forEach(g => {
        tsv += `${g.symbol}\t${g.id}\t${g.chr}\t${g.start}\t${g.end}\t${g.biotype}\t${g.overlap_bp}\t${g.overlap_pct}%\t${g.distance_to_tss}\t${g.strand}\n`;
      });
      downloadFile(tsv, "qgat_genes_list.tsv", "text/tab-separated-values");
    });
  }

  const exportQTLsCSVBtn = document.getElementById("exportQTLsCSV");
  if (exportQTLsCSVBtn) {
    exportQTLsCSVBtn.addEventListener("click", () => {
      let csv = "QTL ID,Species,Chromosome,Start,End,Trait,Category,Overlap %,PubMed ID,Candidate Gene\n";
      const qtls = state.selectedAnimal
        ? MOCK_QTLS.filter(q => q.species === state.selectedAnimal).concat(MOCK_QTLS.filter(q => q.species !== state.selectedAnimal))
        : MOCK_QTLS;
      qtls.forEach(q => {
        csv += `${q.id},"${q.speciesName}",${q.chr},${q.start},${q.end},"${q.trait}","${q.category}",${q.overlap_pct}%,${q.pubmed},${q.candidate}\n`;
      });
      downloadFile(csv, `qgat_${state.selectedAnimal || 'animal'}_qtls.csv`, "text/csv");
    });
  }

  const exportQTLsTSVBtn = document.getElementById("exportQTLsTSV");
  if (exportQTLsTSVBtn) {
    exportQTLsTSVBtn.addEventListener("click", () => {
      let tsv = "QTL ID\tSpecies\tChromosome\tStart\tEnd\tTrait\tCategory\tOverlap %\tPubMed ID\tCandidate Gene\n";
      const qtls = state.selectedAnimal
        ? MOCK_QTLS.filter(q => q.species === state.selectedAnimal).concat(MOCK_QTLS.filter(q => q.species !== state.selectedAnimal))
        : MOCK_QTLS;
      qtls.forEach(q => {
        tsv += `${q.id}\t${q.speciesName}\t${q.chr}\t${q.start}\t${q.end}\t${q.trait}\t${q.category}\t${q.overlap_pct}%\t${q.pubmed}\t${q.candidate}\n`;
      });
      downloadFile(tsv, `qgat_${state.selectedAnimal || 'animal'}_qtls.tsv`, "text/tab-separated-values");
    });
  }

  // ---- Drawer Close Handlers ----
  const drawerCloseBtn = document.getElementById("drawerClose");
  if (drawerCloseBtn) drawerCloseBtn.addEventListener("click", closeGeneDrawer);
  const drawerBackdropEl = document.getElementById("drawerBackdrop");
  if (drawerBackdropEl) drawerBackdropEl.addEventListener("click", closeGeneDrawer);

  // ---- Genome Browser Zoom Controls ----
  const zoomInBtn = document.getElementById("zoomIn");
  if (zoomInBtn) {
    zoomInBtn.addEventListener("click", () => {
      state.genomeZoom = Math.min(state.genomeZoom * 1.5, 20);
      drawGenomeBrowser();
    });
  }

  const zoomOutBtn = document.getElementById("zoomOut");
  if (zoomOutBtn) {
    zoomOutBtn.addEventListener("click", () => {
      state.genomeZoom = Math.max(state.genomeZoom / 1.5, 0.5);
      drawGenomeBrowser();
    });
  }

  const resetZoomBtn = document.getElementById("resetZoom");
  if (resetZoomBtn) {
    resetZoomBtn.addEventListener("click", () => {
      state.genomeZoom = 1;
      state.genomePan = 0;
      drawGenomeBrowser();
    });
  }

  // ---- Genome Canvas Interactions ----
  const genomeCanvas = document.getElementById("genomeCanvas");
  const tooltip = document.getElementById("genomeTooltip");

  if (genomeCanvas && tooltip) {
    genomeCanvas.addEventListener("mousemove", e => {
      const rect = genomeCanvas.getBoundingClientRect();
      const mx = e.clientX - rect.left;
      const my = e.clientY - rect.top;
      const items = genomeCanvas._geneItems || [];
      let found = null;

      for (const item of items) {
        if (mx >= item.pixelX && mx <= item.pixelX + item.pixelW &&
            my >= item.pixelY && my <= item.pixelY + item.pixelH) {
          found = item;
          break;
        }
      }

      if (found) {
        tooltip.classList.add("visible");
        const titleEl = document.getElementById("tooltipTitle");
        if (titleEl) titleEl.textContent = `Gene: ${found.symbol}`;
        const contentEl = document.getElementById("tooltipContent");
        if (contentEl) {
          contentEl.innerHTML = `
            <div class="genome-tooltip__row"><span>Location</span><span>${found.chr || 'Chr14'}: ${found.start.toLocaleString()}–${found.end.toLocaleString()}</span></div>
            <div class="genome-tooltip__row"><span>Width</span><span>${(found.end - found.start + 1).toLocaleString()} bp</span></div>
            <div class="genome-tooltip__row"><span>Candidate</span><span>${found.candidate ? 'Yes ★' : 'No'}</span></div>
          `;
        }
        tooltip.style.left = (e.clientX + 12) + "px";
        tooltip.style.top = (e.clientY - 10) + "px";
        genomeCanvas.style.cursor = "pointer";
      } else {
        tooltip.classList.remove("visible");
        genomeCanvas.style.cursor = "crosshair";
      }
    });

    genomeCanvas.addEventListener("mouseleave", () => {
      tooltip.classList.remove("visible");
    });

    genomeCanvas.addEventListener("click", e => {
      const rect = genomeCanvas.getBoundingClientRect();
      const mx = e.clientX - rect.left;
      const my = e.clientY - rect.top;
      const items = genomeCanvas._geneItems || [];
      for (const item of items) {
        if (mx >= item.pixelX && mx <= item.pixelX + item.pixelW &&
            my >= item.pixelY && my <= item.pixelY + item.pixelH) {
          const fullGene = MOCK_GENES.find(g => g.symbol === item.symbol);
          if (fullGene) openGeneDrawer(fullGene);
          break;
        }
      }
    });

    let isDragging = false;
    let lastDragX = 0;
    genomeCanvas.addEventListener("mousedown", e => {
      isDragging = true;
      lastDragX = e.clientX;
    });
    window.addEventListener("mousemove", e => {
      if (!isDragging) return;
      const dx = e.clientX - lastDragX;
      lastDragX = e.clientX;
      const W = genomeCanvas.offsetWidth;
      const viewRange = 50000000 / state.genomeZoom;
      state.genomePan -= dx / W * viewRange;
      state.genomePan = Math.max(0, state.genomePan);
      drawGenomeBrowser();
    });
    window.addEventListener("mouseup", () => { isDragging = false; });

    genomeCanvas.addEventListener("wheel", e => {
      e.preventDefault();
      if (e.deltaY < 0) {
        state.genomeZoom = Math.min(state.genomeZoom * 1.15, 20);
      } else {
        state.genomeZoom = Math.max(state.genomeZoom / 1.15, 0.5);
      }
      drawGenomeBrowser();
    });
  }

  // ---- Gene Ontology & Functional Enrichment Handlers ----
  const ontologyTextarea = document.getElementById("ontologyGenesTextarea");
  if (ontologyTextarea) {
    ontologyTextarea.addEventListener("input", () => {
      const list = getOntologyGeneSymbols();
      const badge = document.getElementById("ontologyGeneCountBadge");
      if (badge) badge.textContent = `${list.length} gene${list.length === 1 ? '' : 's'} loaded`;
    });
  }

  const btnCopySymbols = document.getElementById("btnCopySymbols");
  if (btnCopySymbols) {
    btnCopySymbols.addEventListener("click", () => {
      const list = getOntologyGeneSymbols();
      if (!list.length) return;
      navigator.clipboard.writeText(list.join(" ")).then(() => {
        const orig = btnCopySymbols.textContent;
        btnCopySymbols.textContent = "✓ Symbols Copied!";
        setTimeout(() => { btnCopySymbols.textContent = orig; }, 2000);
      });
    });
  }

  const btnCopyEnsembl = document.getElementById("btnCopyEnsembl");
  if (btnCopyEnsembl) {
    btnCopyEnsembl.addEventListener("click", () => {
      const symbols = getOntologyGeneSymbols();
      const ids = symbols.map(s => {
        const found = MOCK_GENES.find(g => g.symbol.toUpperCase() === s.toUpperCase());
        return found ? found.id : s;
      });
      navigator.clipboard.writeText(ids.join("\n")).then(() => {
        const orig = btnCopyEnsembl.textContent;
        btnCopyEnsembl.textContent = "✓ Ensembl IDs Copied!";
        setTimeout(() => { btnCopyEnsembl.textContent = orig; }, 2000);
      });
    });
  }

  const btnDownloadGeneList = document.getElementById("btnDownloadGeneList");
  if (btnDownloadGeneList) {
    btnDownloadGeneList.addEventListener("click", () => {
      const list = getOntologyGeneSymbols();
      if (!list.length) return;
      downloadFile(list.join("\n"), "qgat_candidate_genes.txt", "text/plain");
    });
  }

  const btnResetCandidateGenes = document.getElementById("btnResetCandidateGenes");
  if (btnResetCandidateGenes) {
    btnResetCandidateGenes.addEventListener("click", () => {
      state.ontologyGenes = (MOCK_GENES || []).slice();
      const ta = document.getElementById("ontologyGenesTextarea");
      if (ta) ta.value = state.ontologyGenes.map(g => g.symbol).join(" ");
      updateOntologyView();
    });
  }

  const btnLaunchGProfiler = document.getElementById("btnLaunchGProfiler");
  if (btnLaunchGProfiler) btnLaunchGProfiler.addEventListener("click", launchGProfiler);

  const btnLaunchDAVID = document.getElementById("btnLaunchDAVID");
  if (btnLaunchDAVID) btnLaunchDAVID.addEventListener("click", launchDAVID);

  const btnLaunchPanther = document.getElementById("btnLaunchPanther");
  if (btnLaunchPanther) btnLaunchPanther.addEventListener("click", launchPanther);

  const btnLaunchShinyGO = document.getElementById("btnLaunchShinyGO");
  if (btnLaunchShinyGO) btnLaunchShinyGO.addEventListener("click", launchShinyGO);

  document.querySelectorAll("[data-copy-citation]").forEach(btn => {
    btn.addEventListener("click", () => {
      copyCitation(btn.dataset.copyCitation, btn);
    });
  });



  const btnGoToOntology = document.getElementById("btnGoToOntology");
  if (btnGoToOntology) {
    btnGoToOntology.addEventListener("click", () => {
      state.ontologyGenes = (MOCK_GENES || []).slice();
      const ta = document.getElementById("ontologyGenesTextarea");
      if (ta) ta.value = state.ontologyGenes.map(g => g.symbol).join(" ");
      navigateTo("ontology");
    });
  }

  // ---- Resize Handler ----
  window.addEventListener("resize", () => {
    if (state.currentPage === "results") {
      drawGenomeBrowser();
    }
  });
}

// ---- Utility ----
function downloadFile(content, filename, mimeType) {
  const blob = new Blob([content], { type: mimeType });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

// ---- Polyfill for roundRect ----
if (!CanvasRenderingContext2D.prototype.roundRect) {
  CanvasRenderingContext2D.prototype.roundRect = function(x, y, w, h, r) {
    if (typeof r === 'number') r = [r, r, r, r];
    this.moveTo(x + r[0], y);
    this.lineTo(x + w - r[1], y);
    this.quadraticCurveTo(x + w, y, x + w, y + r[1]);
    this.lineTo(x + w, y + h - r[2]);
    this.quadraticCurveTo(x + w, y + h, x + w - r[2], y + h);
    this.lineTo(x + r[3], y + h);
    this.quadraticCurveTo(x, y + h, x, y + h - r[3]);
    this.lineTo(x, y + r[0]);
    this.quadraticCurveTo(x, y, x + r[0], y);
    this.closePath();
  };
}

// ---- Init on DOM ready ----
document.addEventListener("DOMContentLoaded", init);
