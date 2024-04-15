/*!

Translate DNA or RNA sequence into protein.

"X" = invalid amino acid  
"*" = stop codon
You can find the rest of the table
[elsewhere](https://en.wikipedia.org/wiki/DNA_codon_table).

*/

#![deny(missing_docs)]

/// Translate DNA or RNA sequence into a peptide.
///  
/// # Examples
///
/// ```
/// use prot-translate::translate;
///
/// # fn main() { dna_example(); rna_example(); shift_reading_frame(); }
/// fn dna_example() {
///     
///     let dna = b"GCTAGTCGTATCGTAGCTAGTC";
///     let peptide = translate(dna,None);
///     assert_eq!(&peptide, "ASRIVAS");
///
/// }
///
/// fn rna_example() {
///
///     let rna = b"GCUAGUCGUAUCGUAGCUAGUC";
///     let peptide = translate(rna,None);
///     assert_eq!(&peptide, "ASRIVAS");
/// }
///
/// fn shift_reading_frame() {
///
///     // To shift the reading frame, pass in a slice
///     // skipping the first 1-2 nucleotides.
///
///     let dna = b"GCTAGTCGTATCGTAGCTAGTC";
///     let peptide_frame2 = translate(&dna[1..],None);
///     assert_eq!(&peptide_frame2, "LVVS*LV");
///
///     let peptide_frame3 = translate(&dna[2..],"-".chars().next());
///     assert_eq!(&peptide_frame3, "*-S-Y-R-S-*");
///    
/// }
/// ```
pub fn translate(seq: &[u8], delimiter: Option<char>) -> String {
    let mut peptide = String::with_capacity(seq.len().checked_div(3).unwrap());
    let delimit = delimiter.unwrap_or(char::from_u32(0).unwrap());
    'outer: for triplet in seq.chunks_exact(3) {
        for c in triplet {
            if !c.is_ascii() {
                peptide.push('X');
                continue 'outer;
            }
        }
         let c1 = ASCII_TO_INDEX[triplet[0] as usize];
        let c2 = ASCII_TO_INDEX[triplet[1] as usize];
        let c3 = ASCII_TO_INDEX[triplet[2] as usize];
        let amino_acid:char = if c1 == 4 || c2 == 4 || c3 == 4 {
            'X'
        } else {
            AA_TABLE_CANONICAL[c1][c2][c3]
        };

        peptide.push_str(&amino_acid.to_string());
        if delimit != '\0' { peptide.push_str(&delimit.to_string()) };
    }
    peptide.trim_end_matches(delimit).to_string()
}
/// Translate DNA or RNA sequence into a full-letter peptide.
///  
/// # Examples
///
/// ```
/// use prot-translate::translate3;
///
/// # fn main() { dna_example(); rna_example(); shift_reading_frame(); }
/// fn dna_example() {
///     
///     let dna = b"GCTAGTCGTATCGTAGCTAGTC";
///     let peptide = translate3(dna,None);
///     assert_eq!(&peptide, "AlaSerArgIleValAlaSer");
///
/// }
///
/// fn rna_example() {
///
///     let rna = b"GCUAGUCGUAUCGUAGCUAGUC";
///     let peptide = translate3(rna,None);
///     assert_eq!(&peptide, "AlaSerArgIleValAlaSer");
/// }
///
/// fn shift_reading_frame() {
///
///     // To shift the reading frame, pass in a slice
///     // skipping the first 1-2 nucleotides.
///
///     let dna = b"GCTAGTCGTATCGTAGCTAGTC";
///     let peptide_frame2 = translate3(&dna[1..],"-".chars().next());
///     assert_eq!(&peptide_frame2, "Leu-Val-Val-Ser-*-Leu-Val");
///
///     let peptide_frame3 = translate3(&dna[2..],"-".chars().next());
///     assert_eq!(&peptide_frame3, "*-Ser-Tyr-Arg-Ser-*");
///    
/// }
/// ```
pub fn translate3(seq: &[u8], delimiter: Option<char>) -> String {
    let mut peptide = String::with_capacity(seq.len());
    let delimit = delimiter.unwrap_or(char::from_u32(0).unwrap());
    'outer: for triplet in seq.chunks_exact(3) {
        for c in triplet {
            if !c.is_ascii() {
                peptide.push('X');
                continue 'outer;
            }
        }
         let c1 = ASCII_TO_INDEX[triplet[0] as usize];
        let c2 = ASCII_TO_INDEX[triplet[1] as usize];
        let c3 = ASCII_TO_INDEX[triplet[2] as usize];
        let amino_acid:&str = if c1 == 4 || c2 == 4 || c3 == 4 {
            "X"
        } else {
            AA_TABLE_CANONICAL_3[c1][c2][c3]
        };

        peptide.push_str(amino_acid);
        if delimit != '\0' { peptide.push_str(&delimit.to_string()) };
    }
    peptide.trim_end_matches(delimit).to_string()
}
/// Translate DNA or RNA sequence into a 3-letter peptide.
///  
/// # Examples
///
/// ```
/// use prot-translate::translate_full;
///
/// # fn main() { dna_example(); rna_example(); shift_reading_frame(); }
/// fn dna_example() {
///     
///     let dna = b"GCTAGTCGTATCGTAGCTAGTC";
///     let peptide = translate_full(dna,None);
///     assert_eq!(&peptide, "AlanineSerineArginineIsoleucineValineAlanineSerine");
///
/// }
///
/// fn rna_example() {
///
///     let rna = b"GCUAGUCGUAUCGUAGCUAGUC";
///     let peptide = translate_full(rna,None);
///     assert_eq!(&peptide, "AlanineSerineArginineIsoleucineValineAlanineSerine");
/// }
///
/// fn shift_reading_frame() {
///
///     // To shift the reading frame, pass in a slice
///     // skipping the first 1-2 nucleotides.
///
///     let dna = b"GCTAGTCGTATCGTAGCTAGTC";
///     let peptide_frame2 = translate_full(&dna[1..],"-".chars().next());
///     assert_eq!(&peptide_frame2, "Leucine-Valine-Valine-Serine-STOP-Leucine-Valine");
///
///     let peptide_frame3 = translate_full(&dna[2..],"-".chars().next());
///     assert_eq!(&peptide_frame3, "STOP-Serine-Tyrosine-Arginine-Serine-STOP");
///    
/// }
/// ```
pub fn translate_full(seq: &[u8], delimiter: Option<char>) -> String {
    let mut peptide = String::with_capacity(seq.len()*4); //12-length at most
    let delimit = delimiter.unwrap_or(char::from_u32(0).unwrap());
    'outer: for triplet in seq.chunks_exact(3) {
        for c in triplet {
            if !c.is_ascii() {
                peptide.push('X');
                continue 'outer;
            }
        }
         let c1 = ASCII_TO_INDEX[triplet[0] as usize];
        let c2 = ASCII_TO_INDEX[triplet[1] as usize];
        let c3 = ASCII_TO_INDEX[triplet[2] as usize];
        let amino_acid:&str = if c1 == 4 || c2 == 4 || c3 == 4 {
            "X"
        } else {
            AA_TABLE_CANONICAL_FULL[c1][c2][c3]
        };

        peptide.push_str(amino_acid);
        if delimit != '\0' { peptide.push_str(&delimit.to_string()) };
    }
    peptide.trim_end_matches(delimit).to_string()
}
/// https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
/// U is equivalent to T here
///
/// The 1st index picks the 4x4 block
/// The 2nd index picks the row
/// the 3rd index picks the column
static AA_TABLE_CANONICAL: [[[char; 4]; 4]; 4] = [
    [
        ['K', 'N', 'K', 'N'], // AAA, AAC, AAG, AAU/AAT
        ['T', 'T', 'T', 'T'], // ACA, ACC, ACG, ACU/ACT
        ['R', 'S', 'R', 'S'], // AGA, AGC, AGG, AGU/AGT
        ['I', 'I', 'M', 'I'], // AUA/ATA, AUC/ATC, AUG/ATG, AUU/ATT
    ],
    [
        ['Q', 'H', 'Q', 'H'], // CAA, CAC, CAG, CAU/CAT
        ['P', 'P', 'P', 'P'], // CCA, CCC, CCG, CCU/CCT
        ['R', 'R', 'R', 'R'], // CGA, CGC, CGG, CGU/CGT
        ['L', 'L', 'L', 'L'], // CUA/CTA, CUC/CTC, CUG/CTG, CUU/CTT
    ],
    [
        ['E', 'D', 'E', 'D'], // GAA, GAC, GAG, GAU/GAT
        ['A', 'A', 'A', 'A'], // GCA, GCC, GCG, GCU/GCT
        ['G', 'G', 'G', 'G'], // GGA, GGC, GGG, GGU/GGT
        ['V', 'V', 'V', 'V'], // GUA/GTA, GUC/GTC, GUG/GTG, GUU/GTT
    ],
    [
        ['*', 'Y', '*', 'Y'], // UAA/TAA, UAC/TAC, UAG/TAG, UAU/TAT
        ['S', 'S', 'S', 'S'], // UCA/TCA, UCC/TCC, UCG/TCG, UCU/TCT
        ['*', 'C', 'W', 'C'], // UGA/TGA, UGC/TGC, UGG/TGG, UGU/TGT
        ['L', 'F', 'L', 'F'], // UUA/TTA, UUC/TTC, UUG/TTG, UUU/TTT
    ],
];
static AA_TABLE_CANONICAL_3: [[[&str; 4]; 4]; 4] = [
    [
        ["Lys", "Asn", "Lys", "Asn"], // AAA, AAC, AAG, AAU/AAT
        ["Thr", "Thr", "Thr", "Thr"], // ACA, ACC, ACG, ACU/ACT
        ["Arg", "Ser", "Arg", "Ser"], // AGA, AGC, AGG, AGU/AGT
        ["Ile", "Ile", "Met", "Ile"], // AUA/ATA, AUC/ATC, AUG/ATG, AUU/ATT
    ],
    [
        ["Gln", "His", "Gln", "His"], // CAA, CAC, CAGly, CAU/CAT
        ["Pro", "Pro", "Pro", "Pro"], // CCA, CCC, CCGly, CCU/CCT
        ["Arg", "Arg", "Arg", "Arg"], // CGlyA, CGlyC, CGlyGly, CGlyU/CGlyT
        ["Leu", "Leu", "Leu", "Leu"], // CUA/CTA, CUC/CTC, CUGly/CTGly, CUU/CTT
    ],
    [
        ["Glu", "Asp", "Glu", "Asp"], // GlyAA, GlyAC, GlyAGly, GlyAU/GlyAT
        ["Ala", "Ala", "Ala", "Ala"], // GlyCA, GlyCC, GlyCGly, GlyCU/GlyCT
        ["Gly", "Gly", "Gly", "Gly"], // GlyGlyA, GlyGlyC, GlyGlyGly, GlyGlyU/GlyGlyT
        ["Val", "Val", "Val", "Val"], // GlyUA/GlyTA, GlyUC/GlyTC, GlyUGly/GlyTGly, GlyUU/GlyTT
    ],
    [
        ["*", "Tyr", "*", "Tyr"], // UAA/TAA, UAC/TAC, UAGly/TAGly, UAU/TAT
        ["Ser", "Ser", "Ser", "Ser"], // UCA/TCA, UCC/TCC, UCGly/TCGly, UCU/TCT
        ["*", "Cys", "Trp", "Cys"], // UGlyA/TGlyA, UGlyC/TGlyC, UGlyGly/TGlyGly, UGlyU/TGlyT
        ["Leu", "Phe", "Leu", "Phe"], // UUA/TTA, UUC/TTC, UUGly/TTGly, UUU/TTT
    ],
];
static AA_TABLE_CANONICAL_FULL: [[[&str; 4]; 4]; 4] = [
    [
        ["Lysine", "Asparagine", "Lysine", "Asparagine"], // AAA, AAC, AAGly, AAU/AAT
        ["Threonine", "Threonine", "Threonine", "Threonine"], // ACA, ACC, ACGly, ACU/ACT
        ["Arginine", "Serine", "Arginine", "Serine"], // AGlyA, AGlyC, AGlyGly, AGlyU/AGlyT
        ["Isoleucine", "Isoleucine", "Methionine", "Isoleucine"], // AUA/ATA, AUC/ATC, AUGly/ATGly, AUU/ATT
    ],
    [
        ["Glutamine", "Histidine", "Glutamine", "His"], // CAA, CAC, CAGly, CAU/CAT
        ["Proline", "Proline", "Proline", "Proline"], // CCA, CCC, CCGly, CCU/CCT
        ["Arginine", "Arginine", "Arginine", "Arginine"], // CGlyA, CGlyC, CGlyGly, CGlyU/CGlyT
        ["Leucine", "Leucine", "Leucine", "Leucine"], // CUA/CTA, CUC/CTC, CUGly/CTGly, CUU/CTT
    ],
    [
        ["Glutamine", "Aspartic acid", "Glu", "Asp"], // GlyAA, GlyAC, GlyAGly, GlyAU/GlyAT
        ["Alanine", "Alanine", "Alanine", "Alanine"], // GlyCA, GlyCC, GlyCGly, GlyCU/GlyCT
        ["Glycine", "Glycine", "Glycine", "Glycine"], // GlyGlyA, GlyGlyC, GlyGlyGly, GlyGlyU/GlyGlyT
        ["Valine", "Valine", "Valine", "Valine"], // GlyUA/GlyTA, GlyUC/GlyTC, GlyUGly/GlyTGly, GlyUU/GlyTT
    ],
    [
        ["STOP", "Tyrosine", "STOP", "Tyrosine"], // UAA/TAA, UAC/TAC, UAGly/TAGly, UAU/TAT
        ["Serine", "Serine", "Serine", "Serine"], // UCA/TCA, UCC/TCC, UCGly/TCGly, UCU/TCT
        ["", "Cysteine", "Tryptophan", "Cysteine"], // UGlyA/TGlyA, UGlyC/TGlyC, UGlyGly/TGlyGly, UGlyU/TGlyT
        ["Leucine", "Phenylalanine", "Leucine", "Phenylalanine"], // UUA/TTA, UUC/TTC, UUGly/TTGly, UUU/TTT
    ],
];
/// Maps an ASCII character to array index
///
/// A = 65, a = 97  => 0
/// C = 67, c = 99  => 1
/// G = 71, g = 103 => 2
/// T = 84, t = 116 => 3
/// U = 85, u = 117 => 3
static ASCII_TO_INDEX: [usize; 128] = [
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 0-15
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 16-31
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 32-47
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 48-63
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, // 64-79    (65 = A, 67 = C, 71 = G)
    4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 80-95    (84 = T, 85 = U)
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, // 96-111   (97 = a, 99 = c, 103 = g)
    4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 112-127  (116 = t, 117 = u)
];
