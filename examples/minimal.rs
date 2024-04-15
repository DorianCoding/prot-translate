use prot-translate::translate;

fn main() {
    let dna = b"GTGAGTCGTTGAGTCTGATTGCGTATC";
    
    let protein = translate(dna,None);
    assert_eq!("VSR*V*LRI", &protein);

    // To shift reading frame
    let protein_frame2 = translate(&dna[1..],None);
    assert_eq!("*VVESDCV", &protein_frame2);

}