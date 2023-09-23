//! # nounoursmap
//! Trouver le mangeur de choconours en comparant son génome
//! avec de l'ADN issu de la salive laissée sur les lieux du méfait.

use bio::io::fasta;
use clap::Parser;

use crate::pseudomap::u8_to_string;

mod pseudomap;

#[derive(Parser)]
#[clap(version = "0.1.0", author = "Samuel Ortion")]
struct Cli {
    #[clap(short = 'g', long = "genomes", help = "Fichier FASTA d'entrée avec les génomes")]
    fasta_genomes: String,
    #[clap(short = 'r', long = "reads", help = "Fichier FASTA d'entrée avec les reads")]
    fasta_reads: String
}

fn main() {
    let mut nom_coupables: Vec<String> = Vec::new();
    let mut compte_read_coupable: i32 = 0;
    // Récupérer les valeurs des paramètres du programme
    let args = Cli::parse();
    let fasta_genomes = args.fasta_genomes;
    let fasta_reads = args.fasta_reads;
    let genomes_reader = fasta::Reader::from_file(fasta_genomes)
        .expect("Une erreur est survenue lors de la lecture du fichier FASTA des génomes.");
    let reads_reader = fasta::Reader::from_file(fasta_reads)
        .expect("Une erreur est survenue lors de la lecture du fichier FASTA des reads.");
    // Récupération des reads dans un tableau
    let mut reads: Vec<Vec<u8>> = Vec::new();
    for read_result in reads_reader.records() {
        let read_record = read_result.expect("Une erreur est survenue lors de l'analyse du fichier FASTA des reads");
        let read_sequence = read_record.seq().to_vec();
        reads.push(read_sequence);
    }
    for genome_result in genomes_reader.records() {
        let genome_record = genome_result
            .expect("Une erreur est survenue lors de l'analyse du fichier FASTA des génomes");
        let suspect = genome_record.id().to_owned();
        // println!("Le génome de {} va être analyser...", suspect);
        let genome = genome_record.seq();
        let mut total_count: i32 = 0;
        println!("{}", pseudomap::u8_to_string(genome));
        for read in reads.iter() {
            let read_seq = &read;
            let read_count: i32 = pseudomap::compter_read(genome, read_seq);
            let read_positions: Vec<usize> = pseudomap::trouver_positions(genome, read_seq);
            
            println!("{:?}, {:?}", pseudomap::u8_to_string(read), read_positions);
            total_count += read_count;
        }
        println!("{}\t{:?} reads", suspect, total_count);
        // Mettre à jour la liste des coupable les plus probables
        // Si on trouve un suspect avec plus de reads 
        // que ceux de la liste précédente, 
        // on enlève les présumés coupables 
        // et on ajoute le nouveau coupable potentiel.
        if total_count > compte_read_coupable {
            compte_read_coupable = total_count;
            nom_coupables = vec![suspect];
        // Si on trouve un suspect avec le même nombre
        // de reads que le précédent suspect avec le plus de read,
        // on l'ajoute à la liste des coupables potentiels.
        } else if total_count == compte_read_coupable {
            nom_coupables.push(suspect);
        }
    }
    println!("Les suspect les plus louches ont {:?} reads de l'ADN de la salive détectés dans leur génome.", compte_read_coupable);
    // Afficher la liste des suspect louches
    println!("Il s'agit des suspects suivants");
    for nom in nom_coupables {
        println!(" - {}", nom);
    }
}
