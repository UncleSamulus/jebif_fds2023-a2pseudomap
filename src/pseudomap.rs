//! # Pseudo Map
//! Trouver la position d'un read dans un génome, 
//! en supposant qu'il n'y pas de modification
//! dans la chaîne de caractères du read 
//! par rapport au génome.
//! Il s'agit de faire une recherche de sous-séquence.

use std::str;

/// Simple compte du nombre de reads présents dans le génome
pub fn compter_read(genome: &[u8], read: &[u8]) -> i32 {
    let mut count = 0; // On initialise le compteur à 0
    // Pour chaque position du génome, 
    // on extrait une partie de la séquence
    // et on la compare au read
    for i in 0..(genome.len() - read.len()) {
        let sub_sequence = &genome[i..(i+read.len())];
        if eq_simpler(sub_sequence, read) {
            count += 1; // On ajoute 1 au compteur
        }
    }
    return count; // On renvoit le compteur
}

// Comparer deux séquences de même taille
fn eq(seq1: &[u8], seq2: &[u8]) -> bool {
    seq1.iter().zip(seq2.iter()).all(|(a,b)| a == b) 
}

// On pourrait l'écrire aussi comme suit, par exemple
fn eq_simpler(seq1: &[u8], seq2: &[u8]) -> bool {
    // Si les séquences ne sont pas de même taille: elle ne sont pas égale, on retourne Faux.
    if seq1.len() != seq2.len() {
        return false;
    } // Sinon
    for i in 0..seq1.len() {
        // Si le nucleotide à la position i n'est pas le même entre les deux séquences, 
        // les séquences ne sont pas identiques
        if seq1[i] != seq2[i] {
            return false;
        }
    }
    // Si tous les nucléotides sont identiques, on retourne Vrai.
    return true;
}

/// Similaire à `compter_read`, mais cette fois 
/// on garde en mémoire la position des reads qui 
/// sont retrouvés dans le génome
pub fn trouver_positions(genome: &[u8], read: &[u8]) -> Vec<usize> {
    let mut positions: Vec<usize> = Vec::new();
    for i in 0..(genome.len() - read.len()) {
        let sub_sequence = &genome[i..(i+read.len())];
        if eq_simpler(sub_sequence, read) {
            positions.push(i);
        }
    }
    positions
}


pub fn u8_to_string(read: &[u8]) -> &str {
    let s = match str::from_utf8(read) {
        Ok(v) => v,
        Err(e) => panic!("Invalid UTF-8 sequence: {}", e),
    };
    s
}

// Définition de quelques tests des fonctions définies ci-dessus
// Pour s'assurer qu'elle fonctionne comme attendu.
#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_compter_read() {
        let genome: [u8; 10] = [b'A', b'T', b'G', b'C', b'T', b'G', b'C', b'T', b'A', b'T']; // On déclare une séquence de nucleotide sous forme de tableau de nombre (b'A' représente l'encodage ascii du A par un numbre)
        let read: [u8; 2] = [b'G', b'C']; // On définit pareillement une petite séquence lue
        assert_eq!(compter_read(&genome, &read), 2);
    }
}