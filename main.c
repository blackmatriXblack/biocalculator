#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h> // For toupper(), tolower(), isspace()
#include <math.h>  // For log(), fabs(), pow()

// --- Constants ---

// Molecular Weights (Average Daltons)
// Nucleotide Residue Weights (incorporated into chain)
#define MW_A_RESIDUE_DNA 313.21
#define MW_T_RESIDUE_DNA 304.19
#define MW_C_RESIDUE_DNA 289.18
#define MW_G_RESIDUE_DNA 329.21
#define MW_A_RESIDUE_RNA 329.21
#define MW_U_RESIDUE_RNA 306.17
#define MW_C_RESIDUE_RNA 305.18
#define MW_G_RESIDUE_RNA 345.21
#define MW_H2O           18.015

// Amino Acid Residue Weights (minus H2O)
#define MW_ALA 71.0788
#define MW_ARG 156.1875
#define MW_ASN 114.1038
#define MW_ASP 115.0886
#define MW_CYS 103.1388
#define MW_GLN 128.1292
#define MW_GLU 129.1155
#define MW_GLY 57.0519
#define MW_HIS 137.1411
#define MW_ILE 113.1594
#define MW_LEU 113.1594
#define MW_LYS 128.1492
#define MW_MET 131.1926
#define MW_PHE 147.1766
#define MW_PRO 97.1167
#define MW_SER 87.0782
#define MW_THR 101.1051
#define MW_TRP 186.2139
#define MW_TYR 163.1760
#define MW_VAL 99.1357

// Nearest Neighbor Tm Parameters (dH in kcal/mol, dS in cal/mol/K)
// Source: SantaLucia, J. (1998). PNAS, 95(4), 1460-1465.
//       dH     dS
double NN_PARAMS[10][2] = {
    {-7.9, -22.2}, // AA/TT
    {-8.4, -23.6}, // AT/TA (AT)
    {-7.8, -21.9}, // TA/AT (TA)
    {-10.2, -26.2}, // CG/GC
    {-10.6, -27.2}, // GC/CG
    {-8.0, -19.0}, // GG/CC
    {-8.0, -19.9}, // CC/GG - Note: some tables combine GG/CC, others separate. Using slightly different values for CC/GG here.
    {-8.5, -22.7}, // AC/GT
    {-8.4, -22.4}, // GT/AC
    {-7.8, -21.0}  // CA/TG
};
// Index Mapping: 0:AA/TT, 1:AT, 2:TA, 3:CG, 4:GC, 5:GG/CC (usually GG), 6:CC/GG (usually CC), 7:AC/GT, 8:GT/AC, 9:CA/TG

// Gas constant (cal/mol/K)
#define R_GAS 1.987

// Default strand concentration for Tm calculation (Molar)
#define DEFAULT_DNA_CONC 50e-9 // 50 nM

// Standard pKa values for protein pI calculation
// Source: IPC - Isoelectric Point Calculator (Bioinformatics, 2004) - using average/common values
#define PK_N_TERM   8.0
#define PK_C_TERM   3.1
#define PK_ASP      4.1
#define PK_GLU      4.5
#define PK_CYS      6.8 // Note: pKa can vary significantly with environment
#define PK_TYR     10.5
#define PK_HIS      6.0
#define PK_LYS     10.5
#define PK_ARG     12.5

// pI calculation parameters
#define PH_STEP     0.01 // Step size for pH iteration
#define MAX_PH     14.0
#define MIN_PH      0.0

// Molar Extinction Coefficients (Single Strand, at 260nm)
// Source: Based on values from biochemistry texts/online resources (units: M^-1 cm^-1)
#define EPSILON_A_NUCL 15400.0
#define EPSILON_T_NUCL 8700.0
#define EPSILON_C_NUCL 7300.0
#define EPSILON_G_NUCL 11500.0
#define EPSILON_U_NUCL 9900.0

// Molar Extinction Coefficients (Protein, at 280nm in water)
// Source: Pace et al. (1995) Protein Science 4:2411-2423 (units: M^-1 cm^-1)
#define EPSILON_TRP_PROT 5690.0
#define EPSILON_TYR_PROT 1280.0
#define EPSILON_CYS_PROT_REDUCED 120.0 // For reduced cysteine (disulfides have negligible absorbance at 280nm)

// Hydrophobicity values (Kyte & Doolittle Scale)
// Source: Kyte, J., & Doolittle, R. F. (1982). J. Mol. Biol., 157(1), 105-132.
// Order: ARNDCEQGHILKMFPSTWYV
double HYDROPHOBICITY_SCALE[20] = {
    -0.8, -4.5, -3.5, -3.5, 2.5, -3.5, -3.5, -0.4, -3.2, 4.5, 3.8, -3.9, 1.9, 2.8, -1.6, -0.8, -0.7, -0.9, -1.3, 4.2
};

// Standard Genetic Code (RNA codons to Amino Acids)
// This is a lookup table for translation.
// Index mapping: First base (U,C,A,G) -> Second base (U,C,A,G) -> Third base (U,C,A,G)
// Example: CODE[0][0][0] is UUU -> F (Phenylalanine)
char GENETIC_CODE[4][4][4] = {
    {   // First base: U
        {'F', 'F', 'L', 'L'}, // Second base: U (UUU, UUC, UUA, UUG)
        {'S', 'S', 'S', 'S'}, // Second base: C (UCU, UCC, UCA, UCG)
        {'Y', 'Y', '_', '_'}, // Second base: A (UAU, UAC, UAA (Stop), UAG (Stop))
        {'C', 'C', '_', 'W'}  // Second base: G (UGU, UGC, UGA (Stop), UGG)
    },
    {   // First base: C
        {'L', 'L', 'L', 'L'}, // Second base: U
        {'P', 'P', 'P', 'P'}, // Second base: C
        {'H', 'H', 'Q', 'Q'}, // Second base: A
        {'R', 'R', 'R', 'R'}  // Second base: G
    },
    {   // First base: A
        {'I', 'I', 'I', 'M'}, // Second base: U (AUA, AUC, AUU, AUG (Start/Met))
        {'T', 'T', 'T', 'T'}, // Second base: C
        {'N', 'N', 'K', 'K'}, // Second base: A
        {'S', 'S', 'R', 'R'}  // Second base: G
    },
    {   // First base: G
        {'V', 'V', 'V', 'V'}, // Second base: U
        {'A', 'A', 'A', 'A'}, // Second base: C
        {'D', 'D', 'E', 'E'}, // Second base: A
        {'G', 'G', 'G', 'G'}  // Second base: G
    }
};

// Mapping for base character to genetic code index (U->0, C->1, A->2, G->3)
int base_to_index(char base) {
    switch (base) {
        case 'U': return 0;
        case 'C': return 1;
        case 'A': return 2;
        case 'G': return 3;
        default: return -1; // Should not happen with validated RNA
    }
}


// --- Helper Functions ---

// Function to convert sequence to uppercase and remove whitespace/newlines
void clean_sequence(char *seq) {
    if (!seq) return;
    int i, j = 0;
    for (i = 0; seq[i] != '\0'; i++) {
        if (!isspace((unsigned char)seq[i]) && seq[i] != '\n' && seq[i] != '\r') {
            seq[j++] = toupper((unsigned char)seq[i]);
        }
    }
    seq[j] = '\0';
}

// Validate DNA sequence characters
int is_valid_dna(const char *seq) {
    if (!seq || seq[0] == '\0') return 0;
    for (int i = 0; seq[i] != '\0'; i++) {
        if (seq[i] != 'A' && seq[i] != 'T' && seq[i] != 'C' && seq[i] != 'G') {
            return 0; // Invalid character found
        }
    }
    return 1; // All characters are valid
}

// Validate RNA sequence characters
int is_valid_rna(const char *seq) {
    if (!seq || seq[0] == '\0') return 0;
    for (int i = 0; seq[i] != '\0'; i++) {
        if (seq[i] != 'A' && seq[i] != 'U' && seq[i] != 'C' && seq[i] != 'G') {
            return 0; // Invalid character found
        }
    }
    return 1; // All characters are valid
}

// Validate Protein sequence characters (20 standard amino acids)
int is_valid_protein(const char *seq) {
    if (!seq || seq[0] == '\0') return 0;
    const char *valid_chars = "ARNDCEQGHILKMFPSTWYV";
    for (int i = 0; seq[i] != '\0'; i++) {
        if (strchr(valid_chars, seq[i]) == NULL) {
            return 0; // Invalid character found
        }
    }
    return 1; // All characters are valid
}

// Get NN parameter index for a dinucleotide (e.g., "AT")
int get_nn_index(char base1, char base2) {
    if (base1 == 'A' && base2 == 'A') return 0;
    if (base1 == 'T' && base2 == 'T') return 0;

    if (base1 == 'A' && base2 == 'T') return 1;
    if (base1 == 'T' && base2 == 'A') return 2;

    if (base1 == 'C' && base2 == 'G') return 3;
    if (base1 == 'G' && base2 == 'C') return 4;

    if (base1 == 'G' && base2 == 'G') return 5;
    if (base1 == 'C' && base2 == 'C') return 6; // Using separate indices for CC/GG

    if (base1 == 'A' && base2 == 'C') return 7;
    if (base1 == 'G' && base2 == 'T') return 7; // AC/GT are reverse complements

    if (base1 == 'G' && base2 == 'T') return 8; // GT/AC
    if (base1 == 'A' && base2 == 'C') return 8; // GT/AC

    if (base1 == 'C' && base2 == 'A') return 9;
    if (base1 == 'T' && base2 == 'G') return 9;

    return -1; // Should not happen with validated DNA
}


// --- Calculation Functions ---

// Calculate sequence length
size_t calculate_length(const char *seq) {
    return strlen(seq);
}

// Calculate GC content for DNA or RNA
double calculate_gc_content(const char *seq) {
    size_t len = strlen(seq);
    if (len == 0) return 0.0;
    int gc_count = 0;
    for (int i = 0; seq[i] != '\0'; i++) {
        if (seq[i] == 'G' || seq[i] == 'C') {
            gc_count++;
        }
    }
    return ((double)gc_count / len) * 100.0;
}

// Calculate Molecular Weight for DNA or RNA
double calculate_nucleic_mw(const char *seq, const char *type) {
    size_t len = strlen(seq);
    if (len == 0) return 0.0;

    double total_mw = 0.0;

    if (strcasecmp(type, "DNA") == 0) {
        for (int i = 0; seq[i] != '\0'; i++) {
            switch (seq[i]) {
                case 'A': total_mw += MW_A_RESIDUE_DNA; break;
                case 'T': total_mw += MW_T_RESIDUE_DNA; break;
                case 'C': total_mw += MW_C_RESIDUE_DNA; break;
                case 'G': total_mw += MW_G_RESIDUE_DNA; break;
            }
        }
    } else if (strcasecmp(type, "RNA") == 0) {
         for (int i = 0; seq[i] != '\0'; i++) {
            switch (seq[i]) {
                case 'A': total_mw += MW_A_RESIDUE_RNA; break;
                case 'U': total_mw += MW_U_RESIDUE_RNA; break;
                case 'C': total_mw += MW_C_RESIDUE_RNA; break;
                case 'G': total_mw += MW_G_RESIDUE_RNA; break;
            }
        }
    }
    // Add the mass of water molecule for the terminals (common for single linear strand)
    total_mw += MW_H2O;

    return total_mw;
}


// Calculate Molecular Weight for Protein
double calculate_protein_mw(const char *seq) {
    size_t len = strlen(seq);
    if (len == 0) return 0.0;

    double total_mw = 0.0; // Fix: Changed total_tbl to total_mw

    for (int i = 0; seq[i] != '\0'; i++) {
        switch (seq[i]) {
            case 'A': total_mw += MW_ALA; break;
            case 'R': total_mw += MW_ARG; break;
            case 'N': total_mw += MW_ASN; break;
            case 'D': total_mw += MW_ASP; break;
            case 'C': total_mw += MW_CYS; break;
            case 'Q': total_mw += MW_GLN; break;
            case 'E': total_mw += MW_GLU; break;
            case 'G': total_mw += MW_GLY; break;
            case 'H': total_mw += MW_HIS; break;
            case 'I': total_mw += MW_ILE; break;
            case 'L': total_mw += MW_LEU; break;
            case 'K': total_mw += MW_LYS; break;
            case 'M': total_mw += MW_MET; break;
            case 'F': total_mw += MW_PHE; break;
            case 'P': total_mw += MW_PRO; break;
            case 'S': total_mw += MW_SER; break;
            case 'T': total_mw += MW_THR; break; // Fix: Changed total_tbl to total_mw
            case 'W': total_mw += MW_TRP; break;
            case 'Y': total_mw += MW_TYR; break;
            case 'V': total_mw += MW_VAL; break;
            // Validation should prevent other characters
        }
    }
    // Add the mass of water molecule for the terminals (common for linear peptide)
    total_mw += MW_H2O;

    return total_mw;
}

// Calculate Melting Temperature (Tm) for DNA using Nearest Neighbor method
// Formula: Tm = (dH * 1000) / (dS + R * log(Ct / 4)) - 273.15
// dH in kcal/mol, dS in cal/mol/K, R in cal/mol/K, Ct in Molar
double calculate_dna_tm(const char *seq, double concentration) {
    size_t len = strlen(seq);
    if (len < 2) {
         if (len == 1) {
             printf("Error: Sequence too short for Nearest Neighbor Tm calculation (requires >= 2 bases).\n");
             return -1.0; // Indicate error
         }
         return 0.0; // Empty sequence
    }

    double total_dH = 0.0; // kcal/mol
    double total_dS = 0.0; // cal/mol/K

    // Sum parameters for internal dinucleotides
    for (size_t i = 0; i < len - 1; i++) {
        char base1 = seq[i];
        char base2 = seq[i+1];
        int index = get_nn_index(base1, base2);
        if (index != -1) {
            total_dH += NN_PARAMS[index][0];
            total_dS += NN_PARAMS[index][1];
        } else {
            // This case should ideally not be reached if sequence validation is perfect
             printf("Warning: Could not get NN index for dinucleotide %c%c. Skipping.\n", base1, base2);
        }
    }

    // Note: This simplified NN calculation does not include terminal mismatch, loop, or salt corrections beyond the base formula.

    // Calculate Tm in Kelvin
    // Handle case where concentration/4.0 might be non-positive or very small for log
    if (concentration <= 0) {
         printf("Error: DNA concentration must be positive for Tm calculation.\n");
         return -1.0;
    }
    double log_conc_factor = log(concentration / 4.0);

    double tm_k = (total_dH * 1000.0) / (total_dS + R_GAS * log_conc_factor);

    // Convert to Celsius
    double tm_c = tm_k - 273.15;

    return tm_c;
}

// Calculate Reverse Complement for DNA
char* calculate_dna_revcomp(const char *seq) {
    size_t len = strlen(seq);
    if (len == 0) return strdup("");

    char *revcomp_seq = (char *)malloc(len + 1);
    if (!revcomp_seq) {
        perror("Failed to allocate memory for reverse complement");
        return NULL;
    }

    int i, j;
    for (i = len - 1, j = 0; i >= 0; i--, j++) {
        switch (seq[i]) {
            case 'A': revcomp_seq[j] = 'T'; break;
            case 'T': revcomp_seq[j] = 'A'; break;
            case 'C': revcomp_seq[j] = 'G'; break;
            case 'G': revcomp_seq[j] = 'C'; break;
        }
    }
    revcomp_seq[len] = '\0';

    return revcomp_seq; // Remember to free this memory after use!
}

// Transcribe DNA to RNA
char* transcribe_dna_to_rna(const char *seq) {
    size_t len = strlen(seq);
    if (len == 0) return strdup("");

    char *rna_seq = (char *)malloc(len + 1);
    if (!rna_seq) {
        perror("Failed to allocate memory for RNA sequence");
        return NULL;
    }

    for (size_t i = 0; i < len; i++) {
        if (seq[i] == 'T') {
            rna_seq[i] = 'U';
        } else {
            rna_seq[i] = seq[i]; // A, C, G remain the same
        }
    }
    rna_seq[len] = '\0';

    return rna_seq; // Remember to free this memory after use!
}

// Translate RNA to Protein
char* translate_rna_to_protein(const char *rna_seq) {
    size_t rna_len = strlen(rna_seq);
    if (rna_len == 0) return strdup("");

    // Check if length is a multiple of 3
    if (rna_len % 3 != 0) {
        fprintf(stderr, "Warning: RNA sequence length (%zu) is not a multiple of 3. Trailing bases will be ignored.\n", rna_len);
    }

    // Maximum possible protein length (rna_len / 3) + 1 for null terminator
    size_t protein_len = rna_len / 3;
    char *protein_seq = (char *)malloc(protein_len + 1);
     if (!protein_seq) {
        perror("Failed to allocate memory for protein sequence");
        return NULL;
    }

    int protein_idx = 0;
    for (size_t i = 0; i + 2 < rna_len; i += 3) {
        char base1 = rna_seq[i];
        char base2 = rna_seq[i+1];
        char base3 = rna_seq[i+2];

        int idx1 = base_to_index(base1);
        int idx2 = base_to_index(base2);
        int idx3 = base_to_index(base3);

        // Check if bases are valid RNA characters (validation should catch this earlier, but defensive)
        if (idx1 == -1 || idx2 == -1 || idx3 == -1) {
             fprintf(stderr, "Error during translation: Invalid RNA base found.\n");
             free(protein_seq);
             return NULL;
        }

        char aa = GENETIC_CODE[idx1][idx2][idx3];

        if (aa == '_') { // Stop codon
            protein_seq[protein_idx] = '\0'; // Terminate the protein sequence
            break; // Stop translation
        } else if (aa == 'M' && protein_idx == 0) {
             // Start codon (AUG) at the beginning. Keep the 'M'.
             protein_seq[protein_idx++] = aa;
        } else {
            protein_seq[protein_idx++] = aa;
        }
    }

    protein_seq[protein_idx] = '\0'; // Ensure null termination

    return protein_seq; // Remember to free this memory after use!
}


// Calculate Amino Acid Composition for Protein
void calculate_aa_composition(const char *seq, int *composition_counts) {
    // Initialize counts to zero
    for (int i = 0; i < 20; i++) {
        composition_counts[i] = 0;
    }

    const char *aa_order = "ARNDCEQGHILKMFPSTWYV"; // Standard order for output
    size_t len = strlen(seq);

    for (size_t i = 0; i < len; i++) {
        char aa = seq[i];
        for (int j = 0; j < 20; j++) {
            if (aa == aa_order[j]) {
                composition_counts[j]++;
                break; // Found the amino acid, move to the next
            }
        }
    }
}


// Calculate Isoelectric Point (pI) for Protein
// Iterative method: calculate charge at different pH values and find where it's closest to zero.
double calculate_protein_pi(const char *seq) {
    size_t len = strlen(seq);
    if (len == 0) return -1.0; // Indicate error or undefined pI

    // Count ionizable residues
    int count_D = 0; // Aspartic Acid
    int count_E = 0; // Glutamic Acid
    int count_C = 0; // Cysteine
    int count_Y = 0; // Tyrosine
    int count_H = 0; // Histidine
    int count_K = 0; // Lysine
    int count_R = 0; // Arginine

    for (size_t i = 0; i < len; i++) {
        switch (seq[i]) {
            case 'D': count_D++; break;
            case 'E': count_E++; break;
            case 'C': count_C++; break;
            case 'Y': count_Y++; break;
            case 'H': count_H++; break;
            case 'K': count_K++; break;
            case 'R': count_R++; break;
        }
    }

    double best_ph = MIN_PH;
    double min_abs_charge = 1e10; // Initialize with a very large value

    // Iterate through pH values
    for (double ph = MIN_PH; ph <= MAX_PH; ph += PH_STEP) {
        double net_charge = 0.0;

        // N-terminus (positive charge when protonated)
        // pKa is the pH at which 50% is protonated.
        // Charge = +1 * (fraction protonated) = +1 * [H+] / ([H+] + Ka)
        // Charge = +1 / (1 + Ka/[H+]) = +1 / (1 + 10^(pH - pKa))
        net_charge += 1.0 / (1.0 + pow(10.0, ph - PK_N_TERM));

        // C-terminus (negative charge when deprotonated)
        // Charge = -1 * (fraction deprotonated) = -1 * Ka / ([H+] + Ka)
        // Charge = -1 / ([H+]/Ka + 1) = -1 / (10^(pKa - pH) + 1)
        net_charge -= 1.0 / (1.0 + pow(10.0, PK_C_TERM - ph));

        // Acidic side chains (negative charge when deprotonated)
        net_charge -= (double)count_D / (1.0 + pow(10.0, PK_ASP - ph));
        net_charge -= (double)count_E / (1.0 + pow(10.0, PK_GLU - ph));
        net_charge -= (double)count_C / (1.0 + pow(10.0, PK_CYS - ph));
        net_charge -= (double)count_Y / (1.0 + pow(10.0, PK_TYR - ph));

        // Basic side chains (positive charge when protonated)
        net_charge += (double)count_H / (1.0 + pow(10.0, ph - PK_HIS));
        net_charge += (double)count_K / (1.0 + pow(10.0, ph - PK_LYS));
        net_charge += (double)count_R / (1.0 + pow(1.0 + pow(10.0, ph - PK_ARG), 1.0)); // Pow to 1.0 is redundant, but keeps formula structure

        // Find the pH where the absolute charge is minimal
        if (fabs(net_charge) < min_abs_charge) {
            min_abs_charge = fabs(net_charge);
            best_ph = ph;
        }
    }

    return best_ph;
}

// Calculate Molar Extinction Coefficient for Nucleic Acid (at 260nm, single strand)
double calculate_nucleic_extinction(const char *seq, const char *type) {
     size_t len = strlen(seq);
    if (len == 0) return 0.0;

    double extinction_coeff = 0.0;

    if (strcasecmp(type, "DNA") == 0) {
        for (int i = 0; seq[i] != '\0'; i++) {
            switch (seq[i]) {
                case 'A': extinction_coeff += EPSILON_A_NUCL; break;
                case 'T': extinction_coeff += EPSILON_T_NUCL; break;
                case 'C': extinction_coeff += EPSILON_C_NUCL; break;
                case 'G': extinction_coeff += EPSILON_G_NUCL; break;
            }
        }
    } else if (strcasecmp(type, "RNA") == 0) {
         for (int i = 0; seq[i] != '\0'; i++) {
            switch (seq[i]) {
                case 'A': extinction_coeff += EPSILON_A_NUCL; break; // A is same in DNA/RNA at 260nm (approx)
                case 'U': extinction_coeff += EPSILON_U_NUCL; break;
                case 'C': extinction_coeff += EPSILON_C_NUCL; break; // C is same in DNA/RNA at 260nm (approx)
                case 'G': extinction_coeff += EPSILON_G_NUCL; break; // G is same in DNA/RNA at 260nm (approx)
            }
        }
    }

    return extinction_coeff;
}

// Calculate Molar Extinction Coefficient for Protein (at 280nm)
double calculate_protein_extinction(const char *seq) {
    size_t len = strlen(seq);
    if (len == 0) return 0.0;

    int count_W = 0; // Tryptophan
    int count_Y = 0; // Tyrosine
    int count_C = 0; // Cysteine (assuming reduced form)

    for (size_t i = 0; i < len; i++) {
        switch (seq[i]) {
            case 'W': count_W++; break;
            case 'Y': count_Y++; break;
            case 'C': count_C++; break;
        }
    }

    // Formula for extinction coefficient at 280nm
    double extinction_coeff = (double)count_W * EPSILON_TRP_PROT +
                              (double)count_Y * EPSILON_TYR_PROT +
                              (double)count_C * EPSILON_CYS_PROT_REDUCED;

    return extinction_coeff;
}

// Calculate GRAVY index for Protein
double calculate_protein_gravy(const char *seq) {
     size_t len = strlen(seq);
    if (len == 0) return 0.0; // Or return error/NaN? 0.0 seems reasonable for empty.

    double total_hydro = 0.0;
    const char *aa_order = "ARNDCEQGHILKMFPSTWYV"; // Order matches HYDROPHOBICITY_SCALE

    for (size_t i = 0; i < len; i++) {
        char aa = seq[i];
        for (int j = 0; j < 20; j++) {
            if (aa == aa_order[j]) {
                total_hydro += HYDROPHOBICITY_SCALE[j];
                break; // Found the amino acid
            }
        }
    }

    // GRAVY is the average hydrophobicity
    return total_hydro / len;
}


// --- Main Program ---

int main(int argc, char *argv[]) {
    // Check command line arguments
    if (argc < 4) {
        printf("Usage: %s <molecule_type> <calculation> <sequence> [options]\n", argv[0]);
        printf("\nSupported Molecule Types: DNA, RNA, Protein (case-insensitive)\n");
        printf("\nSupported Calculations (case-insensitive):\n");
        printf("  For DNA/RNA: length, gc_content, mw, extinction_coefficient\n");
        printf("  For DNA only: tm, revcomp, transcribe\n");
        printf("  For RNA only: translate (to protein)\n");
        printf("  For Protein only: length, mw, composition, pi, extinction_coefficient, gravy\n");
        printf("\nOptions:\n");
        printf("  -conc <molar_concentration> : Specify DNA concentration for Tm calculation (e.g., -conc 50e-9 for 50 nM)\n");
        return 1; // Indicate error
    }

    char *molecule_type_str = argv[1];
    char *calculation_str = argv[2];
    char *sequence = argv[3];

    // Clean the sequence (uppercase, remove whitespace)
    clean_sequence(sequence);

    // Determine molecule type
    int molecule_type = 0; // 0: None, 1: DNA, 2: RNA, 3: Protein
    if (strcasecmp(molecule_type_str, "DNA") == 0) {
        molecule_type = 1;
    } else if (strcasecmp(molecule_type_str, "RNA") == 0) {
        molecule_type = 2;
    } else if (strcasecmp(molecule_type_str, "Protein") == 0) {
        molecule_type = 3;
    } else {
        printf("Error: Invalid molecule type '%s'. Supported types: DNA, RNA, Protein.\n", molecule_type_str);
        return 1;
    }

    // Validate sequence based on type
    int is_valid = 0;
    if (molecule_type == 1 && is_valid_dna(sequence)) is_valid = 1;
    else if (molecule_type == 2 && is_valid_rna(sequence)) is_valid = 1;
    else if (molecule_type == 3 && is_valid_protein(sequence)) is_valid = 1;

    if (!is_valid) {
        printf("Error: Invalid characters in sequence for molecule type '%s'.\n", molecule_type_str);
        // Optionally print expected characters
        if (molecule_type == 1) printf("Expected characters: A, T, C, G\n");
        else if (molecule_type == 2) printf("Expected characters: A, U, C, G\n");
        else if (molecule_type == 3) printf("Expected characters: %s\n", "ARNDCEQGHILKMFPSTWYV");
        return 1;
    }

    // Process calculation command
    char calculation_command[50]; // Buffer for calculation command
    strncpy(calculation_command, calculation_str, sizeof(calculation_command) - 1);
    calculation_command[sizeof(calculation_command) - 1] = '\0'; // Ensure null termination
    // Convert calculation command to lowercase for case-insensitive matching
    for(int i = 0; calculation_command[i]; i++){
      calculation_command[i] = tolower((unsigned char)calculation_command[i]);
    }


    // --- Perform Calculation ---
    int calculation_found = 0;

    if (strcmp(calculation_command, "length") == 0) {
        calculation_found = 1;
        size_t len = calculate_length(sequence);
        printf("Length of sequence: %zu %s\n", len, (molecule_type == 3 ? "amino acids" : "bases"));

    } else if (strcmp(calculation_command, "gc_content") == 0) {
        if (molecule_type == 1 || molecule_type == 2) {
            calculation_found = 1;
            double gc = calculate_gc_content(sequence);
            printf("GC Content: %.2f%%\n", gc);
        } else {
            printf("Error: 'gc_content' calculation is only supported for DNA and RNA.\n");
            return 1;
        }

    } else if (strcmp(calculation_command, "mw") == 0) {
         calculation_found = 1;
         double mw;
         if (molecule_type == 1) {
             mw = calculate_nucleic_mw(sequence, "DNA");
             printf("Molecular Weight (DNA single strand): %.2f Da\n", mw);
         } else if (molecule_type == 2) {
             mw = calculate_nucleic_mw(sequence, "RNA");
             printf("Molecular Weight (RNA single strand): %.2f Da\n", mw);
         } else if (molecule_type == 3) {
             mw = calculate_protein_mw(sequence);
              printf("Molecular Weight (Protein): %.2f Da\n", mw);
         }

    } else if (strcmp(calculation_command, "tm") == 0) {
        if (molecule_type == 1) { // Only for DNA
            calculation_found = 1;
            double concentration = DEFAULT_DNA_CONC; // Default concentration

            // Check for optional concentration argument
            for (int i = 4; i < argc; i += 2) {
                if (strcmp(argv[i], "-conc") == 0 && i + 1 < argc) {
                    char *endptr;
                    concentration = strtod(argv[i+1], &endptr);
                    if (*endptr != '\0' || concentration <= 0) {
                        printf("Error: Invalid or non-positive concentration value provided for -conc.\n");
                        return 1;
                    }
                    break; // Assume only one -conc option
                }
            }

            double tm = calculate_dna_tm(sequence, concentration);
            if (tm > -273.15) { // Tm in Celsius should be above absolute zero
                 printf("Melting Temperature (Tm) [NN Method, %.2e M DNA]: %.2f °C\n", concentration, tm);
            } else {
                // Error message already printed by calculate_dna_tm if tm was -1.0
                 if (tm == -1.0) return 1;
                 printf("Warning: Calculated Tm is extremely low (possibly below absolute zero), check parameters/sequence.\n");
                 printf("Melting Temperature (Tm) [NN Method, %.2e M DNA]: %.2f °C\n", concentration, tm);
            }

        } else {
            printf("Error: 'tm' calculation is only supported for DNA.\n");
            return 1;
        }

    } else if (strcmp(calculation_command, "revcomp") == 0) {
        if (molecule_type == 1) { // Only for DNA
            calculation_found = 1;
            char *revcomp_seq = calculate_dna_revcomp(sequence);
            if (revcomp_seq) {
                printf("Reverse Complement:\n%s\n", revcomp_seq);
                free(revcomp_seq); // Free allocated memory
            } else {
                 return 1; // Error message printed by function
            }
        } else {
            printf("Error: 'revcomp' calculation is only supported for DNA.\n");
            return 1;
        }

    } else if (strcmp(calculation_command, "transcribe") == 0) {
         if (molecule_type == 1) { // Only for DNA
             calculation_found = 1;
             char *rna_seq = transcribe_dna_to_rna(sequence);
             if (rna_seq) {
                 printf("Transcribed RNA Sequence:\n%s\n", rna_seq);
                 free(rna_seq); // Free allocated memory
             } else {
                  return 1; // Error message printed by function
             }
         } else {
             printf("Error: 'transcribe' calculation is only supported for DNA.\n");
             return 1;
         }

    } else if (strcmp(calculation_command, "translate") == 0) {
         if (molecule_type == 2) { // Only for RNA
             calculation_found = 1;
             char *protein_seq = translate_rna_to_protein(sequence);
             if (protein_seq) {
                 printf("Translated Protein Sequence:\n%s\n", protein_seq);
                 free(protein_seq); // Free allocated memory
             } else {
                 // Error message printed by function
                 return 1;
             }
         } else {
             printf("Error: 'translate' calculation is only supported for RNA.\n");
             return 1;
         }

    } else if (strcmp(calculation_command, "composition") == 0) {
         if (molecule_type == 3) { // Only for Protein
             calculation_found = 1;
             int composition_counts[20]; // 20 standard amino acids
             calculate_aa_composition(sequence, composition_counts);

             printf("Amino Acid Composition:\n");
             const char *aa_order = "ARNDCEQGHILKMFPSTWYV";
             for (int i = 0; i < 20; i++) {
                 printf("  %c: %d\n", aa_order[i], composition_counts[i]);
             }
         } else {
             printf("Error: 'composition' calculation is only supported for Protein.\n");
             return 1;
         }

    } else if (strcmp(calculation_command, "pi") == 0) {
        if (molecule_type == 3) { // Only for Protein
            calculation_found = 1;
            double pi = calculate_protein_pi(sequence);
             if (pi >= MIN_PH) { // calculate_protein_pi returns -1.0 on error
                printf("Isoelectric Point (pI): %.2f\n", pi);
            } else {
                 printf("Error: Could not calculate pI for the given sequence.\n");
                 return 1;
            }
        } else {
            printf("Error: 'pi' calculation is only supported for Protein.\n");
            return 1;
        }

    } else if (strcmp(calculation_command, "extinction_coefficient") == 0) {
         calculation_found = 1;
         double extinction_coeff;
         if (molecule_type == 1) {
             extinction_coeff = calculate_nucleic_extinction(sequence, "DNA");
             printf("Molar Extinction Coefficient (DNA, 260nm, single strand): %.2f M^-1 cm^-1\n", extinction_coeff);
         } else if (molecule_type == 2) {
             extinction_coeff = calculate_nucleic_extinction(sequence, "RNA");
             printf("Molar Extinction Coefficient (RNA, 260nm, single strand): %.2f M^-1 cm^-1\n", extinction_coeff);
         } else if (molecule_type == 3) {
             extinction_coeff = calculate_protein_extinction(sequence);
             printf("Molar Extinction Coefficient (Protein, 280nm): %.2f M^-1 cm^-1\n", extinction_coeff);
         }

    } else if (strcmp(calculation_command, "gravy") == 0) {
         if (molecule_type == 3) { // Only for Protein
             calculation_found = 1;
             double gravy = calculate_protein_gravy(sequence);
             printf("GRAVY (Grand Average of Hydropathy) Index: %.2f\n", gravy);
         } else {
             printf("Error: 'gravy' calculation is only supported for Protein.\n");
             return 1;
         }
    }


    // If calculation command was not recognized
    if (!calculation_found) {
        printf("Error: Invalid calculation command '%s' for molecule type '%s'.\n", calculation_str, molecule_type_str);
         printf("\nSupported Calculations (case-insensitive):\n");
        printf("  For DNA/RNA: length, gc_content, mw, extinction_coefficient\n");
        printf("  For DNA only: tm, revcomp, transcribe\n");
        printf("  For RNA only: translate (to protein)\n");
        printf("  For Protein only: length, mw, composition, pi, extinction_coefficient, gravy\n");
        return 1;
    }

    return 0; // Indicate success
}