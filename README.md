# BioSequence Analyzer: Technical Documentation

## 1. Overview

**BioSequence Analyzer** is a command-line utility written in C designed to perform fundamental bioinformatics calculations on biological sequences. It supports three primary molecule types: **DNA**, **RNA**, and **Proteins**. The tool provides functionalities ranging from basic sequence statistics (length, GC content) to complex biophysical properties such as Melting Temperature ($T_m$ using Nearest-Neighbor methods), Isoelectric Point (pI), and Molar Extinction Coefficients.

This document outlines the installation, usage, supported calculations, and underlying algorithms implemented in `main.c`.

---

## 2. Compilation and Installation

The program requires a standard C compiler (e.g., `gcc`) and links against the math library (`libm`).

### Prerequisites
*   GCC or Clang compiler
*   Standard C Library (`stdio.h`, `stdlib.h`, `string.h`, `ctype.h`, `math.h`)

### Build Instructions
To compile the source code into an executable named `bioanalyzer`:

```bash
gcc -o bioanalyzer main.c -lm
```

*Note: The `-lm` flag is required to link the math library for functions like `log()`, `pow()`, and `fabs()`.*

---

## 3. Usage

### Syntax
```bash
./bioanalyzer <molecule_type> <calculation> <sequence> [options]
```

### Arguments
| Argument | Description | Valid Values |
| :--- | :--- | :--- |
| `molecule_type` | The type of biological sequence provided. | `DNA`, `RNA`, `Protein` (Case-insensitive) |
| `calculation` | The specific analysis to perform. | See Section 4 (Case-insensitive) |
| `sequence` | The raw sequence string. | Alphanumeric characters corresponding to the molecule type. Whitespace is automatically stripped. |
| `options` | Optional flags for specific calculations. | `-conc <molar_concentration>` (For DNA $T_m$ only) |

### Example Commands
1.  **Calculate DNA Melting Temperature:**
    ```bash
    ./bioanalyzer DNA tm ATCGATCGATCG -conc 50e-9
    ```
2.  **Translate RNA to Protein:**
    ```bash
    ./bioanalyzer RNA translate AUGGCCUAA
    ```
3.  **Calculate Protein Isoelectric Point:**
    ```bash
    ./bioanalyzer Protein pi MKWVTFISLLFLFSSAYSR
    ```

---

## 4. Supported Calculations

The available calculations depend on the selected `molecule_type`.

### 4.1. Common Calculations (DNA, RNA, Protein)
*   **`length`**: Returns the number of bases (nucleotides) or amino acids.
*   **`mw`**: Calculates the Molecular Weight in Daltons (Da).
    *   *Nucleic Acids:* Sum of residue weights + mass of one water molecule ($H_2O$) for terminal groups.
    *   *Proteins:* Sum of amino acid residue weights (minus $H_2O$ per peptide bond) + mass of one water molecule for terminals.
*   **`extinction_coefficient`**: Calculates the molar extinction coefficient ($\epsilon$).
    *   *DNA/RNA:* At 260 nm (single-stranded). Based on sum of individual nucleotide contributions.
    *   *Protein:* At 280 nm. Based on Trp, Tyr, and Cys (reduced) content (Pace et al., 1995).

### 4.2. DNA-Specific Calculations
*   **`gc_content`**: Percentage of Guanine and Cytosine bases in the sequence.
*   **`tm`**: Melting Temperature ($^\circ C$) calculated using the **Nearest-Neighbor (NN) Thermodynamic Model** (SantaLucia, 1998).
    *   *Default Concentration:* 50 nM ($50 \times 10^{-9}$ M).
    *   *Custom Concentration:* Use `-conc <value>` (in Molar).
    *   *Formula:* $T_m = \frac{\Delta H \cdot 1000}{\Delta S + R \cdot \ln(C_t / 4)} - 273.15$
*   **`revcomp`**: Generates the Reverse Complement sequence.
*   **`transcribe`**: Converts DNA sequence to RNA (replaces Thymine 'T' with Uracil 'U').

### 4.3. RNA-Specific Calculations
*   **`translate`**: Translates RNA sequence into a Protein sequence using the standard genetic code.
    *   Translation starts at the first codon.
    *   Stops at the first Stop Codon (`UAA`, `UAG`, `UGA`).
    *   If the sequence length is not a multiple of 3, trailing bases are ignored.

### 4.4. Protein-Specific Calculations
*   **`composition`**: Displays the count of each of the 20 standard amino acids.
*   **`pi`**: Isoelectric Point. Calculated via iterative pH scanning (step 0.01) to find the pH where the net charge is closest to zero.
    *   Considers ionizable side chains (D, E, C, Y, H, K, R) and N/C termini.
*   **`gravy`**: Grand Average of Hydropathy Index. Calculated using the Kyte & Doolittle scale. Positive values indicate hydrophobicity; negative values indicate hydrophilicity.

---

## 5. Algorithmic Details & Constants

### 5.1. Nearest-Neighbor Parameters (DNA $T_m$)
The tool uses thermodynamic parameters from *SantaLucia, J. (1998)*.
*   $\Delta H$ (Enthalpy): kcal/mol
*   $\Delta S$ (Entropy): cal/mol/K
*   Gas Constant ($R$): 1.987 cal/mol/K

### 5.2. Genetic Code
Translation uses a standard lookup table mapping RNA codons (UUU, UUC, etc.) to single-letter amino acid codes. Start codon (AUG) translates to Methionine (M).

### 5.3. Protein pI Calculation
*   **Method:** Iterative search from pH 0.0 to 14.0.
*   **Step Size:** 0.01 pH units.
*   **pKa Values Used:**
    *   N-Term: 8.0, C-Term: 3.1
    *   Asp: 4.1, Glu: 4.5, Cys: 6.8, Tyr: 10.5
    *   His: 6.0, Lys: 10.5, Arg: 12.5

### 5.4. Extinction Coefficients
*   **Nucleotides (260 nm):** A=15400, T=8700, C=7300, G=11500, U=9900 $M^{-1}cm^{-1}$.
*   **Proteins (280 nm):** Trp=5690, Tyr=1280, Cys(reduced)=120 $M^{-1}cm^{-1}$.

---

## 6. Error Handling and Validation

*   **Input Cleaning:** The `clean_sequence` function automatically converts input to uppercase and removes whitespace/newlines.
*   **Character Validation:**
    *   **DNA:** Accepts only A, T, C, G.
    *   **RNA:** Accepts only A, U, C, G.
    *   **Protein:** Accepts only the 20 standard amino acid letters (ARNDCEQGHILKMFPSTWYV).
*   **Invalid Inputs:** The program returns exit code `1` and prints an error message to `stdout` or `stderr` if:
    *   Invalid molecule type is specified.
    *   Sequence contains invalid characters for the specified type.
    *   Incompatible calculation is requested (e.g., `tm` for Protein).
    *   Memory allocation fails.

---

## 7. Limitations

1.  **Salt Correction:** The $T_m$ calculation uses the basic Nearest-Neighbor model without explicit salt concentration corrections (beyond the standard formula assumptions).
2.  **Secondary Structure:** $T_m$ and Extinction Coefficient calculations assume linear, single-stranded molecules without secondary structure interactions (hairpins, dimers).
3.  **Protein pI Precision:** The pI is approximated to two decimal places based on a fixed step size iteration. It does not account for local environmental effects on pKa values within a folded protein.
4.  **Memory Management:** Sequences for `revcomp`, `transcribe`, and `translate` are dynamically allocated. The caller (main function) is responsible for freeing this memory.

---

## 8. References

1.  SantaLucia, J. (1998). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. *PNAS*, 95(4), 1460-1465.
2.  Pace, C. N., et al. (1995). How to measure and predict the molar absorption coefficient of a protein. *Protein Science*, 4, 2411-2423.
3.  Kyte, J., & Doolittle, R. F. (1982). A simple method for displaying the hydropathic character of a protein. *J. Mol. Biol.*, 157(1), 105-132.
4.  Bjellqvist, B., et al. (2004). Isoelectric point prediction of proteins. *Bioinformatics*.
