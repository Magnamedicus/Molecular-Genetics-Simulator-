# Molecular Genetics Simulator

## Overview

The **Molecular Genetics Simulator** is a Python-based tool designed to simulate essential molecular biology processes, such as DNA replication, RNA transcription, and protein translation.
 Using a doubly linked list structure to represent nucleotides, the simulator mimics biochemical molecules and their interactions, providing a computational framework for molecular genetics studies.

This project can be used for educational purposes, algorithm development, or even to explore computational representations of genetic processes.

---

## Features

1. **DNA Representation**:
   - DNA strands are modeled as doubly linked lists, where each node represents a nucleotide (`A`, `T`, `C`, `G`).
   - Nucleotides have directional indicators (`5'` and `3'`) to reflect real DNA strand polarity.
   - Complementary strands are generated with antiparallel orientation.

2. **DNA Replication**:
   - Implements the process of complementary strand synthesis via a `DNA_Polymerase` class.

3. **Gene Assignment**:
   - Allows specific nucleotide sequences within a DNA strand to be marked as "genes."
   - Gene start and end positions are tagged, and their index ranges are stored for reference.

4. **RNA Transcription**:
   - Simulates transcription using the `RNA_Polymerase` class, producing mRNA transcripts from DNA templates.
   - Supports gene-specific transcription based on gene markers.

5. **Protein Translation**:
   - Uses the `Ribosome` class to translate mRNA transcripts into polypeptides (primary amino acid sequences).
   - Supports translation using a standard genetic code.

6. **Mutation Handling**:
   - Includes functionality to mutate specific nucleotides within a DNA strand, allowing for simulation of point mutations.

7. **Custom Outputs**:
   - Provides detailed visual representations of DNA, RNA, and polypeptide sequences.

---

## Classes and Responsibilities

### 1. `Nucleotide`
Represents a nucleotide in a DNA or RNA strand.  
Properties:
- `value`: The nitrogenous base (`A`, `T`, `C`, `G` for DNA; `A`, `U`, `C`, `G` for RNA).
- `next`, `previous`: Links to adjacent nucleotides in the strand.
- `complement`: Points to the complementary nucleotide.
- `end_direction`: Indicates the nucleotide's position (`5'` or `3'` end).
- `gene_marker`: Marks the start or end of a gene sequence.
- `promoter_marker`: Optional feature for marking promoters (future enhancement).

---

### 2. `DNAStrand`
Represents a DNA strand as a doubly linked list.  
Key Methods:
- `append(value)`: Adds a nucleotide to the strand.
- `print_sequence()`: Displays the sequence with directional indicators.
- `generate_complement_sequence()`: Produces a list of complementary nucleotides.
- `point_mutation(index, value)`: Introduces a mutation at the specified index.
- `find_sequence(sequence, allow_overlapping=False)`: Locates specific subsequences in the strand.
- `assign_gene(sequence, gene_name)`: Marks a sequence as a gene and stores its index range.

---

### 3. `DNA_Polymerase`
Simulates DNA replication.  
Key Method:
- `DNA_Replication(strand)`: Generates a complementary DNA strand from a given template.

---

### 4. `RNA_Polymerase`
Simulates RNA transcription.  
Key Method:
- `RNA_Transcription(strand, gene)`: Transcribes a specified gene into an mRNA strand.

---

### 5. `RNAStrand`
Represents an mRNA strand as a doubly linked list.  
Key Methods:
- `append(value)`: Adds a nucleotide to the strand.
- `print_sequence()`: Displays the sequence with directional indicators.

---

### 6. `Ribosome`
Simulates protein translation using the genetic code.  
Key Method:
- `translation(mRNA_transcript)`: Translates an mRNA transcript into a polypeptide.

---

### 7. `Polypeptide`
Represents a protein's primary structure as a doubly linked list of amino acids.  
Key Methods:
- `form_peptide_bond(amino_acid)`: Adds an amino acid to the polypeptide chain.
- `print_polypeptide()`: Displays the primary structure of the polypeptide.

---

## How It Works

### Workflow

1. **DNA Strand Initialization**:
   - Create a DNA strand and append nucleotides to build a sequence.

2. **Gene Assignment**:
   - Identify subsequences in the DNA strand and assign gene markers.

3. **DNA Replication**:
   - Use the `DNA_Polymerase` class to generate a complementary strand.

4. **RNA Transcription**:
   - Transcribe a specific gene from the DNA strand into an mRNA transcript using `RNA_Polymerase`.

5. **Protein Translation**:
   - Translate the mRNA transcript into a polypeptide using the `Ribosome` class.

6. **Visual Output**:
   - Print DNA, RNA, and protein sequences to observe the simulated molecular processes.

---

## Example Execution

```python
def main():
    # Instantiate molecular objects
    dna_strand = DNAStrand()
    dna_polymerase = DNA_Polymerase()
    rna_polymerase = RNA_Polymerase()
    ribosome = Ribosome()

    # Define DNA sequence and append nucleotides
    sequence = "AGGCCGGTTAATGCGTATGCGTATGACTA"
    for nucleotide in sequence:
        dna_strand.append(nucleotide)

    # Generate complement strand
    complement_strand = dna_polymerase.DNA_Replication(dna_strand)

    # Print DNA sequences
    dna_strand.print_sequence()
    complement_strand.print_sequence()

    # Assign genes
    dna_strand.assign_gene("CCGGTTAATG", "bronson-149")
    dna_strand.assign_gene("TAATGCGTATGCGTATGAC", "alf-123")

    # Transcribe mRNA
    mRNA_bronson = rna_polymerase.RNA_Transcription(dna_strand, "bronson-149")
    mRNA_alf = rna_polymerase.RNA_Transcription(dna_strand, "alf-123")

    # Print mRNA sequences
    mRNA_bronson.print_sequence()
    mRNA_alf.print_sequence()

    # Translate proteins
    bronson_peptide = ribosome.translation(mRNA_bronson)
    alf_peptide = ribosome.translation(mRNA_alf)

    # Print polypeptides
    bronson_peptide.print_polypeptide()
    alf_peptide.print_polypeptide()

if __name__ == "__main__":
    main()
```

---

## Requirements

- Python 3.7 or higher.
