
#nucleotide class definition
#this is a doubly linked list, each 'nucleotide' has a 'next', 'previous' and 'complement' property
#all initially set to None except for the nitrogenous base value passed at instantiation

class Nucleotide:
    def __init__(self,value, end_direction = None):
        self.value = value
        self.next = None
        self.previous = None
        self.complement = None
        self.end_direction = end_direction
        self.gene_marker = None
        self.promoter_marker = None

class AminoAcid:
    def __init__(self, value):
        self.value = value
        self.next = None
        self.previous = None




class Ribosome:
    def __init__(self):
        self.genetic_code = {
    'UUU': 'Phe', 'UUC': 'Phe',  # Phenylalanine
    'UUA': 'Leu', 'UUG': 'Leu',  # Leucine
    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',  # Leucine
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile',  # Isoleucine
    'AUG': 'Met',  # Methionine (Start Codon)
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',  # Valine
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',  # Serine
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',  # Proline
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',  # Threonine
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',  # Alanine
    'UAU': 'Tyr', 'UAC': 'Tyr',  # Tyrosine
    'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop',  # Stop Codons
    'CAU': 'His', 'CAC': 'His',  # Histidine
    'CAA': 'Gln', 'CAG': 'Gln',  # Glutamine
    'AAU': 'Asn', 'AAC': 'Asn',  # Asparagine
    'AAA': 'Lys', 'AAG': 'Lys',  # Lysine
    'GAU': 'Asp', 'GAC': 'Asp',  # Aspartic acid
    'GAA': 'Glu', 'GAG': 'Glu',  # Glutamic acid
    'UGU': 'Cys', 'UGC': 'Cys',  # Cysteine
    'UGG': 'Trp',  # Tryptophan
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',  # Arginine
    'AGU': 'Ser', 'AGC': 'Ser',  # Serine
    'AGA': 'Arg', 'AGG': 'Arg',  # Arginine
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'  # Glycine
}

    def translation(self, mRNA_transcript):
        base_count = 0
        codon = ""

        polypeptide = Polypeptide()
        polypeptide.name = mRNA_transcript.name
        temp = mRNA_transcript.head
        while temp:
            codon += temp.value
            if len(codon) == 3:
                if codon in self.genetic_code:
                    polypeptide.form_peptide_bond(self.genetic_code[codon])

                codon = ""
            temp = temp.next

        return polypeptide
















class DNA_Polymerase:
    def __init__(self):
        pass

    def DNA_Replication(self, strand):

        complement_sequence = strand.generate_complement_sequence()

        new_strand = DNAStrand()


        for i in range(len(complement_sequence)):
            new_nucleotide = complement_sequence[i]
            if new_strand.head == None:
                new_strand.head = new_nucleotide
                new_strand.tail = new_nucleotide



            else:
                new_strand.tail.next = new_nucleotide
                new_nucleotide.previous = new_strand.tail
                new_strand.tail = new_nucleotide
            new_strand.length +=1


        return new_strand


class RNA_Polymerase:
    def __init__(self):
        pass

    def RNA_Transcription(self, strand, gene):
        RNA_dict = {'A': 'U', 'T': 'A', 'G': 'C', 'C': 'G'}

        # Find the start of the gene
        temp = strand.head
        while temp:
            if temp.gene_marker == f"{gene} start":
                break
            temp = temp.next

        # Check if the gene start was found
        if not temp or temp.gene_marker != f"{gene} start":
            print(f"Gene {gene} start could not be found in this strand")
            return None


        gene_transcriber = temp
        mRNA_transcript = RNAStrand()
        mRNA_transcript.name = gene


        while gene_transcriber:
            if gene_transcriber.gene_marker == f"{gene} end":
                break

            if gene_transcriber.value in RNA_dict:
                new_nucleotide = Nucleotide(RNA_dict[gene_transcriber.value])


                if mRNA_transcript.head is None:
                    mRNA_transcript.head = new_nucleotide
                    mRNA_transcript.tail = new_nucleotide
                else:
                    mRNA_transcript.tail.next = new_nucleotide
                    new_nucleotide.previous = mRNA_transcript.tail
                    mRNA_transcript.tail = new_nucleotide


                mRNA_transcript.length += 1


            gene_transcriber = gene_transcriber.next

        if gene_transcriber is None or gene_transcriber.gene_marker != f"{gene} end":
            print(f"Gene {gene} end could not be found. Transcription incomplete.")
            return None

        return mRNA_transcript

class Polypeptide:
    def __init__(self, value = None):
        self.name = None

        if value:
            new_amino_acid = AminoAcid(value)
            self.head = new_amino_acid
            self.tail = new_amino_acid
            self.length = 1
        else:
            self.head = None
            self.tail = None
            self.length = 0

    def print_polypeptide(self):
        if self.head is None:
            print("This polypeptide contains no residues")
            return False
        else:
            temp = self.head
            print(f"\n\nPrimary Structure for {self.name} Polypeptide: ")
            while temp is not None:
                if temp.next:
                    print(temp.value, end = "--")
                else:
                    print(temp.value, end = "")
                temp = temp.next
            return True

    def form_peptide_bond(self, amino_acid):
        new_amino_acid = AminoAcid(amino_acid)
        if self.head is None:
            self.head = new_amino_acid
            self.tail = new_amino_acid
        else:
            self.tail.next = new_amino_acid
            new_amino_acid.previous = self.tail
            self.tail = new_amino_acid
        self.length +=1
        return True




class RNAStrand:
    def __init__(self, value=None):
        self.name = None

        if value:
            new_nucleotide = Nucleotide(value)
            self.head = new_nucleotide
            self.tail = new_nucleotide
            self.length = 1
        else:
            self.head = None
            self.tail = None
            self.length = 0

    def print_sequence(self):
        if self.head == None:
            print("This strand contains no nucleotides")
            return False
        else:
            if self.head.end_direction == "5'":
                print("5'", end = ' ')
            if self.head.end_direction == "3'":
                print("3'", end=' ')
            temp = self.head
            print("\n")
            print(f"RNA Transcript for {self.name} gene: ")
            while temp is not None:

                if temp.next is not None:
                    print(temp.value, end = "--")
                else:
                    print(temp.value, end = ' ')

                temp = temp.next
            if self.tail.end_direction == "3'":
                print("3'")
            if self.tail.end_direction == "5'":
                print("5'")



            return True

    def append(self,value):
        new_nucleotide = Nucleotide(value)
        if self.head == None:
            self.head,self.tail = new_nucleotide, new_nucleotide
            new_nucleotide.end_direction = "5'"
        else:
            self.tail.next = new_nucleotide
            new_nucleotide.previous = self.tail
            self.tail = new_nucleotide
            if new_nucleotide.next is None:
                new_nucleotide.end_direction = "3'"

        self.length +=1

    def get_nucleotide_byIndex(self,index):
        if self.head is None:
            print("This strand contains no nucleotides")
            return None

        elif index < 0 or index >= self.length:
                print("index out of range")
                return None
        else:
            if index < self.length/2:
                temp = self.head
                for _ in range(index):
                    temp = temp.next
                #print(f"The nucleotide at index {index} is {temp.value}.")
                return temp
            else:
                temp = self.tail
                for _ in range(self.length -1, index, -1):
                    temp = temp.previous
                #print(f"The nucleotide at index {index} is {temp.value}.")
                return temp









class DNAStrand:
    def __init__(self, value = None):



        if value:
            new_nucleotide = Nucleotide(value)
            self.head = new_nucleotide
            self.tail = new_nucleotide
            self.length = 1
        else:
            self.head = None
            self.tail = None
            self.length = 0
        #self.gene_promoter_indices = {}
        self.wildtype_gene_indices = {}


    def generate_complement_sequence(self):
        DNA_dict = {'A':'T',
                    'T':'A',
                    'G':'C',
                    'C':'G'}

        compliment_sequence = []


        if self.head == None:
            print("this strand contains no nucleotides")
            return False

        else:
            temp = self.head
            while temp is not None:
                if temp.value in DNA_dict:
                    new_nucleotide = Nucleotide(DNA_dict[temp.value])
                    temp.complement = new_nucleotide
                    compliment_sequence.append(new_nucleotide)

                temp = temp.next
            compliment_sequence[0].end_direction = "3'"
            compliment_sequence[-1].end_direction = "5'"

        return compliment_sequence





    def print_sequence(self):
        if self.head == None:
            print("This strand contains no nucleotides")
            return False
        else:

            if self.head.end_direction == "5'":
                print("5'", end = ' ')
            if self.head.end_direction == "3'":
                print("3'", end=' ')
            temp = self.head

            while temp is not None:

                if temp.next is not None:
                    print(temp.value, end = "--")
                else:
                    print(temp.value, end = ' ')

                temp = temp.next
            if self.tail.end_direction == "3'":
                print("3'")
            if self.tail.end_direction == "5'":
                print("5'")



            return True

    def append(self,value):
        new_nucleotide = Nucleotide(value)
        if self.head == None:
            self.head,self.tail = new_nucleotide, new_nucleotide
            new_nucleotide.end_direction = "5'"
        else:
            self.tail.next = new_nucleotide
            new_nucleotide.previous = self.tail
            self.tail = new_nucleotide
            if new_nucleotide.next is None:
                new_nucleotide.end_direction = "3'"

        self.length +=1

    def get_nucleotide_byIndex(self,index):
        if self.head is None:
            print("This strand contains no nucleotides")
            return None

        elif index < 0 or index >= self.length:
                print("index out of range")
                return None
        else:
            if index < self.length/2:
                temp = self.head
                for _ in range(index):
                    temp = temp.next
                #print(f"The nucleotide at index {index} is {temp.value}.")
                return temp
            else:
                temp = self.tail
                for _ in range(self.length -1, index, -1):
                    temp = temp.previous
                #print(f"The nucleotide at index {index} is {temp.value}.")
                return temp

    def point_mutation(self,index,value):
        temp = self.get_nucleotide_byIndex(index)
        if temp:
            temp.value = value
            return True
        return False

    def find_sequence(self, sequence, allow_overlapping=False):

        matches = []
        sequence_len = len(sequence)

        if sequence_len == 0 or sequence_len > self.length:
            return []

        pointer_1 = self.head
        current_index = 0
        while pointer_1:

            temp_pointer = pointer_1
            match_index = 0

            while temp_pointer and match_index < sequence_len and temp_pointer.value == sequence[match_index]:
                temp_pointer = temp_pointer.next
                match_index += 1


            if match_index == sequence_len:
                start_index = current_index
                end_index = current_index + sequence_len - 1
                matches.append((start_index, end_index))

                if not allow_overlapping:
                    
                    pointer_1 = self.get_nucleotide_byIndex(end_index + 1)
                    current_index = end_index + 1
                    continue


            pointer_1 = pointer_1.next
            current_index += 1

        return matches if matches else []

    def assign_gene(self, sequence, gene_name):

        matches = self.find_sequence(sequence)

        if not matches:
            print(f"Gene {gene_name} could not be assigned. Sequence not found.")
            return False

        for match in matches:
            gene_start_index, gene_end_index = match

            start_node = self.get_nucleotide_byIndex(gene_start_index)
            end_node = self.get_nucleotide_byIndex(gene_end_index)

            # Assign gene markers
            start_node.gene_marker = f"{gene_name} start"
            end_node.gene_marker = f"{gene_name} end"
            print(f"Gene {gene_name} assigned from index {gene_start_index} to {gene_end_index}")

            # Store gene positions
            if gene_name not in self.wildtype_gene_indices:
                self.wildtype_gene_indices[gene_name] = []
            self.wildtype_gene_indices[gene_name].append((gene_start_index, gene_end_index))

        return True




def main():
    #instantiating molecular objects
    dna_strand = DNAStrand()
    dna_polymerase = DNA_Polymerase()
    rna_polymerase = RNA_Polymerase()
    ribosome = Ribosome()



    #appending 'sequence' to the initial DNA strand object
    sequence = "AGGCCGGTTAATGCGTATGCGTATGACTA"
    for i in sequence:
        dna_strand.append(i)

    #generating a complement strand by passing the DNA strand object to dna_polymerase's DNA_Replication() method
    complement_strand = dna_polymerase.DNA_Replication(dna_strand)


    #printing template and complement strand sequences
    print("\nDemonstration Chromosome:")
    dna_strand.print_sequence()
    complement_strand.print_sequence()

    #Sub-sequences within the DNA template strand are marked as belonging to a specific gene
    #The first and last nucleotide of the gene sequence is tagged with the gene name
    #The index range where the gene is located is then stored in the strand's 'wild-type dictionary'.
    #This is to provide a point of comparison if subjecting the gene to frameshift mutations

    print("\nDemonstration gene assignments: ")
    dna_strand.assign_gene("CCGGTTAATG", "bronson-149")
    dna_strand.assign_gene("TAATGCGTATGCGTATGAC", "alf-123")

    print("\ninitial index ranges of the assigned gene sequences: ")
    print(dna_strand.wildtype_gene_indices)

    #mRNA transcripts generated by passing the template dna strand and the name of an assigned gene
    #to the RNA_Polymerase object's RNA_Transcription method
    mRNA_bronson = rna_polymerase.RNA_Transcription(dna_strand,"bronson-149")
    mRNA_alf = rna_polymerase.RNA_Transcription(dna_strand, "alf-123")

    #Printing mRNA transcripts
    mRNA_bronson.print_sequence()
    mRNA_alf.print_sequence()

    #polypeptide primary structures generated by passing mRNA transcipts to the ribosome object's translation() method
    alf_123_peptide = ribosome.translation(mRNA_alf)
    bronson_149_peptide = ribosome.translation(mRNA_bronson)

    #printing polypeptide primary structure
    alf_123_peptide.print_polypeptide()
    bronson_149_peptide.print_polypeptide()



main()









