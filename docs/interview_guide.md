# Bioinformatics & Python Technical Interview Guide

This guide prepares you to explain and defend every design decision, algorithmic choice, and biological concept implemented in the DNA Sequence Analysis Toolkit during technical interviews.

---

## 1. Python Architecture & Language Questions

### Q1: Why did you use `@dataclass(frozen=True)` instead of plain dictionaries or standard classes?
**Answer**:
- **Immutability and Thread Safety**: In genomic pipelines, sequence objects pass through multiple analysis and reporting stages. `frozen=True` prevents accidental in-place mutation (e.g., modifying `record.sequence` during one operation before passing it to another).
- **Memory & Performance**: Dataclasses with defined attributes are lightweight and provide auto-generated `__init__`, `__repr__`, and `__eq__` without boilerplate.
- **Type Safety**: Paired with type annotations, dataclasses provide explicit contracts that static analyzers (like `mypy`) and IDEs can verify, unlike dictionary keys which can easily suffer from typos (e.g. `'length'` vs `'len'`).

### Q2: How does Biopython's `SeqIO.parse` manage memory, and why does that matter for genome files?
**Answer**:
- `SeqIO.parse()` returns a Python **generator (iterator)** rather than loading all records into an in-memory list at once.
- Whole-genome sequence files (like the human reference genome or multi-gigabyte FASTQ/FASTA files) can easily exceed available RAM. By yielding records on demand, memory consumption remains $\mathcal{O}(1)$ relative to the total file size during streaming.

### Q3: Why is the package structured using the `src` layout instead of placing the package directly in the root directory?
**Answer**:
- The **`src` layout** isolates application code from the current working directory.
- In a flat root layout (`./dna_toolkit`), running `python` or `pytest` from the project root will import the local directory via `sys.path[0]` without exercising the packaging metadata or installation configuration defined in `pyproject.toml`.
- The `src` layout forces tests to import the package as an installed package (e.g. via `pip install -e .`), guaranteeing that entry points and package discovery work as intended.

### Q4: Why create custom exception classes (`FastaFileNotFoundError`, `EmptySequenceError`, `InvalidBaseError`)?
**Answer**:
- Custom exception hierarchies allow callers and CLI coordinators to distinguish between **I/O errors** (file not found), **formatting errors** (not a FASTA file), and **biological data validity errors** (invalid bases).
- It also allows exceptions to carry structured context: for example, `InvalidBaseError` stores `invalid_chars` and `char_positions`, allowing the user-facing error message to point to the exact base and index where corruption occurred.

---

## 2. Molecular Biology & Bioinformatics Questions

### Q5: What is the biological significance of GC-content, and how is it used in genomics?
**Answer**:
- **Hydrogen Bonding**: Guanine (G) and Cytosine (C) form **three hydrogen bonds**, whereas Adenine (A) and Thymine (T) form **two hydrogen bonds**.
- **Thermostability & PCR**: DNA strands with higher GC-content require more thermal energy to denature (separate). In Polymerase Chain Reaction (PCR) and primer design, GC-content determines the annealing and melting temperatures ($T_m$).
- **Genomic Architecture**: Across genomes, GC-content varies systematically:
  - Higher GC-content is frequently correlated with gene-dense regions and CpG islands in vertebrate promoters.
  - Extremophiles (such as thermophilic bacteria living in hydrothermal vents) frequently exhibit elevated GC-content to stabilize their DNA at high temperatures.

### Q6: What is the difference between direct complement and reverse complement?
**Answer**:
- DNA double-stranded helices are **antiparallel**: one strand runs $5' \to 3'$, and its complementary strand runs $3' \to 5'$.
- For sequence $5'\text{-AGTC-}3'$:
  - **Direct Complement**: Forms Watson-Crick pairings ($A \leftrightarrow T$, $G \leftrightarrow C$) yielding $3'\text{-TCAG-}5'$.
  - **Reverse Complement**: Because molecular sequences are universally read and reported $5' \to 3'$, we reverse the direct complement: $5'\text{-GACT-}3'$.
- To study the transcribed sequence or coding sequence of a gene residing on the opposite strand, bioinformaticians must always take the **reverse complement**, not the direct complement.

### Q7: Why are there six reading frames, and why can't we assume a sequence starts at Frame 1?
**Answer**:
- A codon consists of 3 consecutive nucleotides.
- On the forward strand, translation can begin at index 0 (Frame +1), index 1 (Frame +2), or index 2 (Frame +3). An offset of 3 brings you back to Frame +1 shifted by one codon.
- Because DNA is double-stranded and genes can be encoded on either strand, the reverse complement strand also contains three possible reading frames (Frame -1, Frame -2, Frame -3).
- FASTA files frequently store arbitrary sequencing fragments, PCR amplicons, or genomic contigs rather than pre-trimmed mature mRNAs starting with an ATG start codon. Translating without specifying or exploring the reading frame can produce non-functional nonsense peptides full of artificial stop codons.

### Q8: What are IUPAC degenerate nucleotide codes, and why must a validation module support both strict and ambiguous modes?
**Answer**:
- Standard canonical DNA contains only **A, T, G, C**.
- High-throughput sequencing generates ambiguous base calls when signal quality is low (e.g. `N` for any nucleotide).
- Heterozygous diploid organisms contain single nucleotide polymorphisms (SNPs): for example, an individual with an A on one chromosome and a G on the homologous chromosome is designated by IUPAC code `R` (puRine).
- In **strict mode**, an analysis pipeline enforces high quality by discarding or halting on any non-canonical bases.
- In **ambiguous mode**, pipelines permit IUPAC degenerate characters, but the analysis module must track degenerate counts and compute unambiguous GC-content so that `N` or other codes do not artificially distort biological calculations.

### Q9: What is an alternative genetic code, and when does Biopython's translation table parameter matter?
**Answer**:
- While the **Standard Genetic Code (NCBI Table 1)** is conserved across most nuclear eukaryotic and prokaryotic genes, it is not universal.
- Organelles with independent evolutionary histories, such as mitochondria, use alternative code tables:
  - In **Vertebrate Mitochondrial Code (Table 2)**, `TGA` codes for Tryptophan (`W`) instead of acting as a Stop codon (`*`), and `AGA`/`AGG` code for Stop instead of Arginine.
  - In **Bacterial, Archaeal, and Plant Plastid Code (Table 11)**, alternative start codons (like `GTG` and `TTG`) are frequently used.
- A robust bioinformatics toolkit must allow specifying the genetic code table rather than hardcoding the standard nuclear code.

---

## 3. Testing, Verification & Code Quality

### Q10: Why did you write tests using independent biological benchmarks rather than testing against the code's output?
**Answer**:
- Testing against the code's own output creates a **circular test**: if the implementation contains a biological bug (e.g. miscalculating GC percentage or confusing complement with reverse complement), saving the flawed output into the test will cause the buggy code to "pass."
- Instead, tests were written using independently verified, hand-calculated biological ground truth:
  - `ATGC`: length 4, GC% = 50.0%.
  - `GGCC`: length 4, GC% = 100.0%.
  - `ATAT`: length 4, GC% = 0.0%.
  - `GATTACA`: counts A=3, T=2, G=1, C=1, GC% = 28.57%.
  - `AGTC`: reverse complement `GACT`.
  - `ATGGAA`: translated to `ME` (Methionine, Glutamic acid).

### Q11: How did you handle incomplete codons during translation?
**Answer**:
- If a sequence length in a given reading frame is not a multiple of 3 (e.g. 7 nucleotides = 2 codons + 1 leftover base), the leftover bases cannot form a complete amino acid.
- Rather than crashing or silently swallowing warnings, our toolkit:
  1. Computes the remainder: `remainder = len(subseq) % 3`.
  2. Trims the trailing incomplete bases to pass only full triplets to the translation engine.
  3. Records `partial_bases_trimmed` in the `TranslationResult` metadata so researchers know exactly how many nucleotides were omitted.

### Q12: How are newlines and special characters handled in the CSV export?
**Answer**:
- In FASTA files, description lines can contain commas or line breaks.
- In accordance with **RFC 4180**, the toolkit uses Python's standard `csv.DictWriter`, which properly quotes fields containing delimiters or special characters.
- In addition, the export engine normalizes whitespace in sequence strings and descriptions (`" ".join(desc.split())`), ensuring that every FASTA record corresponds to exactly one row in the output spreadsheet.
