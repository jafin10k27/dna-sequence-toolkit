"""
Command-Line Interface (CLI) Module
====================================
Provides an intuitive command-line entry point for the DNA Sequence Analysis Toolkit.

Exit Codes:
    0: Success
    1: Data validation or FASTA parsing error
    2: Command-line syntax or argument error
    3: I/O or filesystem error
"""

import argparse
import sys
from typing import List, Optional

from dna_toolkit import __version__
from dna_toolkit.analysis import analyze_composition
from dna_toolkit.operations import complement, reverse_complement, translate_dna
from dna_toolkit.parser import (
    EmptyFastaError,
    FastaError,
    FastaFileNotFoundError,
    MalformedFastaError,
    parse_fasta,
)
from dna_toolkit.reporting import (
    RecordReport,
    export_csv,
    export_json,
    format_terminal_report,
)
from dna_toolkit.validation import (
    DuplicateRecordError,
    EmptySequenceError,
    InvalidBaseError,
    ValidationError,
    ValidationMode,
    validate_records,
)

EXIT_SUCCESS = 0
EXIT_DATA_ERROR = 1
EXIT_ARG_ERROR = 2
EXIT_IO_ERROR = 3


def build_parser() -> argparse.ArgumentParser:
    """Build and configure the argparse ArgumentParser."""
    parser = argparse.ArgumentParser(
        prog="dna-toolkit",
        description="DNA Sequence Analysis Toolkit: parse, validate, analyze, and translate DNA FASTA files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run complete analysis and display report in terminal
  dna-toolkit data/sample.fasta

  # Run only composition and GC analysis
  dna-toolkit data/sample.fasta --composition

  # Export full analysis to JSON
  dna-toolkit data/sample.fasta --format json -o output/results.json

  # Export to CSV with reading frame 2 and vertebrate mitochondrial code (table 2)
  dna-toolkit data/sample.fasta --translate --frame 2 --table 2 --format csv -o output/results.csv

  # Allow IUPAC ambiguous nucleotides (e.g. N, R, Y)
  dna-toolkit data/ambiguous.fasta --mode ambiguous
""",
    )

    parser.add_argument(
        "fasta_file",
        nargs="?",
        help="Path to the input FASTA sequence file.",
    )
    parser.add_argument(
        "-i", "--input",
        dest="input_file",
        help="Alternative flag to specify input FASTA file path.",
    )
    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    # Operations
    ops_group = parser.add_argument_group("Analysis Operations (defaults to all if none selected)")
    ops_group.add_argument(
        "--all",
        action="store_true",
        help="Run full analysis (composition, GC-content, reverse complement, and translation).",
    )
    ops_group.add_argument(
        "--composition",
        action="store_true",
        help="Calculate sequence length, nucleotide frequencies, and GC percentage.",
    )
    ops_group.add_argument(
        "--reverse-complement",
        action="store_true",
        help="Generate Watson-Crick reverse complement (5' to 3').",
    )
    ops_group.add_argument(
        "--complement",
        action="store_true",
        help="Generate direct strand complement (3' to 5').",
    )
    ops_group.add_argument(
        "--translate",
        action="store_true",
        help="Translate DNA sequence to protein.",
    )

    # Validation options
    val_group = parser.add_argument_group("Validation Options")
    val_group.add_argument(
        "-m", "--mode",
        choices=["strict", "ambiguous"],
        default="strict",
        help="Validation policy: 'strict' allows only A, T, G, C; 'ambiguous' allows IUPAC codes (default: strict).",
    )
    val_group.add_argument(
        "--allow-duplicates",
        action="store_true",
        help="Allow duplicate record identifiers in FASTA file without raising an error.",
    )

    # Translation options
    trans_group = parser.add_argument_group("Translation Options")
    trans_group.add_argument(
        "-f", "--frame",
        type=int,
        choices=[1, 2, 3, -1, -2, -3],
        default=1,
        help="Reading frame for translation: +1, +2, +3 for forward; -1, -2, -3 for reverse strand (default: 1).",
    )
    trans_group.add_argument(
        "-t", "--table",
        default="1",
        help="NCBI genetic code table ID or name (default: '1' for Standard Code).",
    )
    trans_group.add_argument(
        "--to-stop",
        action="store_true",
        help="Halt translation at the first stop codon encountered.",
    )

    # Output options
    out_group = parser.add_argument_group("Output and Export Options")
    out_group.add_argument(
        "--format",
        choices=["text", "json", "csv"],
        default="text",
        help="Output format: 'text' (human-readable), 'json', or 'csv' (default: text).",
    )
    out_group.add_argument(
        "-o", "--output",
        help="Optional path to output file. If omitted, results are printed to stdout.",
    )

    return parser


def run_pipeline(args: argparse.Namespace) -> int:
    """Execute the core analysis workflow based on parsed arguments."""
    input_path = args.fasta_file or args.input_file
    if not input_path:
        sys.stderr.write("Error: Missing input FASTA file. Provide path as argument or with -i/--input.\n")
        sys.stderr.write("Run 'dna-toolkit --help' for usage instructions.\n")
        return EXIT_ARG_ERROR

    # Determine genetic code table (cast to int if digits)
    table = int(args.table) if args.table.isdigit() else args.table

    # Determine active operations
    has_specific_ops = any([
        args.composition,
        args.reverse_complement,
        args.complement,
        args.translate,
    ])
    do_all = args.all or not has_specific_ops

    do_comp = do_all or args.composition
    do_rev_comp = do_all or args.reverse_complement
    do_complement = do_all or args.complement
    do_translate = do_all or args.translate

    # 1. Parse FASTA
    try:
        raw_records = parse_fasta(input_path)
    except FastaFileNotFoundError as err:
        sys.stderr.write(f"File Error: {err}\n")
        return EXIT_IO_ERROR
    except (EmptyFastaError, MalformedFastaError) as err:
        sys.stderr.write(f"FASTA Format Error: {err}\n")
        return EXIT_DATA_ERROR
    except Exception as err:
        sys.stderr.write(f"Unexpected I/O Error reading FASTA: {err}\n")
        return EXIT_IO_ERROR

    # 2. Validate records
    try:
        validated_records = validate_records(
            raw_records,
            mode=args.mode,
            allow_duplicates=args.allow_duplicates,
        )
    except (EmptySequenceError, InvalidBaseError, DuplicateRecordError) as err:
        sys.stderr.write(f"Validation Error: {err}\n")
        return EXIT_DATA_ERROR
    except ValidationError as err:
        sys.stderr.write(f"Validation Error: {err}\n")
        return EXIT_DATA_ERROR

    # 3. Perform analysis and operations
    reports: List[RecordReport] = []
    for rec in validated_records:
        report = RecordReport(
            record_id=rec.id,
            description=rec.description,
            sequence=rec.sequence,
        )

        if do_comp:
            report.composition = analyze_composition(rec, mode=args.mode)

        if do_complement:
            report.complement = complement(rec.sequence, mode=args.mode)

        if do_rev_comp:
            report.reverse_complement = reverse_complement(rec.sequence, mode=args.mode)

        if do_translate:
            try:
                report.translation = translate_dna(
                    rec.sequence,
                    frame=args.frame,
                    table=table,
                    to_stop=args.to_stop,
                    mode=args.mode,
                    record_id=rec.id,
                )
            except Exception as err:
                sys.stderr.write(f"Translation Error for record '{rec.id}': {err}\n")
                return EXIT_DATA_ERROR

        reports.append(report)

    # 4. Generate Output and Export
    metadata = {
        "input_file": input_path,
        "mode": args.mode,
        "frame": args.frame,
        "table": str(table),
    }

    try:
        if args.format == "json":
            result_str = export_json(reports, output_path=args.output, metadata=metadata)
            if not args.output:
                sys.stdout.write(result_str + "\n")
            else:
                sys.stdout.write(f"Successfully exported JSON report to: {args.output}\n")

        elif args.format == "csv":
            result_str = export_csv(reports, output_path=args.output)
            if not args.output:
                sys.stdout.write(result_str)
            else:
                sys.stdout.write(f"Successfully exported CSV report to: {args.output}\n")

        else:  # text
            report_text = format_terminal_report(reports)
            if args.output:
                from dna_toolkit.reporting import _safe_write_file
                _safe_write_file(args.output, report_text)
                sys.stdout.write(f"Successfully wrote report to: {args.output}\n")
            else:
                sys.stdout.write(report_text + "\n")

    except (PermissionError, OSError) as err:
        sys.stderr.write(f"Output File Error: {err}\n")
        return EXIT_IO_ERROR

    return EXIT_SUCCESS


def main(argv: Optional[List[str]] = None) -> int:
    """CLI entry point function."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return run_pipeline(args)


if __name__ == "__main__":
    sys.exit(main())
