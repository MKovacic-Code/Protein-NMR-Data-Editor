# Residue Renumberer & Chemical Shift Referencing Correction Tool

A lightweight Python tool featuring a simple, native system GUI and a headless CLI for batch renumbering residue indices, translating nomenclatures, and applying referencing corrections (offsets) to protein structures and NMR chemical shift files.

The tool preserves exact file alignment, trailing comments, and column spacing across all operations.

---

## Features

- **Native GUI Theme**: Basic, lightweight, and clean system-native theme (no heavy styling or external CSS dependencies) with side-by-side original/modified preview screens.
- **CLI Mode (Headless)**: Fully featured command-line interface supporting headless automation and scripting pipelines.
- **Batch Processing**: Select and process multiple files simultaneously, with a dedicated file-selector dropdown in the GUI to preview any file in the queue.
- **Nomenclature Translation**: Automated mapping between IUPAC and CYANA/DYANA atom naming standards (supporting alpha, beta, and branch-chain protons across standard amino acids).
- **Referencing Corrections**: Apply individual chemical shift offsets in ppm for Proton ($^1$H), Carbon ($^{13}$C), and Nitrogen ($^{15}$N).
- **Out-of-Bounds Validation**: Automatic biological feasibility validation against standard BMRB ranges (warns in GUI logs and CLI if shift values exceed bounds: $^1$H [0-12 ppm], $^{13}$C [10-220 ppm], $^{15}$N [90-140 ppm]).
- **Alignment-Preserving Spacer**: Dynamically adjusts whitespace padding downstream when token lengths change (e.g. `9` -> `109`) to keep table columns perfectly aligned.
- **Zero External Dependencies**: Built entirely using Python's standard library (`tkinter`, `argparse`, and `re`).

---

## File Format Support

1. **NEF (`.nef`)**:
   - Updates `sequence_code` across sequence, chemical shift, and distance restraint lists.
   - Applies chemical shift offsets using `_nef_chemical_shift.element` (with fallback to `_nef_chemical_shift.atom_name` prefixes if elements are unspecified).
2. **NMR-STAR (`.txt`, `.str`, `.star`)**:
   - Updates `_Atom_chem_shift.Comp_index_ID`.
   - Modifies chemical shift values (`_Atom_chem_shift.Val`) and identifies target nuclei from `_Atom_chem_shift.Atom_ID`.
3. **PDB (`.pdb`)**:
   - Renumbers residue sequence numbers in `ATOM`, `HETATM`, `ANISOU`, and `TER` records.
   - Translates atom names according to the selected nomenclature standard.
4. **TALOS (`.tab`, `.talos`)**:
   - Parses VARS columns to renumber `RESID` indices and correct chemical shifts for specified nuclei.
5. **Sparky (`.list`, `.sparky`, `.peaks`)**:
   - Parses multi-dimensional assignment labels (e.g. `G759H-G759N` -> `G859H-G859N`), renumbers them, maps atom nomenclatures, and offsets corresponding shift columns (`w1`, `w2`, `w3`).
6. **Yasara (`.tbl`)**:
   - Fallback regex mode to renumber `resid <num>` and `<num> <amino_acid>` definitions while preserving file spacing.

---

## Getting Started

### Prerequisites
- Python 3.6 or later.

### Running the GUI App
On Windows, you can double-click the `run.bat` script or use the `Residue Renumberer.lnk` shortcut.
Alternatively, launch from your terminal:
```bash
python renumber_residues.py
```

### Running in Headless CLI Mode
Provide arguments directly to the script to run it without launching the GUI:
```bash
python renumber_residues.py -i input1.nef input2.pdb -o ./output_dir -r 100 --proton 0.010 -n IUPAC_TO_CYANA
```

**CLI Options:**
- `-i`, `--input`: One or more input file paths (separated by spaces).
- `-o`, `--output`: One or more output file paths or a single output directory.
- `-r`, `--residue`: Residue index offset integer (e.g., `100` or `-50`).
- `--proton`, `--carbon`, `--nitrogen`: Chemical shift offset values in ppm.
- `--no-detect`: Disable automatic format detection and rely strictly on file extensions.
- `-n`, `--nomenclature`: Nomenclature mapping standard (`NONE`, `IUPAC_TO_CYANA`, or `CYANA_TO_IUPAC`).

---

## Project Structure

- `renumber_residues.py`: Main GUI and CLI script.
- `run.bat`: Batch script shortcut to launch the app on Windows.
- `README.md`: Project documentation.
