import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import re
import os
import sys
import argparse

# Amino acids and nucleic acids pattern for Yasara .tbl fallback parser
AMINO_ACIDS = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", 
               "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "CA",
               "DA", "DT", "DG", "DC", "A", "U", "G", "C", "T", "RA", "RU", "RG", "RC",
               "ADE", "THY", "GUA", "CYT", "URA"]
AA_PATTERN = "|".join(AMINO_ACIDS)
FALLBACK_PATTERN1 = re.compile(rf"([ \t]*)\b(\d+) ({AA_PATTERN})\b")
FALLBACK_PATTERN2 = re.compile(r"\bresid([ \t]+)(\d+)([ \t]*)")

# Standard nomenclature mapping of common differences (IUPAC <-> CYANA/DYANA)
IUPAC_TO_CYANA = {
    # Alpha protons (Glycine HA2/HA3 -> CYANA HA1/HA2)
    "HA2": "HA1", "HA3": "HA2",
    # Methylene protons (IUPAC HB2/HB3 -> CYANA HB1/HB2)
    "HB2": "HB1", "HB3": "HB2",
    "HG2": "HG1", "HG3": "HG2",
    "HD2": "HD1", "HD3": "HD2",
    "HE2": "HE1", "HE3": "HE2",
    # Isoleucine CG1 methylene (IUPAC HG12/HG13 -> CYANA HG11/HG12)
    "HG12": "HG11", "HG13": "HG12",
    # Valine CG1/CG2 methyls (IUPAC HG11/12/13 -> CYANA QG1)
    "HG11": "QG1", "HG12": "QG1", "HG13": "QG1",
    "HG21": "QG2", "HG22": "QG2", "HG23": "QG2",
    # Leucine CD1/CD2 methyls
    "HD11": "QD1", "HD12": "QD1", "HD13": "QD1",
    "HD21": "QD2", "HD22": "QD2", "HD23": "QD2",
    # Isoleucine CD1/CG2 methyls
    "HD11": "QD1", "HD12": "QD1", "HD13": "QD1",
    "HG21": "QG2", "HG22": "QG2", "HG23": "QG2",
    # Amide protons of Arg/Lys etc.
    "HH11": "QH1", "HH12": "QH1", "HH21": "QH2", "HH22": "QH2",
}

# Reverse mapping: CYANA/DYANA -> IUPAC
CYANA_TO_IUPAC = {
    # Alpha protons (CYANA HA1/HA2 -> IUPAC HA2/HA3)
    "HA1": "HA2", "HA2": "HA3",
    # Methylene protons (CYANA HB1/HB2 -> IUPAC HB2/HB3)
    "HB1": "HB2", "HB2": "HB3",
    "HG1": "HG2", "HG2": "HG3",
    "HD1": "HD2", "HD2": "HD3",
    "HE1": "HE2", "HE2": "HE3",
    # Isoleucine CG1 methylene (CYANA HG11/HG12 -> IUPAC HG12/HG13)
    "HG11": "HG12", "HG12": "HG13",
    # Methyl groups
    "QB": "HB1", "HB*": "HB1",
    "QG1": "HG11", "HG1*": "HG11",
    "QG2": "HG21", "HG2*": "HG21",
    "QD1": "HD11", "HD1*": "HD11",
    "QD2": "HD21", "HD2*": "HD21",
    "QH1": "HH11", "QH2": "HH21",
}

ONE_TO_THREE = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
}


def is_nucleic_acid(residue_name, atom_name):
    """
    Checks if a residue-atom pair is in a nucleic acid context.
    """
    res_upper = residue_name.strip().upper() if residue_name else ""
    atom_upper = atom_name.strip().upper() if atom_name else ""
    
    # 1. Uniquely nucleic residue names
    uniquely_nucleic_res = {
        "DA", "DT", "DG", "DC", "RA", "RU", "RG", "RC",
        "ADE", "THY", "GUA", "CYT", "URA", "U"
    }
    if res_upper in uniquely_nucleic_res:
        return True
        
    # 2. Check for prime in atom name (sugar atoms)
    if "'" in atom_upper:
        return True
        
    # 3. Handle ambiguous one-letter codes (A, G, C, T) by checking atom names
    if res_upper in {"A", "G", "C", "T"}:
        # Unique sugar and base atoms for nucleic acids
        nucleic_atoms = {
            "H8", "H6", "H5", "H2", "H1", "H3", "H21", "H22", "H41", "H42", "H61", "H62",
            "N9", "N7", "N3", "N1", "O4", "O6", "O2", "P", "OP1", "OP2", "O5'", "O3'", "O4'",
            "C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''"
        }
        if atom_upper in nucleic_atoms:
            return True
        # Regex check for primed atoms or CYANA-style prime notations (like H2'1, H2'2, H5'1, H5'2)
        if re.search(r"^[C|H|O|N][1-5]'", atom_upper) or re.search(r"^H[25]'[12]$", atom_upper):
            return True
            
    return False


def validate_shift(val, element, residue_name, atom_name, line_num, warnings):
    """
    Validates chemical shifts against biological bounds, adjusting for protein vs nucleic acid context.
    """
    is_nucleic = is_nucleic_acid(residue_name or "", atom_name or "")
    
    if element == 'H':
        max_h = 15.0 if is_nucleic else 12.0
        if val < 0.0 or val > max_h:
            warnings.append(f"Line {line_num}: Proton shift {val:.3f} ppm is out of biological bounds (0.0 - {max_h} ppm)")
    elif element == 'C':
        if val < 10.0 or val > 220.0:
            warnings.append(f"Line {line_num}: Carbon shift {val:.3f} ppm is out of biological bounds (10.0 - 220.0 ppm)")
    elif element == 'N':
        max_n = 260.0 if is_nucleic else 140.0
        if val < 90.0 or val > max_n:
            warnings.append(f"Line {line_num}: Nitrogen shift {val:.3f} ppm is out of biological bounds (90.0 - {max_n} ppm)")


def map_atom_name(atom_name, amino_acid, mode):
    """
    Maps atom name between IUPAC and CYANA/DYANA standards.
    """
    if not mode or mode == "NONE":
        return atom_name

    atom_upper = atom_name.strip().upper()
    aa_upper = amino_acid.strip().upper() if amino_acid else ""
    
    is_nucleic = is_nucleic_acid(aa_upper, atom_upper)
    
    if not is_nucleic and len(aa_upper) == 1:
        aa_upper = ONE_TO_THREE.get(aa_upper, aa_upper)

    if mode == "IUPAC_TO_CYANA":
        if is_nucleic:
            # Sugar proton mappings
            if atom_upper == "H2''": return "H2'2"
            elif atom_upper == "H5'": return "H5'1"
            elif atom_upper == "H5''": return "H5'2"
            elif atom_upper == "H2'":
                # Only map H2' to H2'1 if it is DNA context.
                # RNA has single H2' (so we don't map it to H2'1 if residue starts with 'R' or is 'U')
                if not (aa_upper.startswith('R') or aa_upper == 'U'):
                    return "H2'1"
            return atom_name

        # Check amino-acid specific methyls
        if aa_upper == 'ALA' and atom_upper in ['HB1', 'HB2', 'HB3']:
            return 'QB'
        elif aa_upper == 'VAL':
            if atom_upper in ['HG11', 'HG12', 'HG13']: return 'QG1'
            elif atom_upper in ['HG21', 'HG22', 'HG23']: return 'QG2'
        elif aa_upper == 'LEU':
            if atom_upper in ['HD11', 'HD12', 'HD13']: return 'QD1'
            elif atom_upper in ['HD21', 'HD22', 'HD23']: return 'QD2'
        elif aa_upper == 'ILE':
            if atom_upper in ['HD11', 'HD12', 'HD13']: return 'QD1'
            elif atom_upper in ['HG21', 'HG22', 'HG23']: return 'QG2'
            elif atom_upper == 'HG12': return 'HG11'
            elif atom_upper == 'HG13': return 'HG12'
        elif aa_upper == 'GLY':
            if atom_upper == 'HA2': return 'HA1'
            elif atom_upper == 'HA3': return 'HA2'

        # General/fallback matching
        return IUPAC_TO_CYANA.get(atom_upper, atom_name)

    elif mode == "CYANA_TO_IUPAC":
        if is_nucleic:
            # Sugar proton mappings
            if atom_upper == "H2'2": return "H2''"
            elif atom_upper == "H2'1": return "H2'"
            elif atom_upper == "H5'1": return "H5'"
            elif atom_upper == "H5'2": return "H5''"
            return atom_name

        # Check amino-acid specific methyls
        if aa_upper == 'ALA':
            if atom_upper in ['QB', 'HB*']:
                return 'HB1'
            elif atom_upper in ['HB1', 'HB2', 'HB3']:
                return atom_upper
        elif aa_upper == 'VAL':
            if atom_upper in ['QG1', 'HG1*']: return 'HG11'
            elif atom_upper in ['QG2', 'HG2*']: return 'HG21'
        elif aa_upper == 'LEU':
            if atom_upper in ['QD1', 'HD1*']: return 'HD11'
            elif atom_upper in ['QD2', 'HD2*']: return 'HD21'
        elif aa_upper == 'ILE':
            if atom_upper in ['QD1', 'HD1*']: return 'HD11'
            elif atom_upper in ['QG2', 'HG2*']: return 'HG21'
            elif atom_upper == 'HG11': return 'HG12'
            elif atom_upper == 'HG12': return 'HG13'
        elif aa_upper == 'GLY':
            if atom_upper == 'HA1': return 'HA2'
            elif atom_upper == 'HA2': return 'HA3'

        return CYANA_TO_IUPAC.get(atom_upper, atom_name)

    return atom_name


def detect_file_format(filepath):
    """
    Scans the header and first few lines of a file to auto-detect its NMR/biological format.
    """
    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            head = [f.readline() for _ in range(50)]
        content = "".join(head)
        
        # PDB detection
        if any(line.startswith(('ATOM  ', 'HETATM', 'CRYST1', 'HEADER', 'SEQRES', 'ANISOU')) for line in head):
            return 'pdb'
            
        # TALOS detection
        if 'DATA SEQUENCE' in content or ('VARS' in content and 'FORMAT' in content):
            return 'talos'
            
        # NEF detection
        if 'data_nef_' in content or '_nef_' in content:
            return 'nef'
            
        # NMR-STAR detection
        if 'data_' in content and ('_Atom_chem_shift' in content or 'loop_' in content):
            return 'txt'
            
        # Sparky detection
        if any(re.search(r'Assignment\s+w\d+\s+w\d+', line) for line in head) or 'Assignment' in content:
            return 'sparky'
            
        # Yasara TBL detection
        if any('assign' in line and 'resid' in line for line in head):
            return 'fallback'
            
    except Exception:
        pass
    
    # Fallback to extension if content scanning fails
    _, ext = os.path.splitext(filepath.lower())
    if ext == '.nef': return 'nef'
    elif ext in ('.txt', '.str', '.star'): return 'txt'
    elif ext == '.pdb': return 'pdb'
    elif ext in ('.tab', '.talos'): return 'talos'
    elif ext in ('.list', '.sparky', '.peak'): return 'sparky'
    elif ext == '.tbl': return 'fallback'
    else: return 'fallback'


def apply_fallback_renumbering(line, offset):
    """
    Applies regex renumbering logic for Yasara .tbl or other formats.
    """
    def repl_pattern1(match):
        spaces = match.group(1)
        num_str = match.group(2)
        aa = match.group(3)
        num = int(num_str)
        new_num = num + offset
        
        original_len = len(spaces) + len(num_str)
        new_num_str = str(new_num)
        padded_num = new_num_str.rjust(original_len)
        if spaces and not (padded_num.startswith(' ') or padded_num.startswith('\t')):
            padded_num = " " + new_num_str
        return f"{padded_num} {aa}"
        
    def repl_pattern2(match):
        num_str = match.group(2)
        spaces_after = match.group(3)
        num = int(num_str)
        new_num = num + offset
        new_num_str = str(new_num)
        
        len_diff = len(new_num_str) - len(num_str)
        if len_diff > 0:
            remove_count = min(len_diff, len(spaces_after))
            new_spaces_after = spaces_after[remove_count:]
        elif len_diff < 0:
            new_spaces_after = spaces_after + " " * (-len_diff)
        else:
            new_spaces_after = spaces_after
            
        return f"resid {new_num_str}{new_spaces_after}"

    line = FALLBACK_PATTERN1.sub(repl_pattern1, line)
    line = FALLBACK_PATTERN2.sub(repl_pattern2, line)
    return line


def process_pdb_line(line, residue_offset, nomenclature_mode):
    """
    Processes a single line of a PDB file, applying residue offsets and nomenclature mapping.
    """
    if line.startswith(('ATOM  ', 'HETATM', 'ANISOU', 'TER   ')):
        res_num_str = line[22:26]
        res_name = line[17:20].strip()
        atom_name = line[12:16].strip()
        
        # 1. Update Residue Sequence Number
        try:
            res_num = int(res_num_str.strip())
            new_res_num = res_num + residue_offset
            new_res_str = f"{new_res_num:4d}"
            if len(new_res_str) > 4:
                new_res_str = new_res_str[-4:]  # standard PDB limits
        except ValueError:
            new_res_str = res_num_str
            
        # 2. Update Atom Nomenclature if needed
        if nomenclature_mode != 'NONE':
            new_atom_name = map_atom_name(atom_name, res_name, nomenclature_mode)
            if new_atom_name != atom_name:
                if len(new_atom_name) >= 4:
                    atom_str = f"{new_atom_name[:4]}"
                else:
                    atom_str = f" {new_atom_name:<3}"
                line = line[:12] + atom_str + line[16:]
                
        line = line[:22] + new_res_str + line[26:]
    return line


def process_talos(lines, residue_offset, proton_offset, carbon_offset, nitrogen_offset, nomenclature_mode, warnings):
    """
    Processes TALOS .tab files.
    """
    processed_lines = []
    vars_headers = []
    resid_idx = -1
    resname_idx = -1
    shift_cols = {}
    
    for i, line in enumerate(lines):
        stripped = line.strip()
        if not stripped or stripped.startswith('#'):
            processed_lines.append(line)
            continue
            
        parts = stripped.split()
        if parts[0] == 'VARS':
            vars_headers = parts[1:]
            for idx, h in enumerate(vars_headers):
                h_upper = h.upper()
                if h_upper == 'RESID':
                    resid_idx = idx
                elif h_upper == 'RESNAME':
                    resname_idx = idx
                elif h_upper in ['HN', 'HA', 'H', 'HA2', 'HA3', 'HB', 'HD', 'HE']:
                    shift_cols[idx] = ('H', 3)
                elif h_upper in ['CA', 'CB', 'C', 'CO']:
                    shift_cols[idx] = ('C', 3)
                elif h_upper in ['N', 'NH']:
                    shift_cols[idx] = ('N', 3)
            processed_lines.append(line)
            continue
            
        if parts[0] == 'FORMAT' or parts[0].startswith('DATA'):
            processed_lines.append(line)
            continue
            
        if vars_headers and len(parts) >= len(vars_headers):
            lead_ws = ""
            match_lead = re.match(r'^(\s*)', line)
            if match_lead:
                lead_ws = match_lead.group(1)
                
            remain = line[len(lead_ws):]
            match_nl = re.search(r'(\r?\n)$', remain)
            newline = match_nl.group(1) if match_nl else ""
            if newline:
                remain = remain[:-len(newline)]
                
            comment = ""
            comment_idx = remain.find('#')
            if comment_idx != -1:
                comment = remain[comment_idx:]
                remain = remain[:comment_idx]
                
            line_parts = re.split(r'(\s+)', remain)
            tokens_info = []
            for idx, part in enumerate(line_parts):
                if idx % 2 == 0:
                    if part != '':
                        tokens_info.append((part, idx))
                        
            # Apply offsets and nomenclature
            resname = ""
            if resname_idx != -1 and resname_idx < len(tokens_info):
                resname = tokens_info[resname_idx][0]

            for col_idx, (token, part_idx) in enumerate(tokens_info):
                new_token = token
                if col_idx == resid_idx:
                    try:
                        new_token = str(int(token) + residue_offset)
                    except ValueError:
                        pass
                elif col_idx == resname_idx and nomenclature_mode != 'NONE':
                    # Residue name nomenclature map not commonly required in TALOS, but maps if applicable
                    pass
                elif col_idx in shift_cols:
                    nuc, dec = shift_cols[col_idx]
                    atom_name = vars_headers[col_idx] if vars_headers and col_idx < len(vars_headers) else ""
                    offset = 0.0
                    if nuc == 'H': offset = proton_offset
                    elif nuc == 'C': offset = carbon_offset
                    elif nuc == 'N': offset = nitrogen_offset
                    
                    try:
                        val = float(token)
                        new_val = val + offset
                        new_token = f"{new_val:.{dec}f}"
                        
                        # Validation checks
                        validate_shift(new_val, nuc, resname, atom_name, i + 1, warnings)
                    except ValueError:
                        pass
                        
                if new_token != token:
                    len_diff = len(new_token) - len(token)
                    line_parts[part_idx] = new_token
                    if part_idx + 1 < len(line_parts):
                        whitespace = line_parts[part_idx + 1]
                        if len_diff > 0:
                            remove_count = min(len_diff, len(whitespace) - 1)
                            if remove_count > 0:
                                line_parts[part_idx + 1] = whitespace[remove_count:]
                        elif len_diff < 0:
                            pad_char = whitespace[0] if whitespace else ' '
                            line_parts[part_idx + 1] = whitespace + pad_char * (-len_diff)
                            
            processed_lines.append(lead_ws + "".join(line_parts) + comment + newline)
        else:
            processed_lines.append(line)
            
    return processed_lines


def process_sparky_assignment_token(token, residue_offset, nomenclature_mode):
    """
    Parses and renumbers Sparky assignment strings (e.g. G759H-G759N -> G859H-G859N).
    """
    match = re.match(r'^([A-Za-z]*)(\d+)([A-Za-z0-9#%*]*)$', token)
    if match:
        aa = match.group(1)
        res_num_str = match.group(2)
        atom = match.group(3)
        
        try:
            res_num = int(res_num_str)
            new_res_num = res_num + residue_offset
            new_res_str = str(new_res_num)
        except ValueError:
            new_res_str = res_num_str
            
        if nomenclature_mode != 'NONE':
            new_atom = map_atom_name(atom, aa, nomenclature_mode)
        else:
            new_atom = atom
            
        return f"{aa}{new_res_str}{new_atom}"
    return token


def process_sparky(lines, residue_offset, proton_offset, carbon_offset, nitrogen_offset, nomenclature_mode, warnings):
    """
    Processes Sparky peak lists (.list).
    """
    processed_lines = []
    assign_col = -1
    dim_cols = {}
    
    for i, line in enumerate(lines):
        stripped = line.strip()
        if not stripped or stripped.startswith('#'):
            processed_lines.append(line)
            continue
            
        parts = stripped.split()
        if any(w in parts for w in ['Assignment', 'w1', 'w2']):
            for idx, token in enumerate(parts):
                if token.lower() == 'assignment':
                    assign_col = idx
                elif token.lower() == 'w1':
                    dim_cols[0] = idx
                elif token.lower() == 'w2':
                    dim_cols[1] = idx
                elif token.lower() == 'w3':
                    dim_cols[2] = idx
            processed_lines.append(line)
            continue
            
        if assign_col != -1 and len(parts) > assign_col:
            lead_ws = ""
            match_lead = re.match(r'^(\s*)', line)
            if match_lead:
                lead_ws = match_lead.group(1)
                
            remain = line[len(lead_ws):]
            match_nl = re.search(r'(\r?\n)$', remain)
            newline = match_nl.group(1) if match_nl else ""
            if newline:
                remain = remain[:-len(newline)]
                
            line_parts = re.split(r'(\s+)', remain)
            tokens_info = []
            for idx, part in enumerate(line_parts):
                if idx % 2 == 0:
                    if part != '':
                        tokens_info.append((part, idx))
                        
            if assign_col < len(tokens_info):
                orig_assign, assign_part_idx = tokens_info[assign_col]
                
                sep = '-'
                if '/' in orig_assign:
                    sep = '/'
                
                sub_parts = orig_assign.split(sep)
                new_sub_parts = []
                nuclei_types = []
                sparky_atoms_info = []
                
                for part in sub_parts:
                    new_sub_parts.append(process_sparky_assignment_token(part, residue_offset, nomenclature_mode))
                    match = re.match(r'^([A-Za-z]*)(\d+)([A-Za-z0-9#%*]*)$', part)
                    nuc = None
                    aa = ""
                    atom = ""
                    if match:
                        aa = match.group(1)
                        atom = match.group(3)
                        atom_clean = atom.upper()
                        atom_clean = re.sub(r'[%*xy]$', '', atom_clean)
                        if atom_clean.startswith('H') or atom_clean.startswith('Q') or atom_clean.startswith('M'):
                            nuc = 'H'
                        elif atom_clean.startswith('C'):
                            nuc = 'C'
                        elif atom_clean.startswith('N'):
                            nuc = 'N'
                    nuclei_types.append(nuc)
                    sparky_atoms_info.append((aa, atom))
                    
                new_assign = sep.join(new_sub_parts)
                
                # Update assignment token
                len_diff = len(new_assign) - len(orig_assign)
                line_parts[assign_part_idx] = new_assign
                if assign_part_idx + 1 < len(line_parts):
                    whitespace = line_parts[assign_part_idx + 1]
                    if len_diff > 0:
                        remove_count = min(len_diff, len(whitespace) - 1)
                        if remove_count > 0:
                            line_parts[assign_part_idx + 1] = whitespace[remove_count:]
                    elif len_diff < 0:
                        pad_char = whitespace[0] if whitespace else ' '
                        line_parts[assign_part_idx + 1] = whitespace + pad_char * (-len_diff)
                
                # Update chemical shift columns
                for dim_idx, col_idx in dim_cols.items():
                    if col_idx < len(tokens_info) and dim_idx < len(nuclei_types):
                        nuc = nuclei_types[dim_idx]
                        if nuc:
                            offset = 0.0
                            if nuc == 'H': offset = proton_offset
                            elif nuc == 'C': offset = carbon_offset
                            elif nuc == 'N': offset = nitrogen_offset
                            
                            val_token, val_part_idx = tokens_info[col_idx]
                            try:
                                val = float(val_token)
                                new_val = val + offset
                                decimals = 3
                                if '.' in val_token:
                                    decimals = len(val_token.split('.')[1])
                                new_token = f"{new_val:.{decimals}f}"
                                
                                # Validation checks
                                aa, atom = sparky_atoms_info[dim_idx] if dim_idx < len(sparky_atoms_info) else ("", "")
                                validate_shift(new_val, nuc, aa, atom, i + 1, warnings)
                                    
                                len_diff_val = len(new_token) - len(val_token)
                                line_parts[val_part_idx] = new_token
                                if val_part_idx + 1 < len(line_parts):
                                    whitespace = line_parts[val_part_idx + 1]
                                    if len_diff_val > 0:
                                        remove_count = min(len_diff_val, len(whitespace) - 1)
                                        if remove_count > 0:
                                            line_parts[val_part_idx + 1] = whitespace[remove_count:]
                                    elif len_diff_val < 0:
                                        pad_char = whitespace[0] if whitespace else ' '
                                        line_parts[val_part_idx + 1] = whitespace + pad_char * (-len_diff_val)
                            except ValueError:
                                pass
                                
            processed_lines.append(lead_ws + "".join(line_parts) + newline)
        else:
            processed_lines.append(line)
            
    return processed_lines


def process_data_line_tokens(original_line, line_num, mod_cols, residue_offset, proton_offset, carbon_offset, nitrogen_offset, nomenclature_mode, warnings):
    """
    Applies updates to tokens inside NEF and STAR lines while maintaining column spacings.
    """
    match_lead = re.match(r'^(\s*)', original_line)
    lead_ws = match_lead.group(1) if match_lead else ""
    
    remain = original_line[len(lead_ws):]
    match_nl = re.search(r'(\r?\n)$', remain)
    newline = match_nl.group(1) if match_nl else ""
    if newline:
        remain = remain[:-len(newline)]
        
    comment = ""
    comment_idx = remain.find('#')
    if comment_idx != -1:
        comment = remain[comment_idx:]
        remain = remain[:comment_idx]
        
    parts = re.split(r'(\s+)', remain)
    
    tokens_info = []
    for idx, part in enumerate(parts):
        if idx % 2 == 0:
            if part != '':
                tokens_info.append((part, idx))
                
    for col_idx, mod in mod_cols.items():
        if col_idx < len(tokens_info):
            old_token, part_idx = tokens_info[col_idx]
            mod_type = mod.get('type')
            
            if mod_type == 'residue':
                try:
                    val = int(old_token)
                    new_token = str(val + residue_offset)
                except ValueError:
                    new_token = old_token
            elif mod_type == 'atom_name':
                resname_col = mod.get('resname_col', -1)
                resname = ""
                if resname_col != -1 and resname_col < len(tokens_info):
                    resname = tokens_info[resname_col][0]
                new_token = map_atom_name(old_token, resname, nomenclature_mode)
            elif mod_type == 'chem_shift':
                el_col = mod.get('el_col', -1)
                atom_col = mod.get('atom_col', -1)
                
                element = ""
                atom_id = ""
                if el_col != -1 and el_col < len(tokens_info):
                    element = tokens_info[el_col][0]
                if atom_col != -1 and atom_col < len(tokens_info):
                    atom_id = tokens_info[atom_col][0]
                    
                element_clean = element.strip().upper()
                if (not element_clean or element_clean == '.') and atom_id:
                    atom_clean = atom_id.strip().upper()
                    atom_clean = re.sub(r'[%*xy]$', '', atom_clean)
                    if atom_clean.startswith('H') or atom_clean.startswith('Q') or atom_clean.startswith('M'):
                        element_clean = 'H'
                    elif atom_clean.startswith('C'):
                        element_clean = 'C'
                    elif atom_clean.startswith('N'):
                        element_clean = 'N'
                        
                offset = 0.0
                if element_clean == 'H':
                    offset = proton_offset
                elif element_clean == 'C':
                    offset = carbon_offset
                elif element_clean == 'N':
                    offset = nitrogen_offset
                    
                try:
                    val = float(old_token)
                    new_val = val + offset
                    decimals = 3
                    if '.' in old_token:
                        decimals = len(old_token.split('.')[1])
                    new_token = f"{new_val:.{decimals}f}"
                    
                    # Validation Checks
                    resname_col = mod.get('resname_col', -1)
                    resname = ""
                    if resname_col != -1 and resname_col < len(tokens_info):
                        resname = tokens_info[resname_col][0]
                    validate_shift(new_val, element_clean, resname, atom_id, line_num, warnings)
                except ValueError:
                    new_token = old_token
            else:
                new_token = old_token
                
            if new_token != old_token:
                len_diff = len(new_token) - len(old_token)
                parts[part_idx] = new_token
                
                if part_idx + 1 < len(parts):
                    whitespace = parts[part_idx + 1]
                    if len_diff > 0:
                        remove_count = min(len_diff, len(whitespace) - 1)
                        if remove_count > 0:
                            parts[part_idx + 1] = whitespace[remove_count:]
                    elif len_diff < 0:
                        pad_char = whitespace[0] if whitespace else ' '
                        parts[part_idx + 1] = whitespace + pad_char * (-len_diff)
                        
    modified_remain = "".join(parts)
    return lead_ws + modified_remain + comment + newline


def process_lines(lines, residue_offset, proton_offset, carbon_offset, nitrogen_offset, file_type, nomenclature_mode, warnings):
    """
    Main state machine dispatcher for all file types.
    """
    processed_lines = []
    state = 'NORMAL'
    loop_headers = []
    mod_cols = {}
    
    if file_type == 'pdb':
        for i, line in enumerate(lines):
            processed_lines.append(process_pdb_line(line, residue_offset, nomenclature_mode))
        return processed_lines
    elif file_type == 'talos':
        return process_talos(lines, residue_offset, proton_offset, carbon_offset, nitrogen_offset, nomenclature_mode, warnings)
    elif file_type == 'sparky':
        return process_sparky(lines, residue_offset, proton_offset, carbon_offset, nitrogen_offset, nomenclature_mode, warnings)
        
    for i, line in enumerate(lines):
        stripped = line.strip()
        
        if file_type == 'fallback':
            processed_lines.append(apply_fallback_renumbering(line, residue_offset))
            continue
            
        if state == 'NORMAL':
            if stripped == 'loop_':
                state = 'HEADER'
                loop_headers = []
                processed_lines.append(line)
            else:
                processed_lines.append(line)
                
        elif state == 'HEADER':
            if stripped.startswith('_'):
                loop_headers.append(stripped)
                processed_lines.append(line)
            elif stripped == '' or stripped.startswith('#'):
                processed_lines.append(line)
            else:
                state = 'DATA'
                mod_cols = {}
                
                if file_type == 'nef':
                    seq_indices = []
                    val_idx = -1
                    el_idx = -1
                    atom_indices = []
                    resname_indices = []
                    
                    for idx, h in enumerate(loop_headers):
                        if re.search(r'\b_nef_.*sequence_code(?:_\d+)?\b', h):
                            seq_indices.append(idx)
                        elif h.endswith('.value'):
                            val_idx = idx
                        elif h.endswith('.element'):
                            el_idx = idx
                        elif h.endswith('.atom_name') or h.endswith('.atom_name_1') or h.endswith('.atom_name_2'):
                            atom_indices.append(idx)
                        elif h.endswith('.residue_name') or h.endswith('.residue_name_1') or h.endswith('.residue_name_2'):
                            resname_indices.append(idx)
                            
                    for idx in seq_indices:
                        mod_cols[idx] = {'type': 'residue'}
                    if val_idx != -1:
                        loop_atom_idx = -1
                        for idx, h in enumerate(loop_headers):
                            if h.endswith('.atom_name'):
                                loop_atom_idx = idx
                        loop_resname_idx = resname_indices[0] if resname_indices else -1
                        mod_cols[val_idx] = {'type': 'chem_shift', 'el_col': el_idx, 'atom_col': loop_atom_idx, 'resname_col': loop_resname_idx}
                        
                    for atom_col in atom_indices:
                        suffix = ""
                        atom_name_str = loop_headers[atom_col]
                        if '_' in atom_name_str:
                            suffix_parts = atom_name_str.split('_')
                            if suffix_parts[-1].isdigit():
                                suffix = "_" + suffix_parts[-1]
                        
                        matching_resname_col = -1
                        for resname_col in resname_indices:
                            resname_str = loop_headers[resname_col]
                            if resname_str.endswith(suffix):
                                matching_resname_col = resname_col
                                break
                        mod_cols[atom_col] = {'type': 'atom_name', 'resname_col': matching_resname_col}
                        
                elif file_type == 'txt':
                    seq_indices = []
                    val_idx = -1
                    atom_idx = -1
                    comp_idx = -1
                    
                    for idx, h in enumerate(loop_headers):
                        if h.endswith('.Comp_index_ID'):
                            seq_indices.append(idx)
                        elif h.endswith('.Val'):
                            val_idx = idx
                        elif h.endswith('.Atom_ID'):
                            atom_idx = idx
                        elif h.endswith('.Comp_ID'):
                            comp_idx = idx
                            
                    for idx in seq_indices:
                        mod_cols[idx] = {'type': 'residue'}
                    if val_idx != -1:
                        mod_cols[val_idx] = {'type': 'chem_shift', 'el_col': -1, 'atom_col': atom_idx, 'resname_col': comp_idx}
                    if atom_idx != -1:
                        mod_cols[atom_idx] = {'type': 'atom_name', 'resname_col': comp_idx}
                
                if mod_cols:
                    processed_line = process_data_line_tokens(line, i + 1, mod_cols, residue_offset, proton_offset, carbon_offset, nitrogen_offset, nomenclature_mode, warnings)
                    processed_lines.append(processed_line)
                else:
                    processed_lines.append(line)
                    
        elif state == 'DATA':
            if stripped == 'stop_':
                state = 'NORMAL'
                processed_lines.append(line)
            elif stripped.startswith('save_') or stripped.startswith('_'):
                state = 'NORMAL'
                processed_lines.append(line)
            elif stripped == '' or stripped.startswith('#'):
                processed_lines.append(line)
            else:
                if mod_cols:
                    processed_line = process_data_line_tokens(line, i + 1, mod_cols, residue_offset, proton_offset, carbon_offset, nitrogen_offset, nomenclature_mode, warnings)
                    processed_lines.append(processed_line)
                else:
                    processed_lines.append(line)
                    
    return processed_lines


class RenumberApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Residue Renumbering & Referencing Correction Tool")
        self.root.geometry("950x800")
        
        self.input_files = []
        
        self.setup_ui()

    def setup_ui(self):
        # Top panel for configurations
        top_frame = tk.Frame(self.root, padx=10, pady=10)
        top_frame.pack(fill=tk.X)
        
        top_frame.columnconfigure(1, weight=1)
        
        # 1. Inputs Browser
        tk.Label(top_frame, text="Input File(s):").grid(row=0, column=0, sticky="e", pady=5)
        self.input_entry = tk.Entry(top_frame)
        self.input_entry.grid(row=0, column=1, padx=5, pady=5, sticky="we")
        tk.Button(top_frame, text="Browse...", command=self.browse_input).grid(row=0, column=2, pady=5)
        
        # 2. Outputs / Suffix
        tk.Label(top_frame, text="Output Directory/File:").grid(row=1, column=0, sticky="e", pady=5)
        self.output_entry = tk.Entry(top_frame)
        self.output_entry.grid(row=1, column=1, padx=5, pady=5, sticky="we")
        tk.Button(top_frame, text="Browse...", command=self.browse_output).grid(row=1, column=2, pady=5)
        
        # 3. Offsets panel
        offsets_frame = tk.LabelFrame(top_frame, text="Configurations", padx=10, pady=5)
        offsets_frame.grid(row=2, column=0, columnspan=3, sticky="we", pady=10)
        
        tk.Label(offsets_frame, text="Residue Offset:").grid(row=0, column=0, sticky="e", pady=5)
        self.offset_entry = tk.Entry(offsets_frame, width=8)
        self.offset_entry.insert(0, "0")
        self.offset_entry.grid(row=0, column=1, sticky="w", padx=5, pady=5)
        
        tk.Label(offsets_frame, text="Proton (1H) Offset (ppm):").grid(row=0, column=2, sticky="e", pady=5, padx=(15, 0))
        self.offset_h_entry = tk.Entry(offsets_frame, width=8)
        self.offset_h_entry.insert(0, "0.0")
        self.offset_h_entry.grid(row=0, column=3, sticky="w", padx=5, pady=5)
        
        tk.Label(offsets_frame, text="Carbon (13C) Offset (ppm):").grid(row=0, column=4, sticky="e", pady=5, padx=(15, 0))
        self.offset_c_entry = tk.Entry(offsets_frame, width=8)
        self.offset_c_entry.insert(0, "0.0")
        self.offset_c_entry.grid(row=0, column=5, sticky="w", padx=5, pady=5)
        
        tk.Label(offsets_frame, text="Nitrogen (15N) Offset (ppm):").grid(row=0, column=6, sticky="e", pady=5, padx=(15, 0))
        self.offset_n_entry = tk.Entry(offsets_frame, width=8)
        self.offset_n_entry.insert(0, "0.0")
        self.offset_n_entry.grid(row=0, column=7, sticky="w", padx=5, pady=5)
        
        # 4. Advanced Settings Panel
        settings_frame = tk.Frame(top_frame)
        settings_frame.grid(row=3, column=0, columnspan=3, sticky="we", pady=5)
        
        self.autodetect_var = tk.BooleanVar(value=True)
        self.autodetect_cb = tk.Checkbutton(settings_frame, text="Auto-detect format", variable=self.autodetect_var, command=self.toggle_force_format)
        self.autodetect_cb.grid(row=0, column=0, sticky="w", pady=5)
        
        tk.Label(settings_frame, text="Force Format:").grid(row=0, column=1, sticky="e", padx=(15, 0))
        self.force_format_combo = ttk.Combobox(settings_frame, values=["NEF", "NMR-STAR (txt)", "Yasara TBL", "PDB", "TALOS", "Sparky"], width=15, state="disabled")
        self.force_format_combo.current(0)
        self.force_format_combo.grid(row=0, column=2, sticky="w", padx=5)
        
        tk.Label(settings_frame, text="Nomenclature Mapping:").grid(row=0, column=3, sticky="e", padx=(25, 0))
        self.nomenclature_combo = ttk.Combobox(settings_frame, values=["NONE", "IUPAC_TO_CYANA", "CYANA_TO_IUPAC"], width=18, state="readonly")
        self.nomenclature_combo.current(0)
        self.nomenclature_combo.grid(row=0, column=4, sticky="w", padx=5)
        
        # 5. File preview selector (for batch mode)
        self.preview_label = tk.Label(settings_frame, text="Preview File:")
        self.preview_label.grid(row=0, column=5, sticky="e", padx=(25, 0))
        self.preview_combo = ttk.Combobox(settings_frame, width=25, state="readonly")
        self.preview_combo.grid(row=0, column=6, sticky="w", padx=5)
        self.preview_combo.bind("<<ComboboxSelected>>", lambda e: self.load_preview())
        
        tk.Button(top_frame, text="Load Preview", command=self.load_preview).grid(row=3, column=2, sticky="we", pady=5, padx=5)

        # Middle panel for file previews
        middle_frame = tk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        middle_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        # Left Panel (Original)
        left_lf = tk.LabelFrame(middle_frame, text="Original View")
        left_scroll_y = tk.Scrollbar(left_lf)
        left_scroll_y.pack(side=tk.RIGHT, fill=tk.Y)
        left_scroll_x = tk.Scrollbar(left_lf, orient=tk.HORIZONTAL)
        left_scroll_x.pack(side=tk.BOTTOM, fill=tk.X)
        self.text_original = tk.Text(left_lf, wrap=tk.NONE, yscrollcommand=left_scroll_y.set, xscrollcommand=left_scroll_x.set)
        self.text_original.pack(fill=tk.BOTH, expand=True)
        left_scroll_y.config(command=self.text_original.yview)
        left_scroll_x.config(command=self.text_original.xview)
        middle_frame.add(left_lf)
        
        # Right Panel (Modified Preview)
        right_lf = tk.LabelFrame(middle_frame, text="Modified Preview")
        right_scroll_y = tk.Scrollbar(right_lf)
        right_scroll_y.pack(side=tk.RIGHT, fill=tk.Y)
        right_scroll_x = tk.Scrollbar(right_lf, orient=tk.HORIZONTAL)
        right_scroll_x.pack(side=tk.BOTTOM, fill=tk.X)
        self.text_preview = tk.Text(right_lf, wrap=tk.NONE, yscrollcommand=right_scroll_y.set, xscrollcommand=right_scroll_x.set)
        self.text_preview.pack(fill=tk.BOTH, expand=True)
        right_scroll_y.config(command=self.text_preview.yview)
        right_scroll_x.config(command=self.text_preview.xview)
        middle_frame.add(right_lf)
        
        # Bottom Validation Log Pane
        bottom_lf = tk.LabelFrame(self.root, text="Validation Log / Alerts", height=120)
        bottom_lf.pack(fill=tk.X, padx=10, pady=5)
        
        log_scroll_y = tk.Scrollbar(bottom_lf)
        log_scroll_y.pack(side=tk.RIGHT, fill=tk.Y)
        self.text_log = tk.Text(bottom_lf, height=5, yscrollcommand=log_scroll_y.set)
        self.text_log.pack(fill=tk.BOTH, expand=True)
        log_scroll_y.config(command=self.text_log.yview)
        
        # Action Bar at bottom
        action_frame = tk.Frame(self.root, padx=10, pady=10)
        action_frame.pack(fill=tk.X)
        
        tk.Button(action_frame, text="Process and Save File(s)", command=self.process_files, font=("", 10, "bold"), padx=10).pack(side=tk.LEFT)
        self.status_label = tk.Label(action_frame, text="Ready", fg="blue")
        self.status_label.pack(side=tk.LEFT, padx=15)

    def toggle_force_format(self):
        if self.autodetect_var.get():
            self.force_format_combo.config(state="disabled")
        else:
            self.force_format_combo.config(state="readonly")

    def browse_input(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        filenames = filedialog.askopenfilenames(
            title="Select Input File(s)",
            initialdir=script_dir,
            filetypes=(("All Supported Files", "*.nef;*.txt;*.str;*.star;*.pdb;*.tab;*.talos;*.list;*.sparky;*.tbl"),
                       ("NEF Files", "*.nef"),
                       ("NMR-STAR Files", "*.txt;*.str;*.star"),
                       ("PDB Files", "*.pdb"),
                       ("TALOS Files", "*.tab;*.talos"),
                       ("Sparky Files", "*.list;*.sparky"),
                       ("Yasara TBL Files", "*.tbl"),
                       ("All Files", "*.*"))
        )
        if filenames:
            self.input_files = list(filenames)
            self.input_entry.delete(0, tk.END)
            
            if len(self.input_files) == 1:
                self.input_entry.insert(0, self.input_files[0])
                # Suggest output filename
                dir_name = os.path.dirname(self.input_files[0])
                base_name = os.path.basename(self.input_files[0])
                name, ext = os.path.splitext(base_name)
                self.output_entry.delete(0, tk.END)
                self.output_entry.insert(0, os.path.join(dir_name, f"{name}_renumbered{ext}"))
            else:
                self.input_entry.insert(0, f"{len(self.input_files)} files selected")
                # Suggest output directory
                dir_name = os.path.dirname(self.input_files[0])
                self.output_entry.delete(0, tk.END)
                self.output_entry.insert(0, dir_name)
                
            # Update preview dropdown selector
            self.preview_combo.config(values=[os.path.basename(f) for f in self.input_files])
            self.preview_combo.current(0)
            self.load_preview()

    def browse_output(self):
        if len(self.input_files) <= 1:
            initial = self.output_entry.get()
            dir_init = os.path.dirname(initial) if initial else ""
            file_init = os.path.basename(initial) if initial else ""
            filename = filedialog.asksaveasfilename(
                title="Select Output File",
                initialdir=dir_init,
                initialfile=file_init,
                filetypes=(("All Files", "*.*"),)
            )
            if filename:
                self.output_entry.delete(0, tk.END)
                self.output_entry.insert(0, filename)
        else:
            # Batch mode: select directory
            dirname = filedialog.askdirectory(title="Select Output Directory", initialdir=self.output_entry.get())
            if dirname:
                self.output_entry.delete(0, tk.END)
                self.output_entry.insert(0, dirname)

    def get_offset(self):
        try:
            return int(self.offset_entry.get())
        except ValueError:
            messagebox.showerror("Invalid Input", "Please enter a valid integer for Residue Offset.")
            return None

    def get_float_offset(self, entry, name):
        try:
            return float(entry.get())
        except ValueError:
            messagebox.showerror("Invalid Input", f"Please enter a valid number for {name}.")
            return None

    def get_active_format(self, filepath):
        if self.autodetect_var.get():
            return detect_file_format(filepath)
        else:
            mapping = {
                "NEF": "nef",
                "NMR-STAR (txt)": "txt",
                "Yasara TBL": "fallback",
                "PDB": "pdb",
                "TALOS": "talos",
                "Sparky": "sparky"
            }
            return mapping.get(self.force_format_combo.get(), "fallback")

    def load_preview(self):
        if not self.input_files:
            return
            
        # Identify which file in the batch we are previewing
        active_idx = self.preview_combo.current()
        if active_idx == -1:
            active_idx = 0
        filepath = self.input_files[active_idx]
        
        if not os.path.exists(filepath):
            messagebox.showwarning("File Error", f"File does not exist: {filepath}")
            return
            
        residue_offset = self.get_offset()
        if residue_offset is None: return
        proton_offset = self.get_float_offset(self.offset_h_entry, "Proton Offset")
        if proton_offset is None: return
        carbon_offset = self.get_float_offset(self.offset_c_entry, "Carbon Offset")
        if carbon_offset is None: return
        nitrogen_offset = self.get_float_offset(self.offset_n_entry, "Nitrogen Offset")
        if nitrogen_offset is None: return
        
        try:
            ftype = self.get_active_format(filepath)
            nomenclature_mode = self.nomenclature_combo.get()
            
            with open(filepath, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                
            self.text_original.delete(1.0, tk.END)
            self.text_original.insert(tk.END, "".join(lines))
            
            warnings = []
            processed_lines = process_lines(
                lines, 
                residue_offset, 
                proton_offset, 
                carbon_offset, 
                nitrogen_offset, 
                ftype, 
                nomenclature_mode, 
                warnings
            )
            
            self.text_preview.delete(1.0, tk.END)
            self.text_preview.insert(tk.END, "".join(processed_lines))
            
            # Populate validation log
            self.text_log.delete(1.0, tk.END)
            if warnings:
                self.text_log.insert(tk.END, f"Found {len(warnings)} validation warnings in {os.path.basename(filepath)}:\n")
                for w in warnings:
                    self.text_log.insert(tk.END, f" - {w}\n")
            else:
                self.text_log.insert(tk.END, "Validation complete: No chemical shift out-of-bounds alerts found.\n")
                
            self.status_label.config(text=f"Loaded preview for {os.path.basename(filepath)}", fg="blue")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load preview: {str(e)}")

    def process_files(self):
        if not self.input_files:
            messagebox.showwarning("File Error", "Please select input file(s).")
            return
            
        residue_offset = self.get_offset()
        if residue_offset is None: return
        proton_offset = self.get_float_offset(self.offset_h_entry, "Proton Offset")
        if proton_offset is None: return
        carbon_offset = self.get_float_offset(self.offset_c_entry, "Carbon Offset")
        if carbon_offset is None: return
        nitrogen_offset = self.get_float_offset(self.offset_n_entry, "Nitrogen Offset")
        if nitrogen_offset is None: return
        
        output_target = self.output_entry.get().strip()
        if not output_target:
            messagebox.showwarning("Output Error", "Please specify an output location.")
            return
            
        nomenclature_mode = self.nomenclature_combo.get()
        self.status_label.config(text="Processing...", fg="orange")
        self.root.update()
        
        self.text_log.delete(1.0, tk.END)
        all_warnings = []
        success_count = 0
        
        try:
            # Determine output paths
            output_paths = []
            if len(self.input_files) == 1:
                output_paths = [output_target]
            else:
                # Batch mode: output_target must be a directory
                if not os.path.isdir(output_target):
                    try:
                        os.makedirs(output_target, exist_ok=True)
                    except Exception as ex:
                        messagebox.showerror("Directory Error", f"Could not create output directory: {str(ex)}")
                        return
                for f in self.input_files:
                    base = os.path.basename(f)
                    name, ext = os.path.splitext(base)
                    output_paths.append(os.path.join(output_target, f"{name}_renumbered{ext}"))
            
            for inf, outf in zip(self.input_files, output_paths):
                ftype = self.get_active_format(inf)
                
                with open(inf, 'r', encoding='utf-8') as f:
                    lines = f.readlines()
                    
                warnings = []
                processed_lines = process_lines(
                    lines, 
                    residue_offset, 
                    proton_offset, 
                    carbon_offset, 
                    nitrogen_offset, 
                    ftype, 
                    nomenclature_mode, 
                    warnings
                )
                
                if warnings:
                    all_warnings.append((os.path.basename(inf), warnings))
                    
                with open(outf, 'w', encoding='utf-8') as f:
                    f.writelines(processed_lines)
                success_count += 1
                
            # Log results
            self.text_log.insert(tk.END, f"Processed {success_count} file(s) successfully.\n")
            if all_warnings:
                self.text_log.insert(tk.END, "\nValidation Alerts Log:\n")
                for fname, warnings in all_warnings:
                    self.text_log.insert(tk.END, f"[{fname}]\n")
                    for w in warnings:
                        self.text_log.insert(tk.END, f"  - {w}\n")
                        
            self.status_label.config(text=f"Processed {success_count} file(s) successfully.", fg="green")
            messagebox.showinfo("Success", f"Successfully processed {success_count} files.")
            
        except Exception as e:
            self.status_label.config(text="Processing failed.", fg="red")
            messagebox.showerror("Processing Error", f"An error occurred: {str(e)}")


def run_cli():
    parser = argparse.ArgumentParser(description="Residue Renumberer & Referencing Correction Tool (CLI Mode)")
    parser.add_argument("-i", "--input", nargs="+", required=True, help="Input file path(s)")
    parser.add_argument("-o", "--output", nargs="+", help="Output file/dir path(s)")
    parser.add_argument("-r", "--residue", type=int, default=0, help="Residue offset (integer)")
    parser.add_argument("--proton", type=float, default=0.0, help="Proton (1H) chemical shift offset (ppm)")
    parser.add_argument("--carbon", type=float, default=0.0, help="Carbon (13C) chemical shift offset (ppm)")
    parser.add_argument("--nitrogen", type=float, default=0.0, help="Nitrogen (15N) chemical shift offset (ppm)")
    parser.add_argument("--no-detect", action="store_true", help="Disable automatic format detection and rely on extensions")
    parser.add_argument("-n", "--nomenclature", choices=["NONE", "IUPAC_TO_CYANA", "CYANA_TO_IUPAC"], default="NONE", help="Nomenclature mapping standard")
    
    args = parser.parse_args()
    
    input_files = args.input
    output_files = args.output
    
    if output_files:
        if len(output_files) == 1 and (os.path.isdir(output_files[0]) or not os.path.splitext(output_files[0])[1]):
            # Output is a directory
            out_dir = output_files[0]
            os.makedirs(out_dir, exist_ok=True)
            output_paths = []
            for f in input_files:
                base = os.path.basename(f)
                name, ext = os.path.splitext(base)
                output_paths.append(os.path.join(out_dir, f"{name}_renumbered{ext}"))
        elif len(output_files) == len(input_files):
            output_paths = output_files
        else:
            print("Error: The number of output paths must match the number of input paths, or be a single directory.", file=sys.stderr)
            sys.exit(1)
    else:
        # Default suffixing
        output_paths = []
        for f in input_files:
            dir_name = os.path.dirname(f)
            base = os.path.basename(f)
            name, ext = os.path.splitext(base)
            output_paths.append(os.path.join(dir_name, f"{name}_renumbered{ext}"))
            
    for inf, outf in zip(input_files, output_paths):
        if not os.path.exists(inf):
            print(f"Error: File not found: {inf}", file=sys.stderr)
            continue
            
        if args.no_detect:
            _, ext = os.path.splitext(inf.lower())
            if ext == '.nef': ftype = 'nef'
            elif ext in ('.txt', '.str', '.star'): ftype = 'txt'
            elif ext == '.pdb': ftype = 'pdb'
            elif ext in ('.tab', '.talos'): ftype = 'talos'
            elif ext in ('.list', '.sparky'): ftype = 'sparky'
            else: ftype = 'fallback'
        else:
            ftype = detect_file_format(inf)
            
        print(f"\nProcessing: {inf} -> {outf} (Format: {ftype.upper()})")
        
        try:
            with open(inf, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                
            warnings = []
            processed = process_lines(
                lines, 
                args.residue, 
                args.proton, 
                args.carbon, 
                args.nitrogen, 
                ftype, 
                args.nomenclature, 
                warnings
            )
            
            if warnings:
                print("Validation Alerts:")
                for w in warnings:
                    print(f"  [WARNING] {w}")
                    
            with open(outf, 'w', encoding='utf-8') as f:
                f.writelines(processed)
                
            print("Completed successfully.")
            
        except Exception as e:
            print(f"Error processing {inf}: {str(e)}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        run_cli()
    else:
        root = tk.Tk()
        app = RenumberApp(root)
        root.mainloop()
