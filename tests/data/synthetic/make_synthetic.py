"""Write ``SYNPHOS.pdb``, the synthetic structure of the modified-residue tests.

Run from this directory: ``python make_synthetic.py``.  The output is fully
determined by this script (no randomness, no network), so the committed file
can be regenerated and compared byte for byte.

The geometry is not a real protein: residues are laid out along two parallel
paths 4.6 A apart and each atom is placed at a fixed offset from its residue's
centre, with a small deterministic perturbation so no three atoms are
collinear.  What the file is for is its *records*: phosphothreonine (TPO),
phosphoserine (SEP) and phosphotyrosine (PTR) written as HETATM inside chain
A, as RCSB files write them; a modified cysteine (CME) that BioTSR does not
read; two residues with insertion codes (10A, 10B); a second chain (B) with
its own TPO; and, after TER, a small ligand (LIG) and a calcium ion whose atom
is also named CA.
"""

import hashlib
import math
import os

ATOMS = {
    "GLY": "N CA C O",
    "ALA": "N CA C O CB",
    "SER": "N CA C O CB OG",
    "THR": "N CA C O CB OG1 CG2",
    "VAL": "N CA C O CB CG1 CG2",
    "LEU": "N CA C O CB CG CD1 CD2",
    "LYS": "N CA C O CB CG CD CE NZ",
    "ASP": "N CA C O CB CG OD1 OD2",
    "TYR": "N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH",
    "TPO": "N CA CB CG2 OG1 P O1P O2P O3P C O",
    "SEP": "N CA CB OG P O1P O2P O3P C O",
    "PTR": "N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH P O1P O2P O3P",
    "CME": "N CA CB SG SD CE CZ OH C O",
    "LIG": "C1 C2 C3 C4 C5 C6 N1 O1",
}
HETERO = {"TPO", "SEP", "PTR", "CME", "LIG", "CA"}

CHAIN_A = [("1", "GLY"), ("2", "ALA"), ("3", "TPO"), ("4", "SER"), ("5", "SEP"),
           ("6", "TYR"), ("7", "PTR"), ("8", "CME"), ("9", "LEU"), ("10", "LYS"),
           ("10A", "GLY"), ("10B", "ALA"), ("11", "THR"), ("12", "VAL")]
CHAIN_B = [("1", "GLY"), ("2", "TYR"), ("3", "TPO"), ("4", "LYS"), ("5", "ASP"),
           ("6", "LEU")]


def _jitter(*key):
    digest = hashlib.sha256(":".join(key).encode()).digest()
    return [(b / 255.0 - 0.5) * 0.6 for b in digest[:3]]


def _offset(index):
    angle = 2.0 * math.pi * index / 7.0
    radius = 1.2 + 0.35 * (index % 4)
    return (0.45 * index - 1.5, radius * math.cos(angle), radius * math.sin(angle))


def _line(record, serial, name, resname, chain, resid, xyz, element):
    seq = "".join(c for c in resid if c.isdigit() or c == "-")
    icode = resid[len(seq):]
    field = (" %-3s" % name) if len(element) == 1 and len(name) < 4 else "%-4s" % name
    return ("%-6s%5d %-4s %3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n"
            % (record, serial, field, resname, chain, int(seq), icode,
               xyz[0], xyz[1], xyz[2], 1.0, 20.0, element))


def build():
    lines = ["HEADER    SYNTHETIC TEST STRUCTURE                 01-OCT-26   SYNP\n",
             "REMARK   1 WRITTEN BY tests/data/synthetic/make_synthetic.py; NOT A\n",
             "REMARK   1 REAL STRUCTURE. GEOMETRY IS ARBITRARY.\n"]
    serial = 0
    for chain, residues, y0 in (("A", CHAIN_A, 0.0), ("B", CHAIN_B, 4.6)):
        for position, (resid, resname) in enumerate(residues):
            centre = (3.8 * position, y0 + 0.8 * math.sin(position), 0.8 * math.cos(position))
            for index, name in enumerate(ATOMS[resname].split()):
                dx, dy, dz = _offset(index)
                jx, jy, jz = _jitter(chain, resid, name)
                xyz = (centre[0] + dx + jx, centre[1] + dy + jy, centre[2] + dz + jz)
                serial += 1
                record = "HETATM" if resname in HETERO else "ATOM"
                lines.append(_line(record, serial, name, resname, chain, resid, xyz,
                                   name[0]))
        serial += 1
        last_id, last_name = residues[-1]
        seq = "".join(c for c in last_id if c.isdigit())
        lines.append("TER   %5d      %3s %1s%4d\n" % (serial, last_name, chain, int(seq)))
    for index, name in enumerate(ATOMS["LIG"].split()):
        dx, dy, dz = _offset(index)
        jx, jy, jz = _jitter("A", "101", name)
        serial += 1
        lines.append(_line("HETATM", serial, name, "LIG", "A", "101",
                           (11.0 + dx + jx, -3.5 + dy + jy, 2.0 + dz + jz), name[0]))
    serial += 1
    lines.append(_line("HETATM", serial, "CA", "CA", "A", "102", (20.0, -3.0, 0.0), "CA"))
    lines.append("END\n")
    return "".join(lines)


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "SYNPHOS.pdb"), "w", newline="\n") as fh:
        fh.write(build())
