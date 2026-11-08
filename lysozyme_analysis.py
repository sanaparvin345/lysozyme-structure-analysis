"""
Structural analysis of hen egg-white lysozyme (PDB: 1LYZ) with Biopython.

What it does
  1. Downloads the structure from the RCSB PDB
  2. Calculates solvent accessibility (SASA) for every residue
  3. Counts residue contacts (CA atoms within 8 A) and draws a contact map
  4. Draws a Ramachandran plot (phi/psi angles)
  5. Finds residues around the catalytic site (Glu35 and Asp52)
  6. Detects disulfide bonds (Cys SG-SG < 2.5 A)

Requirements:  pip install biopython numpy pandas matplotlib
Run:           python lysozyme_analysis.py
"""

import os
import urllib.request

import numpy as np 
import pandas as pd
import matplotlib

matplotlib.use("Agg")  # save figures to files, no pop-up windows
import matplotlib.pyplot as plt

from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB.SASA import ShrakeRupley
import math
import pandas as pd

rows = []
for pp in ppb.build_peptides(structure):
    for res, (phi, psi) in zip(pp, pp.get_phi_psi_list()):
        if phi is not None and psi is not None:
            rows.append({
                "resnum": res.id[1],
                "resname": res.get_resname(),
                "phi": round(math.degrees(phi), 1),
                "psi": round(math.degrees(psi), 1),
            })
rama = pd.DataFrame(rows)
rama.to_csv("results/phi_psi.csv", index=False)
nongly = rama[rama.resname != "GLY"]
print(nongly[(nongly.phi.between(-105, -90)) & (nongly.psi.between(-150, -135))])
print(nongly[nongly.phi > 0])
# ----------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------
PDB_ID = "1LYZ"
CATALYTIC_RESIDUES = [35, 52]   # Glu35 and Asp52 in hen lysozyme
CONTACT_CUTOFF = 8.0            # Angstrom, CA-CA distance
SITE_CUTOFF = 5.0               # Angstrom, around catalytic residues
EXPOSED_CUTOFF = 0.25           # relative SASA >= 0.25 -> "exposed"

DATA_DIR = "data"
OUT_DIR = "results"
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

# Theoretical maximum accessible surface area per residue type (A^2),
# from Tien et al. (2013), used to convert SASA into relative SASA.
MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}

# ----------------------------------------------------------------------
# 1. Download and load the structure
# ----------------------------------------------------------------------
pdb_path = os.path.join(DATA_DIR, f"{PDB_ID}.pdb")
if not os.path.exists(pdb_path):
    url = f"https://files.rcsb.org/download/{PDB_ID}.pdb"
    print(f"Downloading {url} ...")
    urllib.request.urlretrieve(url, pdb_path)

parser = PDBParser(QUIET=True)
structure = parser.get_structure(PDB_ID, pdb_path)
model = structure[0]
chain = next(model.get_chains())

# Remove water and ligands so they don't affect the surface calculation
for res in list(chain):
    if res.id[0] != " ":
        chain.detach_child(res.id)

residues = list(chain)
print(f"Loaded {PDB_ID}: chain {chain.id}, {len(residues)} residues")

# ----------------------------------------------------------------------
# 2. Solvent accessibility
# ----------------------------------------------------------------------
sr = ShrakeRupley()
sr.compute(model, level="R")  # stores the result in residue.sasa

rows = []
for r in residues:
    name = r.get_resname()
    rows.append({
        "resnum": r.id[1],
        "resname": name,
        "aa": protein_letters_3to1[name.capitalize()],
        "sasa": round(r.sasa, 2),
        "rsasa": round(min(r.sasa / MAX_ASA[name], 1.0), 3),
    })
df = pd.DataFrame(rows)
df["location"] = np.where(df["rsasa"] >= EXPOSED_CUTOFF, "exposed", "buried")

# ----------------------------------------------------------------------
# 3. Residue contacts and contact map
# ----------------------------------------------------------------------
ca = np.array([r["CA"].coord for r in residues])
dist = np.linalg.norm(ca[:, None, :] - ca[None, :, :], axis=-1)
df["contacts"] = (dist < CONTACT_CUTOFF).sum(axis=1) - 1  # minus itself

fig, ax = plt.subplots(figsize=(6, 5))
im = ax.imshow(dist, cmap="viridis_r", origin="lower")
ax.set_xlabel("Residue number")
ax.set_ylabel("Residue number")
ax.set_title(f"{PDB_ID} CA-CA distance map")
fig.colorbar(im, label="Distance (A)")
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "contact_map.png"), dpi=300)
plt.close(fig)

# ----------------------------------------------------------------------
# 4. Ramachandran plot
# ----------------------------------------------------------------------
angles = []
for pp in PPBuilder().build_peptides(chain):
    for res, (phi, psi) in zip(pp, pp.get_phi_psi_list()):
        if phi is not None and psi is not None:
            angles.append({
                "resnum": res.id[1],
                "resname": res.get_resname(),
                "phi": np.degrees(phi),
                "psi": np.degrees(psi),
            })
rama = pd.DataFrame(angles)
gly = rama[rama["resname"] == "GLY"]
other = rama[rama["resname"] != "GLY"]

fig, ax = plt.subplots(figsize=(5.5, 5.5))
ax.scatter(other["phi"], other["psi"], s=18, alpha=0.7, label="Non-Gly")
ax.scatter(gly["phi"], gly["psi"], s=30, marker="^", color="tab:red", label="Gly")
ax.set_xlim(-180, 180)
ax.set_ylim(-180, 180)
ax.axhline(0, color="grey", lw=0.5)
ax.axvline(0, color="grey", lw=0.5)
ax.set_xlabel("Phi (degrees)")
ax.set_ylabel("Psi (degrees)")
ax.set_title(f"{PDB_ID} Ramachandran plot")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "ramachandran.png"), dpi=300)
plt.close(fig)

# ----------------------------------------------------------------------
# 5. Residues around the catalytic site
# ----------------------------------------------------------------------
ns = NeighborSearch(list(chain.get_atoms()))
site = set(CATALYTIC_RESIDUES)
for resnum in CATALYTIC_RESIDUES:
    for atom in chain[resnum]:
        for neighbour in ns.search(atom.coord, SITE_CUTOFF, level="R"):
            site.add(neighbour.id[1])

df["catalytic_site"] = df["resnum"].isin(site)
site_df = df[df["catalytic_site"]]
site_df.to_csv(os.path.join(OUT_DIR, "catalytic_site_residues.csv"), index=False)

# ----------------------------------------------------------------------
# 6. Disulfide bonds
# ----------------------------------------------------------------------
cys = [r for r in residues if r.get_resname() == "CYS"]
bonds = []
for i in range(len(cys)):
    for j in range(i + 1, len(cys)):
        d = cys[i]["SG"] - cys[j]["SG"]
        if d < 2.5:
            bonds.append({"cys_1": cys[i].id[1], "cys_2": cys[j].id[1],
                          "SG_SG_distance": round(float(d), 2)})
ss_df = pd.DataFrame(bonds)
ss_df.to_csv(os.path.join(OUT_DIR, "disulfide_bonds.csv"), index=False)

# ----------------------------------------------------------------------
# 7. Relative SASA along the sequence
# ----------------------------------------------------------------------
colors = np.where(df["location"] == "exposed", "tab:orange", "tab:blue")
fig, ax = plt.subplots(figsize=(10, 4))
ax.bar(df["resnum"], df["rsasa"], color=colors, width=0.9)
ax.axhline(EXPOSED_CUTOFF, color="black", ls="--", lw=0.8)
for resnum in CATALYTIC_RESIDUES:
    ax.axvline(resnum, color="red", lw=1)
    ax.text(resnum, 1.02, str(resnum), color="red", ha="center", fontsize=8)
ax.set_ylim(0, 1.1)
ax.set_xlabel("Residue number")
ax.set_ylabel("Relative SASA")
ax.set_title(f"{PDB_ID} solvent accessibility (orange = exposed, blue = buried)")
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "rsasa_profile.png"), dpi=300)
plt.close(fig)

# ----------------------------------------------------------------------
# Save table and print summary
# ----------------------------------------------------------------------
df.to_csv(os.path.join(OUT_DIR, "residue_table.csv"), index=False)

corr = df["rsasa"].corr(df["contacts"])
print("\n--- Summary ---")
print(f"Residues analysed        : {len(df)}")
print(f"Exposed / buried         : {(df.location == 'exposed').sum()} / "
      f"{(df.location == 'buried').sum()}")
print(f"Correlation rSASA vs contacts: {corr:.2f}")
print(f"Catalytic-site residues  : {sorted(site)}")
print(f"Disulfide bonds found    : {len(ss_df)}")
if len(ss_df):
    print(ss_df.to_string(index=False))
print(f"\nResults saved in '{OUT_DIR}/'")
