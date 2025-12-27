import re
import streamlit as st
from urllib.parse import unquote

from molecule_helpers import (
    _maybe_imports,
    handle_polymer_name,
    get_compound_from_name,
    get_compound_from_iupac_web,
    get_chemical_properties,
    generate_3d_mol,
    smiles_to_2d_image,
    add_functional_group_labels,
)

st.set_page_config(page_title="TwinSpec • Structure", layout="wide", initial_sidebar_state="collapsed")

# ----------------------------
# Helpers
# ----------------------------

# Keep curated keys LOWERCASE only; normalize user input before lookup.
# Commodity polymer shorthands map to canonical names (resolved via PubChem/OPSIN).
# Conjugated polymer shorthands map to a curated SMILES proxy.
CURATED_SMILES = {
    # simple abbreviations -> canonical names
    "pe": "polyethylene",
    "pp": "polypropylene",
    "ps": "polystyrene",
    "pmma": "poly(methyl methacrylate)",
    "pet": "poly(ethylene terephthalate)",

    # conjugated polymer shorthand -> curated SMILES proxy
    "p(ndi2od-t2)": "CCCCCCCCCCC(CCCCCCCC)Cn3c(=O)c4ccc5c(=O)n(CC(CCCCCCC)CCCCCCCCCC)c(=O)c6c(c2ccc(c1cccs1)s2)cc(c3=O)c4c56",
    "pndi2odt2":    "CCCCCCCCCCC(CCCCCCCC)Cn3c(=O)c4ccc5c(=O)n(CC(CCCCCCC)CCCCCCCCCC)c(=O)c6c(c2ccc(c1cccs1)s2)cc(c3=O)c4c56",
    "ndi2od-t2":    "CCCCCCCCCCC(CCCCCCCC)Cn3c(=O)c4ccc5c(=O)n(CC(CCCCCCC)CCCCCCCCCC)c(=O)c6c(c2ccc(c1cccs1)s2)cc(c3=O)c4c56",
    "ndi2odt2":     "CCCCCCCCCCC(CCCCCCCC)Cn3c(=O)c4ccc5c(=O)n(CC(CCCCCCC)CCCCCCCCCC)c(=O)c6c(c2ccc(c1cccs1)s2)cc(c3=O)c4c56",
    "n2200":        "CCCCCCCCCCC(CCCCCCCC)Cn3c(=O)c4ccc5c(=O)n(CC(CCCCCCC)CCCCCCCCCC)c(=O)c6c(c2ccc(c1cccs1)s2)cc(c3=O)c4c56",
    "pndi2od-t2":   "CCCCCCCCCCC(CCCCCCCC)Cn3c(=O)c4ccc5c(=O)n(CC(CCCCCCC)CCCCCCCCCC)c(=O)c6c(c2ccc(c1cccs1)s2)cc(c3=O)c4c56",
}

def normalize_user_input(raw: str) -> str:
    s = (raw or "").strip()
    if not s:
        return ""

    # normalize whitespace + unicode-ish dashes
    s = s.replace("–", "-").replace("—", "-")
    s = re.sub(r"\s+", " ", s).strip()

    # common polymer wrappers:
    # poly(X) -> X, poly X -> X, poly-X -> X
    s_l = s.lower()
    if s_l.startswith("poly(") and s.endswith(")"):
        s = s[5:-1].strip()
    elif s_l.startswith("poly "):
        s = s[5:].strip()
    elif s_l.startswith("poly-"):
        s = s[5:].strip()

    # collapse some punctuation that breaks lookups
    s = s.strip().strip('"').strip("'")
    return s

def looks_like_smiles(s: str) -> bool:
    if " " in s:
        return False
    return bool(re.search(r"[\[\]\(\)=#@\\/%]|Br|Cl|Si|Na|Li|Mg|Al|Ca", s))

def _curated_lookup(raw_norm: str):
    """
    Returns (curated_value, matched_key) or (None, None)
    Tries:
      - lowercase exact
      - whitespace-stripped
      - if wrapped as p(...), try inner as well
      - also try re-wrapping inner as p(inner)
    """
    low = (raw_norm or "").strip().lower()
    if not low:
        return None, None

    key = low.replace(" ", "")

    keys_to_try = [low, key]

    # If user typed P(...), try inner
    if key.startswith("p(") and key.endswith(")"):
        inner = key[2:-1]
        keys_to_try.append(inner)
        # also try normalized wrapped form (in case user had spaces)
        keys_to_try.append(f"p({inner})")
    else:
        # If user typed inner, also try wrapped key
        keys_to_try.append(f"p({key})")

    for k in keys_to_try:
        v = CURATED_SMILES.get(k)
        if v:
            return v, k
    return None, None

def resolve_to_smiles(user_input: str):
    """
    Returns:
      (resolved_label, smiles, polymer_info)
    where polymer_info is either None or dict with keys:
      { "monomer_name", "monomer_smiles", "repeat_unit", "source" }
    """
    raw = (user_input or "").strip()
    norm = normalize_user_input(raw)
    if not norm:
        return ("unresolved", "", None)

    # 1) curated alias
    curated, matched_key = _curated_lookup(norm)
    if curated:
        st.caption(f"Resolved via curated alias: `{matched_key}`")

        # If curated maps to a name (e.g., 'polyethylene'), treat as name and keep going
        if not looks_like_smiles(curated):
            norm = curated
        else:
            # curated SMILES proxy: treat as polymer-ish but provide required fields (no KeyError)
            polymer_info = {
                "source": "curated_map",
                "monomer_name": matched_key,   # display-only
                "monomer_smiles": curated,
                "repeat_unit": None,
            }
            return (f"curated:{matched_key}", curated, polymer_info)

    # 2) polymer handler (your existing logic; works for poly(NAME) cases)
    monomer_name, monomer_compound, repeat_unit = handle_polymer_name(norm)
    if monomer_compound:
        mon_smiles = monomer_compound.smiles
        return (
            f"polymer:{norm}",
            mon_smiles,
            {
                "source": "polymer_name",
                "monomer_name": monomer_name,
                "monomer_smiles": mon_smiles,
                "repeat_unit": repeat_unit,
            },
        )

    # 3) name resolution (PubChem / OPSIN)
    compound = get_compound_from_name(norm) or get_compound_from_iupac_web(norm)
    if compound and getattr(compound, "smiles", None):
        return (f"name:{norm}", compound.smiles, None)

    # 4) direct SMILES input fallback
    if looks_like_smiles(norm):
        return ("smiles", norm, None)

    # 5) give up cleanly (but show what we tried)
    return ("unresolved", norm, None)

def _render_3d(mods, smiles: str, width: int, height: int):
    rd_mol = generate_3d_mol(smiles)
    if not rd_mol:
        st.error("Could not generate 3D geometry (invalid SMILES or RDKit issue).")
        return

    mol_block = mods["Chem"].MolToMolBlock(rd_mol)

    view = mods["py3Dmol"].view(width=width, height=height)
    view.addModel(mol_block, "mol")
    view.setStyle({"stick": {}})
    view.setBackgroundColor("0xeeeeee")

    # Add labels BEFORE final camera framing (labels can change bounds)
    view = add_functional_group_labels(view, rd_mol)

    # Final camera + render LAST
    view.zoomTo()
    view.render()

    mods["showmol"](view, height=height, width=width)

    # mol_block = mods["Chem"].MolToMolBlock(rd_mol)
    # view = mods["py3Dmol"].view(width=width, height=height)
    # view.addModel(mol_block, "mol")
    # view.setStyle({"stick": {}})
    # view.setBackgroundColor("0xeeeeee")
    # view.zoomTo()
    # view = add_functional_group_labels(view, rd_mol)
    # # labels can change bounding box slightly; re-center
    # view.zoomTo()
    # mods["showmol"](view, height=height, width=width)

# ----------------------------
# UI
# ----------------------------

# Pull query param
qp = st.query_params
default = qp.get("smiles", "") if isinstance(qp.get("smiles", ""), str) else ""
default = unquote(default) if default else ""

# Embed mode for iframe usage (Mode B)
embed = str(qp.get("embed", "0")).lower() in ("1", "true", "yes", "y")

if embed:
    st.markdown(
        """
        <style>
          #MainMenu {visibility: hidden;}
          header {visibility: hidden;}
          footer {visibility: hidden;}
          .block-container {padding-top: 0.6rem; padding-bottom: 0.6rem;}
        </style>
        """,
        unsafe_allow_html=True,
    )

if not embed:
    st.markdown("### TwinSpec • Chemical structure (placeholder)")

chem = st.text_input("SMILES / name", value=default)

mods = _maybe_imports()
missing = [k for k, v in {
    "RDKit": mods["Chem"],
    "py3Dmol": mods["py3Dmol"],
    "stmol": mods["showmol"],
    "PubChemPy": mods["pcp"],
}.items() if v is None]

if missing:
    st.warning("Missing packages: " + ", ".join(missing))
    st.stop()

if not chem.strip():
    st.info("Enter a SMILES or chemical name.")
    st.stop()

resolved_label, smiles, polymer_info = resolve_to_smiles(chem)

# Tune viewer size for embed vs full page
VIEW_W = 380 if embed else 900
VIEW_H = 380 if embed else 560

# ----------------------------
# Polymer path
# ----------------------------
if polymer_info is not None:
    mon_name = polymer_info.get("monomer_name") or "monomer (proxy)"
    mon_smiles = polymer_info.get("monomer_smiles") or smiles
    repeat_unit = polymer_info.get("repeat_unit")

    if not embed:
        st.info(f"Detected polymer. Showing monomer **{mon_name}** and representative repeat unit.")

        c1, c2 = st.columns(2)
        with c1:
            st.subheader("Monomer (2D)")
            img = smiles_to_2d_image(mon_smiles, legend=f"{mon_name} (Monomer)")
            if img:
                st.image(img, width=400)
            st.code(mon_smiles)

        with c2:
            st.subheader("Repeat unit (2D)")
            if repeat_unit:
                img2 = smiles_to_2d_image(repeat_unit, legend=f"{chem} (Repeat Unit)")
                if img2:
                    st.image(img2, width=400)
                st.code(repeat_unit)
            else:
                st.caption("No repeat unit provided for this curated entry.")
                st.write("—")

        st.subheader("3D (Monomer)")
        _render_3d(mods, mon_smiles, width=VIEW_W, height=VIEW_H)

    else:
        # EMBED MODE:
        # 1) show 3D first so it's immediately visible inside the iframe
        st.subheader("3D (Monomer)")
        _render_3d(mods, mon_smiles, width=VIEW_W, height=VIEW_H)

        # 2) tuck the rest under an expander to avoid iframe scrolling
        with st.expander("2D + repeat unit", expanded=False):
            c1, c2 = st.columns(2)
            with c1:
                st.subheader("Monomer (2D)")
                img = smiles_to_2d_image(mon_smiles, legend=f"{mon_name} (Monomer)")
                if img:
                    st.image(img, use_container_width=True)
                st.code(mon_smiles)

            with c2:
                st.subheader("Repeat unit (2D)")
                if repeat_unit:
                    img2 = smiles_to_2d_image(repeat_unit, legend=f"{chem} (Repeat Unit)")
                    if img2:
                        st.image(img2, use_container_width=True)
                    st.code(repeat_unit)
                else:
                    st.caption("No repeat unit provided for this curated entry.")
                    st.write("—")

# ----------------------------
# Non-polymer path
# ----------------------------
else:
    if not embed and resolved_label == "unresolved":
        st.warning(
            "Could not resolve this as a known name or SMILES. "
            "If this is a polymer abbreviation (e.g., P(NDI2OD-T2)), add a curated alias or provide SMILES."
        )

    if not embed:
        st.caption(f"Using SMILES: `{smiles}`" if resolved_label != "unresolved" else f"Using input: `{smiles}`")

        c1, c2 = st.columns(2)
        with c1:
            st.subheader("Properties")
            props = get_chemical_properties(smiles)
            if props:
                for k, v in props.items():
                    st.markdown(f"**{k}:** {v}")
            else:
                st.write("—")

        with c2:
            st.subheader("2D")
            img = smiles_to_2d_image(smiles)
            if img:
                st.image(img, width=400)

        st.subheader("3D")
        _render_3d(mods, smiles, width=VIEW_W, height=VIEW_H)

    else:
        # EMBED MODE:
        # 1) show 3D first
        st.subheader("3D")
        _render_3d(mods, smiles, width=VIEW_W, height=VIEW_H)

        # 2) tuck properties + 2D away
        with st.expander("2D + properties", expanded=False):
            c1, c2 = st.columns(2)
            with c1:
                st.subheader("Properties")
                props = get_chemical_properties(smiles)
                if props:
                    for k, v in props.items():
                        st.markdown(f"**{k}:** {v}")
                else:
                    st.write("—")

            with c2:
                st.subheader("2D")
                img = smiles_to_2d_image(smiles)
                if img:
                    st.image(img, use_container_width=True)