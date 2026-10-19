from pymol import cmd, stored
import numpy as np
from scipy.spatial.distance import pdist, squareform

def AFcurate(
    object_name="all",
    min_plddt=70.0,
    max_lur_frac=0.30,
    max_loop_frac=0.65,
    min_len=40,
    min_sse=3,
    min_contact_density=1.8
):
    """
    Curates AlphaFold structures, leaving well-folded/passing regions untouched,
    and only coloring/selecting regions that fail quality criteria.
    
    Usage in PyMOL:
        AFcurate [object_name]   <-- Press TAB to auto-complete loaded objects!
        AFcurate my_protein, min_plddt=70.0
    """
    min_plddt = float(min_plddt)
    max_lur_frac = float(max_lur_frac)
    max_loop_frac = float(max_loop_frac)
    min_len = int(min_len)
    min_sse = int(min_sse)
    min_contact_density = float(min_contact_density)

    # 1. Update PyMOL secondary structure assignments
    cmd.dss(object_name)

    # 2. Extract C-alpha atom data
    stored.ca_data = []
    ca_selector = f"({object_name}) and name CA and polymer.protein"
    
    cmd.iterate(
        ca_selector,
        "stored.ca_data.append((chain, resi, resn, b, ss))"
    )

    if not stored.ca_data:
        print(f"[AF2-Curator] Error: No protein C-alpha atoms found in selection '{object_name}'.")
        return

    chains = [x[0] for x in stored.ca_data]
    resis = [x[1] for x in stored.ca_data]
    plddts = np.array([x[3] for x in stored.ca_data])
    ss_list = [x[4] for x in stored.ca_data]
    coords = cmd.get_coords(ca_selector)
    n_res = len(plddts)

    # --- Parameter 1: Low pLDDT Residues (< 70) ---
    low_plddt_mask = plddts < min_plddt
    avg_plddt = float(np.mean(plddts))
    p1_pass = avg_plddt >= min_plddt

    # --- Parameter 2: Long Unordered Regions (LURs: >= 5 continuous res with pLDDT < 70) ---
    lur_indices = []
    current_run = []
    for i, is_low in enumerate(low_plddt_mask):
        if is_low:
            current_run.append(i)
        else:
            if len(current_run) >= 5:
                lur_indices.extend(current_run)
            current_run = []
    if len(current_run) >= 5:
        lur_indices.extend(current_run)

    max_lur_len = max([len(lur_indices)] if not lur_indices else [
        len(run) for run in np.split(lur_indices, np.where(np.diff(lur_indices) != 1)[0] + 1)
    ], default=0)
    max_lur_ratio = max_lur_len / n_res
    p2_pass = max_lur_ratio <= max_lur_frac

    # --- Parameter 3 & 5: Loop Content & SSE Count ---
    loop_indices = [i for i, ss in enumerate(ss_list) if ss not in ['H', 'S']]
    loop_ratio = len(loop_indices) / n_res
    p3_pass = loop_ratio <= max_loop_frac

    # Count distinct contiguous Helices (>=3 res) and Strands (>=2 res)
    sse_count = 0
    curr_type = None
    curr_len = 0
    for ss in ss_list:
        stype = 'H' if ss == 'H' else ('S' if ss == 'S' else 'L')
        if stype == curr_type:
            curr_len += 1
        else:
            if curr_type == 'H' and curr_len >= 3:
                sse_count += 1
            elif curr_type == 'S' and curr_len >= 2:
                sse_count += 1
            curr_type = stype
            curr_len = 1
    if curr_type == 'H' and curr_len >= 3:
        sse_count += 1
    elif curr_type == 'S' and curr_len >= 2:
        sse_count += 1

    p4_pass = n_res >= min_len
    p5_pass = sse_count >= min_sse

    # --- Parameter 6: Inter-residue Contact Density & Underpacked Residues ---
    dist_matrix = squareform(pdist(coords))
    per_res_contacts = np.zeros(n_res)
    total_contacts = 0

    for i in range(n_res):
        for j in range(i + 3, n_res):
            if dist_matrix[i, j] <= 8.0:
                total_contacts += 1
                per_res_contacts[i] += 1
                per_res_contacts[j] += 1

    contact_density = total_contacts / n_res
    p6_pass = contact_density >= min_contact_density
    underpacked_indices = np.where(per_res_contacts < 2)[0]

    # --- Helper: Convert index list to PyMOL selection string ---
    def build_selection_query(indices):
        if len(indices) == 0:
            return "none"
        res_selectors = [f"(chain {chains[i]} and resi {resis[i]})" if chains[i] else f"resi {resis[i]}" for i in indices]
        return f"({object_name}) and ({' or '.join(res_selectors)})"

    # --- Create PyMOL Named Selections for Flagged Regions ---
    clean_prefix = str(object_name).replace(" ", "_").replace("(", "").replace(")", "")
    obj_prefix = clean_prefix if clean_prefix != "all" else "af2"
    
    cmd.select(f"{obj_prefix}_LURs", build_selection_query(lur_indices))
    cmd.select(f"{obj_prefix}_low_plddt", build_selection_query(np.where(low_plddt_mask)[0]))
    cmd.select(f"{obj_prefix}_underpacked", build_selection_query(underpacked_indices))

    # --- Color ONLY the Problematic Regions ---
    # 1. Underpacked residues (Yellow)
    cmd.color("yellow", f"{obj_prefix}_underpacked and not {obj_prefix}_low_plddt")
    # 2. Isolated low-pLDDT residues (Orange)
    cmd.color("orange", f"{obj_prefix}_low_plddt and not {obj_prefix}_LURs")
    # 3. Long Unordered Regions / major disorder (Dark Red)
    cmd.color("firebrick", f"{obj_prefix}_LURs")

    # --- Print Console Diagnostic & Color Legend Report ---
    passed_all = all([p1_pass, p2_pass, p3_pass, p4_pass, p5_pass, p6_pass])
    
    print("\n" + "=" * 85)
    print(f" AF2 STRUCTURE QUALITY & CURATION REPORT: {object_name}")
    print("=" * 85)
    print(f"{'Parameter':<28} | {'Value':<12} | {'Threshold':<14} | {'Status'}")
    print("-" * 85)
    print(f"{'1. Mean pLDDT':<28} | {avg_plddt:<12.2f} | >= {min_plddt:<10.1f} | {'PASS' if p1_pass else 'FAIL'}")
    print(f"{'2. Max LUR Fraction':<28} | {max_lur_ratio:<12.2%} | <= {max_lur_frac:<10.1%} | {'PASS' if p2_pass else 'FAIL'}")
    print(f"{'3. Loop Fraction':<28} | {loop_ratio:<12.2%} | <= {max_loop_frac:<10.1%} | {'PASS' if p3_pass else 'FAIL'}")
    print(f"{'4. Sequence Length':<28} | {n_res:<12} | >= {min_len:<10} | {'PASS' if p4_pass else 'FAIL'}")
    print(f"{'5. SSE Count':<28} | {sse_count:<12} | >= {min_sse:<10} | {'PASS' if p5_pass else 'FAIL'}")
    print(f"{'6. Contact Density':<28} | {contact_density:<12.2f} | >= {min_contact_density:<10.2f} | {'PASS' if p6_pass else 'FAIL'}")
    print("-" * 85)
    print(f" OVERALL CURATION DECISION: {'PASSED (Fold is reliable)' if passed_all else 'FAILED (Contains artifacts / disorder)'}")
    print("=" * 85)
    
    # Report showing ONLY flagged failure regions
    print(" FLAGGED REGIONS COLOR & SELECTION LEGEND:")
    print("-" * 85)
    print(f" {'Color':<14} | {'Selection Name':<24} | {'Condition / Flaw'}")
    print("-" * 85)
    print(f" {'🔴 Dark Red':<14} | {f'{obj_prefix}_LURs':<24} | Long Unordered Regions (>= 5 continuous res with pLDDT < {min_plddt:.0f})")
    print(f" {'🟠 Orange':<14} | {f'{obj_prefix}_low_plddt':<24} | Isolated low-confidence residues (pLDDT < {min_plddt:.0f})")
    print(f" {'🟡 Yellow':<14} | {f'{obj_prefix}_underpacked':<24} | Underpacked residues (< 2 non-local CA contacts within 8 Å)")
    print("-" * 85)
    print(" Note: All well-folded core residues retain their original colors untouched.")
    print("=" * 85 + "\n")

# Extend PyMOL CLI with the custom command
cmd.extend("AFcurate", AFcurate)

# Enable PyMOL Tab Auto-Completion for loaded objects and selections
cmd.auto_arg[0]["AFcurate"] = [cmd.selection_sc, "object or selection", ", "]