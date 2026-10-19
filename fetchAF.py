import os
import re
import json
import ssl
import colorsys
import urllib.request
import urllib.error
import tempfile
from pymol import cmd

SSL_CTX = ssl._create_unverified_context()
HTTP_HEADERS = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/122.0.0.0 Safari/537.36',
    'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,application/json,*/*;q=0.8'
}

CATH_NAME_CACHE = {}

def _generate_distinct_rgb_palette(num_colors):
    """Generates N distinct, vibrant RGB colors using Golden Ratio hue spacing."""
    palette = []
    golden_ratio_conjugate = 0.618033988749895
    h = 0.12
    for _ in range(num_colors):
        h = (h + golden_ratio_conjugate) % 1.0
        r, g, b = colorsys.hls_to_rgb(h, 0.55, 0.85)
        palette.append([round(r, 3), round(g, 3), round(b, 3)])
    return palette

def _parse_boundary_to_selection(boundary_str, object_name):
    """
    Parses single segments ('12-76') and multi-segment discontinuous domains
    ('12-50_80-120' or '12-50,80-120') into valid PyMOL selection syntax.
    """
    clean_str = str(boundary_str).replace(';', '_').replace(',', '_').replace(' ', '_')
    segments = clean_str.split('_')
    sub_queries = []
    
    for seg in segments:
        seg = seg.strip()
        if '-' in seg:
            parts = seg.split('-')
            if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
                sub_queries.append(f"resi {parts[0]}-{parts[1]}")
        elif seg.isdigit():
            sub_queries.append(f"resi {seg}")
            
    if not sub_queries:
        return "none"
    return f"({object_name}) and ({' or '.join(sub_queries)})"

def _resolve_cath_info(cath_str):
    """
    Identifies CATH hierarchy depth.
    - Curated names from REST API are displayed normally.
    - Automated external mappings from Web HTML have ' ⚠' appended.
    """
    if not cath_str or str(cath_str).strip() in ["-", "", "None", "null", "No match"]:
        return {"id": "No match", "level": "Unclassified", "name": "No match"}

    cath_clean = str(cath_str).strip()

    # Handle comma-separated consensus topology strings (e.g., '2.60.40,2.60.40')
    if ',' in cath_clean:
        items = list(dict.fromkeys([x.strip() for x in cath_clean.split(',') if x.strip()]))
        resolved_items = [_resolve_cath_info(item) for item in items]
        return {
            "id": "/".join([r["id"] for r in resolved_items]),
            "level": resolved_items[0]["level"],
            "name": " / ".join([r["name"] for r in resolved_items])
        }

    parts = cath_clean.split('.')
    num_parts = len(parts)

    levels = {
        1: "Class (C)",
        2: "Architecture (A)",
        3: "Topology (T)",
        4: "Homology (H)"
    }
    level_name = levels.get(num_parts, f"Level-{num_parts}")

    if cath_clean in CATH_NAME_CACHE:
        return CATH_NAME_CACHE[cath_clean]

    name_found = cath_clean

    if num_parts == 4:
        # 1. Check REST API for curated CATH classification name
        api_url = f"https://www.cathdb.info/version/v4_3_0/api/rest/superfamily/{cath_clean}"
        is_curated = False
        try:
            req = urllib.request.Request(api_url, headers=HTTP_HEADERS)
            with urllib.request.urlopen(req, context=SSL_CTX, timeout=8) as resp:
                data_obj = json.loads(resp.read().decode('utf-8')).get("data", {})
                if isinstance(data_obj, dict):
                    curated = data_obj.get("classification_name")
                    if curated and str(curated).strip() not in ["", "-", "None", "null"]:
                        name_found = str(curated).strip()
                        is_curated = True
        except Exception:
            pass

        # 2. If uncurated (null in API), scrape automated name from Web HTML & append ⚠
        if not is_curated:
            web_url = f"https://www.cathdb.info/version/v4_3_0/superfamily/{cath_clean}"
            try:
                req = urllib.request.Request(web_url, headers=HTTP_HEADERS)
                with urllib.request.urlopen(req, context=SSL_CTX, timeout=8) as resp:
                    html = resp.read().decode('utf-8', errors='ignore')
                    titles = re.findall(r'<h[1-3][^>]*>(.*?)</h[1-3]>', html, flags=re.IGNORECASE | re.DOTALL)
                    cleaned_titles = [re.sub(r'<[^>]+>', '', t).strip() for t in titles if t.strip()]
                    
                    # The second extracted heading is the Superfamily Title
                    if len(cleaned_titles) >= 2:
                        candidate_name = cleaned_titles[1]
                        if candidate_name and not candidate_name.startswith("CATH Superfamily") and "Functional" not in candidate_name:
                            name_found = f"{candidate_name} ⚠"
            except Exception:
                pass

    result = {"id": cath_clean, "level": level_name, "name": name_found}
    CATH_NAME_CACHE[cath_clean] = result
    return result

def _fetch_ted_domains(uniprot_id):
    """Fetches all TED domain definitions from CATH-TED API."""
    url = f"https://ted.cathdb.info/api/v1/uniprot/summary/{uniprot_id}"
    try:
        req = urllib.request.Request(url, headers=HTTP_HEADERS)
        with urllib.request.urlopen(req, context=SSL_CTX, timeout=15) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            domain_list = data.get("data", [])
            if isinstance(domain_list, list) and len(domain_list) > 0:
                parsed = []
                for idx, d in enumerate(domain_list, 1):
                    bounds = d.get("chopping") or "N/A"
                    cath_raw = d.get("cath_label")
                    cath_info = _resolve_cath_info(cath_raw)
                    plddt = d.get("plddt") or "N/A"
                    consensus = d.get("consensus_level") or "N/A"
                    
                    if bounds and str(bounds).strip() not in ["-", "", "None"]:
                        parsed.append({
                            "domain_id": idx,
                            "boundaries": str(bounds).strip(),
                            "cath_id": cath_info["id"],
                            "cath_level": cath_info["level"],
                            "cath_name": cath_info["name"],
                            "plddt": plddt,
                            "consensus": consensus
                        })
                return parsed
    except Exception as e:
        print(f"[fetchAF] Note while querying TED API: {e}")
    return []

def fetchAF(uniprot_id, color="domains", show_domains=True, save_dir=None):
    """
    Fetches AlphaFold model, annotates header with TED domains and CATH names,
    and applies domain or pLDDT coloring.
    """
    uniprot_id = str(uniprot_id).strip().upper()
    object_name = f"{uniprot_id}"
    
    print(f"\n[fetchAF] Fetching live AlphaFold & TED data for {uniprot_id}...")
    
    # 1. Fetch TED domains & CATH annotations
    ted_domains = _fetch_ted_domains(uniprot_id)
    num_domains = len(ted_domains)
    print(f"[fetchAF] Successfully retrieved {num_domains} TED domain(s) from CATH-TED API.")

    # 2. Download AlphaFold coordinates
    candidate_urls = [
        f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v6.cif",
        f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.cif",
        f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v6.pdb"
    ]
    
    raw_content = None
    for url in candidate_urls:
        try:
            req = urllib.request.Request(url, headers=HTTP_HEADERS)
            with urllib.request.urlopen(req, context=SSL_CTX, timeout=15) as resp:
                raw_content = resp.read().decode('utf-8')
                if raw_content:
                    print(f"[fetchAF] Downloaded coordinates from: {url}")
                    break
        except Exception:
            continue

    if not raw_content:
        print(f"[fetchAF] Error: Failed to download model for UniProt ID {uniprot_id}.")
        return

    # 3. Format TED Domain Header Block
    header_block = [
        "# " + "=" * 125,
        f"# ALPHAFOLD2 / TED DOMAIN ANNOTATIONS FOR {uniprot_id}",
        "# " + "=" * 125,
        f"# Total TED Domains Detected: {num_domains}",
        "# (Note: '⚠' denotes automated external mapping, uncurated in CATH)",
        "# " + "-" * 125,
        f"# {'Dom #':<7} | {'Boundaries':<20} | {'CATH Level':<16} | {'CATH Classification / Name':<45} | {'pLDDT':<8} | {'Quality'}",
        "# " + "-" * 125,
    ]
    
    for d in ted_domains:
        plddt_val = d['plddt']
        plddt_str = f"{plddt_val:.2f}" if isinstance(plddt_val, (int, float)) else str(plddt_val)
        cath_display = (d['cath_name'][:42] + '...') if len(d['cath_name']) > 45 else d['cath_name']
        header_block.append(f"# {d['domain_id']:<7} | {d['boundaries']:<20} | {d['cath_level']:<16} | {cath_display:<45} | {plddt_str:<8} | {d['consensus']}")
    header_block.append("# " + "=" * 125 + "\n")

    if raw_content.startswith("data_"):
        first_newline = raw_content.find("\n")
        annotated_content = raw_content[:first_newline+1] + "\n".join(header_block) + raw_content[first_newline+1:]
    else:
        annotated_content = "\n".join(header_block) + raw_content

    # 4. Save and Load inside PyMOL
    with tempfile.TemporaryDirectory() as tmpdirname:
        cif_file_path = os.path.join(tmpdirname, f"{uniprot_id}.cif")
        with open(cif_file_path, "w", encoding="utf-8") as f:
            f.write(annotated_content)
        
        cmd.load(cif_file_path, object_name)

        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
            perm_path = os.path.join(save_dir, f"AF-{uniprot_id}-F1-model_annotated.cif")
            with open(perm_path, "w", encoding="utf-8") as pf:
                pf.write(annotated_content)
            print(f"[fetchAF] Saved annotated CIF to: {perm_path}")

    loaded_objects = cmd.get_names("objects")
    if object_name not in loaded_objects:
        for name in loaded_objects:
            if uniprot_id in name or "tmp" in name:
                cmd.set_name(name, object_name)
                break

    # 5. Create Selections and Apply Dynamic Domain Colors
    if ted_domains and show_domains:
        palette = _generate_distinct_rgb_palette(num_domains)
        domain_selections = []

        if str(color).lower() in ["domains", "domain", "ted"]:
            cmd.color("gray80", object_name)

        print("\n" + "=" * 135)
        print(f" TED DOMAIN ARCHITECTURE: {object_name}")
        print("=" * 135)
        print(f"{'Domain Selection':<30} | {'Boundaries':<20} | {'CATH Level':<16} | {'CATH Name / ID':<48} | {'pLDDT'}")
        print("-" * 135)

        for idx, dom in enumerate(ted_domains, 1):
            bounds = dom["boundaries"]
            cath_clean_sel = dom["cath_id"].replace(".", "_").replace(",", "_").replace("/", "_").replace(" ", "_").replace(":", "_")
            plddt_str = f"{dom['plddt']:.2f}" if isinstance(dom['plddt'], (int, float)) else str(dom['plddt'])
            
            sel_name = f"{object_name}_D{idx}_{cath_clean_sel}"
            sel_query = _parse_boundary_to_selection(bounds, object_name)
            
            cmd.select(sel_name, sel_query)
            domain_selections.append(sel_name)
            
            color_name = f"ted_col_{object_name}_{idx}"
            cmd.set_color(color_name, palette[idx - 1])
            
            if str(color).lower() in ["domains", "domain", "ted"]:
                cmd.color(color_name, sel_name)

            cath_display = (dom['cath_name'][:45] + '...') if len(dom['cath_name']) > 48 else dom['cath_name']
            print(f"{sel_name:<30} | {bounds:<20} | {dom['cath_level']:<16} | {cath_display:<48} | {plddt_str}")

        print("=" * 135)
        print("Note: '⚠' denotes an automated external annotation (uncurated in CATH).")

        group_name = f"{object_name}_TED_Domains"
        cmd.group(group_name, " ".join(domain_selections))
        cmd.order(f"{object_name}_TED_Domains", "yes")

    # 6. Apply Sequential pLDDT 4-Tier Coloring
    if str(color).lower() in ["true", "plddt", "b"]:
        cmd.set_color('AF_darkblue', [0.0, 0.325490, 0.839216])
        cmd.set_color('AF_lightblue', [0.396078, 0.796078, 0.952941])
        cmd.set_color('AF_yellow', [1.0, 0.858824, 0.074510])
        cmd.set_color('AF_orange', [1.0, 0.490196, 0.270588])

        cmd.color("AF_orange", f"{object_name}")
        cmd.color("AF_yellow", f"{object_name} and b > 50")
        cmd.color("AF_lightblue", f"{object_name} and b > 70")
        cmd.color("AF_darkblue", f"{object_name} and b > 90")
        print(f"[fetchAF] Applied official AlphaFold pLDDT 4-color scheme.")

    cmd.deselect()
    cmd.orient(object_name)

cmd.extend("fetchAF", fetchAF)