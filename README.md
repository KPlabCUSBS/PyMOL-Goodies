fetchAF.py allows you to directly fetch protein structures from AlphaFold using the uniprot ID in addition to exposing multiple useful arguments in pymol.
Usage
fetchAF uniprot_id [, color [, show_domains [, save_dir ]]]

Here the argument color can be 'color=domains' or 'color=plddt'.

When 'plddt' is used as the argument, the fetched structure is colored according to pLDDT values using the standard AlphaFold coloring code (Dark Blue (> 90, Very high confidence); Light Blue (70 – 90, High confidence); Yellow (50 – 70; Low confidence); Orange (< 50, Very low confidence).

When 'domains' is used as the argument, the information from fetched structure is used to talk to CATH and TED databases to identify domain boundaries and color each domain with a different color. It also outputs the different domains with boundaries as a table.


af2_curate.py is similar to fetchAF.py in action but it works on any protein already loaded in the pymol object list.
