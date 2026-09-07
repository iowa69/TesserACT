#!/usr/bin/env python3
"""Page the NCBI Datasets v2 API for every closed genome of a taxon.

Writes a TSV of one row per assembly. Complete genomes only, atypical excluded.
Kept separate from the download step so the panel can be chosen -- and the
choice argued about -- before 500 genomes are pulled over the wire.
"""
import json, sys, time, urllib.parse, urllib.request

API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon"

def get(url, tries=5):
    for a in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except Exception as e:
            if a == tries - 1:
                raise
            time.sleep(2 * (a + 1))

def main(taxon, source, out):
    fields = ",".join([
        "accession", "organism", "assembly_info", "assembly_stats", "checkm_info",
    ])
    tok, n, rows = None, 0, []
    while True:
        q = {
            "filters.assembly_level": "complete_genome",
            "filters.assembly_source": source,
            "filters.exclude_atypical": "true",
            "page_size": "1000",
            "returned_content": "COMPLETE",
        }
        if tok:
            q["page_token"] = tok
        d = get(f"{API}/{taxon}/dataset_report?" + urllib.parse.urlencode(q))
        for r in d.get("reports", []):
            ai = r.get("assembly_info", {}) or {}
            st = r.get("assembly_stats", {}) or {}
            bs = ai.get("biosample", {}) or {}
            rows.append({
                "accession": r.get("accession", ""),
                "organism": (r.get("organism", {}) or {}).get("organism_name", ""),
                "taxid": (r.get("organism", {}) or {}).get("tax_id", ""),
                "strain": (r.get("organism", {}) or {}).get("infraspecific_names", {}).get("strain", ""),
                "biosample": bs.get("accession", ""),
                "level": ai.get("assembly_level", ""),
                "refseq_category": ai.get("refseq_category", ""),
                "release_date": ai.get("release_date", ""),
                "submitter": ai.get("submitter", ""),
                "total_len": st.get("total_sequence_length", ""),
                "n_contigs": st.get("number_of_contigs", ""),
                "gc": st.get("gc_percent", ""),
                # Clinical provenance -- the isolation source is how a clinical
                # isolate is told from an environmental or veterinary one, and
                # the whole clonality argument rests on the clinical subset.
                "host": next((a.get("value","") for a in (bs.get("attributes") or [])
                              if a.get("name") in ("host","specific_host")), ""),
                "isolation_source": next((a.get("value","") for a in (bs.get("attributes") or [])
                              if a.get("name") in ("isolation_source","isolation-source")), ""),
                "geo": next((a.get("value","") for a in (bs.get("attributes") or [])
                              if a.get("name") in ("geo_loc_name","geographic location")), ""),
                "collection_date": next((a.get("value","") for a in (bs.get("attributes") or [])
                              if a.get("name") in ("collection_date",)), ""),
            })
        n += len(d.get("reports", []))
        tok = d.get("next_page_token")
        sys.stderr.write(f"  {n} assemblies\n")
        if not tok:
            break
    cols = list(rows[0].keys())
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]).replace("\t", " ").replace("\n", " ") for c in cols) + "\n")
    sys.stderr.write(f"wrote {len(rows)} rows -> {out}\n")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
