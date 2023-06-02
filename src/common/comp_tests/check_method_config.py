import yaml
import requests


## VIASH START

meta = {
    "config" : "foo"
}

## VIASH END

NAME_MAXLEN = 50

SUMMARY_MAXLEN = 400

DESCRIPTION_MAXLEN = 1000

_MISSING_DOIS = ["vandermaaten2008visualizing", "hosmer2013applied"]


def assert_dict(dict, functionality):

    arg_names = []
    args = functionality["arguments"]

    for i in args:
        arg_names.append(i["name"].replace("--",""))
    
    info = functionality["info"]
    if dict:
        for key in dict:
            assert key in arg_names or info, f"{key} is not a defined argument or .functionality.info field"

def _load_bib():
    bib_path = meta["resources_dir"]+"/library.bib"
    with open(bib_path, "r") as file:
        return file.read()

def check_url(url):
    get = requests.get(url)

    assert get.status_code is (200 or 429), f"{url} is not reachable, {get.status_code}." # 429 rejected, too many requests

def search_ref_bib(reference):
    import re
    bib = _load_bib()
    
    entry_pattern =  r"(@\w+{[^}]*" + reference + r"[^}]*}(.|\n)*?)(?=@)"

    bib_entry = re.search(entry_pattern, bib)

    if bib_entry:

        type_pattern = r"@(.*){" + reference
        doi_pattern = r"(?=doi\s*=\s*{([^,}]+)})"

        entry_type = re.search(type_pattern, bib_entry.group(1))

        if not (entry_type.group(1) == "misc" or reference in _MISSING_DOIS):
            entry_doi = re.search(doi_pattern, bib_entry.group(1))
            assert entry_doi.group(1), "doi not found in bibtex reference"
            check_url(f"https://doi.org/{entry_doi.group(1)}")

        return True

    else:
        return False

print("Load config data", flush=True)
with open(meta["config"], "r") as file:
    config = yaml.safe_load(file)


print("Check general fields", flush=True)
assert len(config["functionality"]["name"]) <= NAME_MAXLEN, f"Component id (.functionality.name) should not exceed {NAME_MAXLEN} characters."
assert "namespace" in config["functionality"] is not None, "namespace not a field or is empty"

print("Check info fields", flush=True)
info = config['functionality']['info']
assert "type" in info, "type not an info field"
info_types = ["method", "control_method"]
assert info["type"] in info_types , f"got {info['type']} expected one of {info_types}"
assert "pretty_name" in info is not None, "pretty_name not an info field or is empty"
assert "summary" in info is not None, "summary not an info field or is empty"
assert "FILL IN:" not in info["summary"], "Summary not filled in"
assert len(info["summary"]) <= SUMMARY_MAXLEN, f"Component id (.functionality.info.summary) should not exceed {SUMMARY_MAXLEN} characters."
assert "description" in info is not None, "description not an info field or is empty"
assert "FILL IN:" not in info["description"], "description not filled in"
assert len(info["description"]) <= DESCRIPTION_MAXLEN, f"Component id (.functionality.info.description) should not exceed {DESCRIPTION_MAXLEN} characters."
if ("control" not in info["type"]):
    assert "reference" in info, "reference not an info field"
    bib = _load_bib()
    if info["reference"]:
        assert search_ref_bib(info["reference"]), f"reference {info['reference']} not added to library.bib"
    assert "documentation_url" in info is not None, "documentation_url not an info field or is empty"
    assert "repository_url" in info is not None, "repository_url not an info field or is empty"
    check_url(info["documentation_url"])
    check_url(info["repository_url"])


if "variants" in info:
    for key in info["variants"]:
        assert_dict(info["variants"][key], config['functionality'])

assert "preferred_normalization" in info, "preferred_normalization not an info field"
norm_methods = ["log_cpm", "counts", "log_scran_pooling", "sqrt_cpm", "l1_sqrt"]
assert info["preferred_normalization"] in norm_methods, "info['preferred_normalization'] not one of '" + "', '".join(norm_methods) + "'."



print("All checks succeeded!", flush=True)
