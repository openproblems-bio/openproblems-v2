import yaml


## VIASH START

meta = {
    "config" : "foo"
}

## VIASH END

NAME_MAXLEN = 50

SUMMARY_MAXLEN = 400

DESCRIPTION_MAXLEN = 1000


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
assert len(info["summary"]) <= SUMMARY_MAXLEN, "Summary is too long"
assert "description" in info is not None, "description not an info field or is empty"
assert "FILL IN:" not in info["description"], "description not filled in"
assert len(info["description"]) <= DESCRIPTION_MAXLEN, "description is too long"
if ("control" not in info["type"]):
    assert "reference" in info, "reference not an info field"
    assert "documentation_url" in info is not None, "documentation_url not an info field or is empty"
    assert "repository_url" in info is not None, "repository_url not an info field or is empty"
# assert "variants" in info,  "variants not an info field"
# TODO: if variants is in info, check whether it's a dictionary with correct arguments
assert "preferred_normalization" in info, "preferred_normalization not an info field"
norm_methods = ["log_cpm", "counts", "log_scran_pooling", "sqrt_cpm", "l1_sqrt"]
assert info["preferred_normalization"] in norm_methods, "info['preferred_normalization'] not one of '" + "', '".join(norm_methods) + "'."



print("All checks succeeded!", flush=True)
