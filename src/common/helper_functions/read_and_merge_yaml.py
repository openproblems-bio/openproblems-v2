def read_and_merge_yaml(path):
    """Read a Viash YAML
    
    If the YAML contains a "__merge__" key anywhere in the yaml,
    the path specified in that YAML will be read and the two
    lists will be merged. This is a recursive procedure.
    
    Arguments:
    path -- Path to the Viash YAML"""
    import ruamel.yaml as yaml
    print("imported ruayaml", flush=True)
    # check if path is pathlib.Path
    if hasattr(path, 'read_text'):
        data = yaml.safe_load(path.read_text())
    else:
        with open(path, 'r') as stream:
            data = yaml.safe_load(stream)
    print("loaded yaml", flush=True)
    return _ram_process_merge(data, path)

def _ram_deep_merge(dict1, dict2):
    if isinstance(dict1, dict) and isinstance(dict2, dict):
        keys = set(list(dict1.keys()) + list(dict2.keys()))
        out = {}
        for key in keys:
            if key in dict1:
                if key in dict2:
                    out[key] = _ram_deep_merge(dict1[key], dict2[key])
                else:
                    out[key] = dict1[key]
            else:
                out[key] = dict2[key]
        return out
    elif isinstance(dict1, list) and isinstance(dict2, list):
        return dict1 + dict2
    else:
        return dict2

def _ram_process_merge(data, path):
    print("processing merge", flush=True)
    import os
    print("imported os", flush=True)
    if isinstance(data, dict):
        print("is dict", flush=True)
        processed_data = {k: _ram_process_merge(v, path) for k, v in data.items()}
        print("processed data", flush=True)

        if "__merge__" in processed_data:
            print("merge key found", flush=True)
            new_data_path = os.path.join(os.path.dirname(path), processed_data["__merge__"])
            print("new data path", flush=True)
            new_data = read_and_merge_yaml(new_data_path)
            print("new data", flush=True)
        else:
            new_data = {}
        
        print("returning merge", flush=True)

        return _ram_deep_merge(new_data, processed_data)
    elif isinstance(data, list):
        return [_ram_process_merge(dat, path) for dat in data]
    else:
        return data

