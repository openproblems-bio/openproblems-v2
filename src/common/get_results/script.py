import yaml
import json
from pandas import read_csv
from datetime import timedelta
import re
import numpy as np

## VIASH START

par = {
    'input_scores': 'output/v2/batch_integration/output.run_benchmark.output_scores.yaml',
    'input_execution': 'output/v2/batch_integration/trace.txt',
    'task_id': 'batch_integration',
    'output': 'output/temp/results.json'
}

meta = {
}

## VIASH END

print('Loading scores', flush=True)
with open(par['input_scores'], 'r') as f:
    scores = yaml.safe_load(f)

print('Loading execution trace', flush=True)
execution = read_csv(par['input_execution'], sep='\t')

def organise_score (scores):
    '''
    combine all the metric values into one dictionary per method, dataset and normalization
    '''
    score_temp = {}
    for score in scores:
        score_id = score["dataset_id"] + "_" + score["method_id"] + "_" + score["normalization_id"]
        
        if score.get("metric_values") is None:
            score["metric_values"] = [None] * len(score["metric_ids"])
        for i, value in enumerate(score["metric_values"]):
            if np.isnan(value):
                score["metric_values"][i] = None
        comb_metric = zip(score["metric_ids"], score["metric_values"])
        score["metric_values"] = dict(comb_metric)
        score["task_id"] = par["task_id"]
        del score["metric_ids"]
        if score_temp.get(score_id) is None:
            score_temp[score_id] = score
        else:
            score_temp[score_id]["metric_values"].update(score["metric_values"])
    
    return score_temp


def convert_size (df, col):
    '''
        Convert the size to MB and to float type
    '''
    mask_kb = df[col].str.contains("KB")
    mask_mb = df[col].str.contains("MB")
    mask_gb = df[col].str.contains("GB")
    if mask_kb.any():
        df.loc[mask_kb, col] = df.loc[mask_kb, col].str.replace(" KB", "").astype(float)/1024

    if mask_mb.any():
        df.loc[mask_mb, col] = df.loc[mask_mb, col].str.replace(" MB", "").astype(float)

    
    if mask_gb.any():
        df.loc[mask_gb, col] = df.loc[mask_gb, col].str.replace(" GB", "").astype(float)*1024
    return df

def convert_duration(duration_str):
    '''
        Convert the duration to seconds
    '''
    components = duration_str.split(" ")
    hours = 0
    minutes = 0
    seconds = 0
    for component in components:
            milliseconds = 0
    for component in components:
        if "h" in component:
            hours = int(component[:-1])
        elif "ms" in component:
            milliseconds = float(component[:-2])
        elif "m" in component:
            minutes = int(component[:-1])
        elif "s" in component:
            seconds = float(component[:-1])
    duration = timedelta(hours=hours, minutes=minutes, seconds=seconds).total_seconds()
    return duration

def join_trace (trace, result):
    '''
        Join the Seqera trace with the scores
    '''
    id = trace["name"]
    dataset_id = None
    method_id = None
    match = re.search(r'\((.*?)\)', id)
    id_split = id.split(":")
    if len(id_split)>4:
        method_id = id_split[4]
    if match:
        group = match.group(1)
        split_group = group.split(".")
        if len(split_group)>1:
            dataset_id = split_group[0]
    for score in result:
        if method_id and method_id == result[score]["method_id"] and dataset_id in result[score]["dataset_id"]:
            result[score]["resources"] = {
                "duration_sec": trace["realtime"],
                "cpu_pct": trace["%cpu"],
                "peak_memory_mb": trace["peak_vmem"],
                "disk_read_mb": trace["rchar"],
                "disk_write_mb": trace["wchar"]
            }
    return result

print('Organising scores', flush=True)
org_scores = organise_score(scores)

print('Cleaning execution trace', flush=True)
execution = convert_size(execution, "rchar")
execution = convert_size(execution, "wchar")
execution = convert_size(execution, "peak_vmem")
execution["%cpu"].replace("%", "", regex=True, inplace=True)
execution["realtime"] = execution["realtime"].apply(convert_duration)

print('Joining traces and scores', flush=True)
traces = execution.to_dict(orient="records")
for trace in traces:
    org_scores = join_trace(trace, org_scores)

print('Writing results', flush=True)
result = []
for id in org_scores:
        result.append(org_scores[id])

with open (par['output'], 'w') as f:
    json.dump(result, f, indent=4)