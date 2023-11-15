import sys
import os
import pytest
import anndata as ad


def test_cellxgene_source_download(run_component):

    run_component([
        "--dataset_id", "ea0a934f-799c-4b2d-b66b-9d04ca72356a",
        "--cellxgene_release", "2023-07-25",
        "--output", "output.h5ad",
    ])

    # check whether file exists
    assert os.path.exists("output.h5ad"), "Output file does not exist"

if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
