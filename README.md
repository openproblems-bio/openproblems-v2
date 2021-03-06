opsca-viash
================

-   [Requirements](#requirements)
-   [Quick start](#quick-start)
-   [Project structure](#project-structure)
-   [Adding a viash component](#adding-a-viash-component)
-   [Running a component from CLI](#running-a-component-from-cli)
-   [Building a component](#building-a-component)
-   [Unit testing a component](#unit-testing-a-component)
-   [Frequently asked questions](#frequently-asked-questions)
-   [Benefits of using Nextflow +
    viash](#benefits-of-using-nextflow--viash)

Proof Of Concept in adapting [Open Problems
repository](https://github.com/openproblems-bio/openproblems) with
Nextflow and viash. Documentation for viash is available at
[viash.io](https://viash.io).

## Requirements

To use this repository, please install the following dependencies:

-   Bash
-   Java (Java 8 or higher)
-   Docker (Instructions [here](https://docs.docker.com/get-docker/))
-   Nextflow (Optional, though [very easy to
    install](https://www.nextflow.io/index.html#GetStarted))

## Quick start

The `src/` folder contains modular software components for running a
modality alignment benchmark. Running the full pipeline is quite easy.

**Step 0, fetch viash and nextflow:** run the `bin/init` executable.

``` bash
bin/init
```

    > Using tag develop
    > Cleanup
    > Install viash develop under /viash_automount/home/rcannood/workspace/opsca/opsca-viash/bin/
    > Fetching components sources (version develop)
    > Building components
    > Done, happy viash-ing!
    > Nextflow installation completed.

**Step 1, build all the components:** in the `src/` folder as standalone
executables in the `target/` folder. Use the `-q 'xxx'` parameter to
build only several components of the repository.

``` bash
bin/viash_build -q 'modality_alignment|common'
```

    Exporting src/modality_alignment/metrics/knn_auc/ (modality_alignment/metrics) =nextflow=> target/nextflow/modality_alignment/metrics/knn_auc
    Exporting src/modality_alignment/methods/scot/ (modality_alignment/methods) =nextflow=> target/nextflow/modality_alignment/methods/scot
    Exporting src/modality_alignment/methods/mnn/ (modality_alignment/methods) =docker=> target/docker/modality_alignment/methods/mnn
    ...

These standalone executables you can give to somebody else, and they
will be able to run it, provided that they have Bash and Docker
installed. The command might take a while to run, since it is building a
docker container for each of the components. If you???re interested in
building only a subset of components, you can apply a regex to the
selected components. For example:
`bin/viash_build -q 'common|modality_alignment'`.

**Step 2, run the pipeline with nextflow.** To do so, run the bash
script located at `src/modality_alignment/workflows/run_nextflow.sh`:

``` bash
src/modality_alignment/workflows/run_nextflow.sh 
```

    [15/84d27c] process > get_scprep_csv_datasets:scprep_csv:scprep_csv_process (CBMC_8K_13AB_10x) [100%] 1 of 1 ???
    [2f/318ad9] process > mnn:mnn_process (CBMC_8K_13AB_10x)                                       [100%] 1 of 1 ???
    [6f/dd22c1] process > scot:scot_process (CBMC_8K_13AB_10x)                                     [100%] 1 of 1 ???
    [c6/6e0999] process > knn_auc:knn_auc_process (CBMC_8K_13AB_10x.scot)                          [100%] 2 of 2 ???
    [73/f6dffa] process > mse:mse_process (CBMC_8K_13AB_10x.scot)                                  [100%] 2 of 2 ???
    [54/c89b95] process > extract_scores:extract_scores_process (combined)                         [100%] 1 of 1 ???
    Completed at: 20-Apr-2021 06:52:19
    Duration    : 3m 40s
    CPU hours   : 0.1
    Succeeded   : 8

## Project structure

    bin/                     Helper scripts for building the project and developing a new component.
    src/                     Source files for each component in the pipeline.
      modality_alignment/    Source files related to the 'Modality alignment' task.
        datasets/            Dataset downloader components.
        methods/             Modality alignment method components.
        metrics/             Modality alignment metric components.
        utils/               Utils functions.
        workflow/            The pipeline workflow for this task.
      common/                 Helper files.
    target/                  Executables generated by viash based on the components listed under `src/`.
      docker/                Bash executables which can be used from a terminal.
      nextflow/              Nextflow modules which can be used in a Nextflow pipeline.
    work/                    A working directory used by Nextflow.
    output/                  Output generated by the pipeline.

## Adding a viash component

[`viash`](https://github.com/data-intuitive/viash) allows you to create
pipelines in Bash or Nextflow by wrapping Python, R, or Bash scripts
into reusable components. You can start creating a new component using
`bin/skeleton`. For example to create a new Python-based viash component
in the `src/modality_alignment/methods/foo` folder, run: You can start
creating a new component by using the `bin/skeleton` command:

``` bash
bin/viash_skeleton --name foo --namespace "modality_alignment/methods" --language python
```

This should create a few files in this folder:

    script.py                A python script for you to edit.
    config.vsh.yaml          Metadata for the script containing info on the input/output arguments of the component.
    test.py                  A python script with which you can start unit testing your component.

The [Getting started](http://www.data-intuitive.com/viash_docs/) page on
the viash documentation site provides some information on how a basic
viash component works, or on the specifications of the `config.vsh.yaml`
[config file](http://www.data-intuitive.com/viash_docs/config/).

## Running a component from CLI

You can view the interface of the executable by running the executable
with the `-h` parameter.

``` bash
viash run src/modality_alignment/methods/foo/config.vsh.yaml -- -h
```

    foo 0.0.1
    Replace this with a (multiline) description of your component.

    Options:
        -i, --input
            type: file, required parameter
            Describe the input file.

        -o, --output
            type: file, required parameter, output
            Describe the output file.

        --option
            type: string
            default: default-
            Describe an optional parameter.

You can **run the component** as follows:

``` bash
viash run src/modality_alignment/methods/foo/config.vsh.yaml -- -i LICENSE -o foo_output.txt
```

    This is a skeleton component
    The arguments are:
     - input:  /viash_automount/home/rcannood/workspace/opsca/opsca-viash/LICENSE
     - output:  /viash_automount/home/rcannood/workspace/opsca/opsca-viash/foo_output.txt
     - option:  default-

## Building a component

`viash` has several helper functions to help you quickly develop a
component.

With **`viash build`**, you can turn the component into a standalone
executable. This standalone executable you can give to somebody else,
and they will be able to run it, provided that they have Bash and Docker
installed.

``` bash
viash build src/modality_alignment/methods/foo/config.vsh.yaml \
  -o target/docker/modality_alignment/methods/foo
```

Note that the `bin/viash_build` component does a much better job of
setting up a collection of components. You can filter which components
will be built by providing a regex to the `-q` parameter,
e.g.??`bin/viash_build -q 'common|modality_alignment'`.

You can now view the same interface of the executable by running the
executable with the `-h` parameter.

``` bash
target/docker/modality_alignment/methods/foo/foo -h
```

    foo 0.0.1
    Replace this with a (multiline) description of your component.

    Options:
        -i, --input
            type: file, required parameter
            Describe the input file.

        -o, --output
            type: file, required parameter, output
            Describe the output file.

        --option
            type: string
            default: default-
            Describe an optional parameter.

Or **run the component** as follows:

``` bash
target/docker/modality_alignment/methods/foo/foo -i LICENSE -o foo_output.txt
```

    This is a skeleton component
    The arguments are:
     - input:  /viash_automount/home/rcannood/workspace/opsca/opsca-viash/LICENSE
     - output:  /viash_automount/home/rcannood/workspace/opsca/opsca-viash/foo_output.txt
     - option:  default-

## Unit testing a component

Provided that you wrote a script that allows you to test the
functionality of a component, you can run the tests by using the
**`viash test`** command.

``` bash
viash test src/modality_alignment/methods/foo/config.vsh.yaml
```

    Running tests in temporary directory: '/home/rcannood/workspace/viash_temp/viash_test_foo18431554913355206711'
    ====================================================================
    +/home/rcannood/workspace/viash_temp/viash_test_foo18431554913355206711/build_executable/foo --verbosity 6 ---setup cachedbuild
    [notice] Running 'docker build -t modality_alignment/methods_foo:da7Qhr5nLi6A /home/rcannood/workspace/viash_temp/viashsetupdocker-foo-4h6urx'
    Sending build context to Docker daemon  22.53kB

    Step 1/2 : FROM python:3.9.3-buster
     ---> 05034335a2e3
    Step 2/2 : RUN pip install --upgrade pip &&   pip install --no-cache-dir "numpy"
     ---> Using cache
     ---> 45db33ebb9de
    Successfully built 45db33ebb9de
    Successfully tagged modality_alignment/methods_foo:da7Qhr5nLi6A
    ====================================================================
    +/home/rcannood/workspace/viash_temp/viash_test_foo18431554913355206711/test_test.py/test.py
    >> Writing test file
    >> Running component
    >> Checking whether output file exists
    >> Checking contents of output file
    >> All tests succeeded successfully!
    ====================================================================
    [32mSUCCESS! All 1 out of 1 test scripts succeeded![0m
    Cleaning up temporary directory

To run all the unit tests of all the components in the repository, use
`bin/viash_test`.

## Frequently asked questions

### My component doesn???t work!

Debugging your component based on the output from a Nextflow pipeline is
easier than you might realise. For example, the error message below
tells you that the ???mse??? component failed:

``` bash
src/modality_alignment/workflows/run_nextflow.sh 
```

    N E X T F L O W  ~  version 20.10.0
    [f8/2acb9f] process > get_scprep_csv_datasets:scprep_csv:scprep_csv_process (CBMC_8K_13AB_10x) [100%] 1 of 1 ???
    [6e/0cb81b] process > mnn:mnn_process (CBMC_8K_13AB_10x)                                       [100%] 1 of 1 ???
    [43/edc9a1] process > scot:scot_process (CBMC_8K_13AB_10x)                                     [100%] 1 of 1 ???
    [00/41ee55] process > knn_auc:knn_auc_process (CBMC_8K_13AB_10x.scot)                          [100%] 2 of 2 ???
    [3d/0d6afe] process > mse:mse_process (CBMC_8K_13AB_10x.scot)                                  [100%] 2 of 2, failed: 2 ???
    [22/5899a9] process > extract_scores:extract_scores_process (combined)                         [100%] 1 of 1 ???
    [3d/0d6afe] NOTE: Process `mse:mse_process (CBMC_8K_13AB_10x.scot)` terminated with an error exit status (1) -- Error is ignored
    Completed at: 19-Apr-2021 20:09:22
    Duration    : 3m 46s
    CPU hours   : 0.1 (2.7% failed)
    Succeeded   : 6
    Ignored     : 2
    Failed      : 2

Looking at this output reveals in which step of the pipeline the ???mse???
component failed, namely `3d/0d6afe`. This means we should check a
folder called `work/3d/0d6afe...`:

``` bash
ls -la work/3d/0d6afe9c27ab68d3f10551c3d3104c/
```

    total 28
    drwxrwxr-x. 1 rcannood rcannood  216 Apr 19 20:09 .
    drwxrwxr-x. 1 rcannood rcannood   60 Apr 19 20:09 ..
    lrwxrwxrwx. 1 rcannood rcannood  108 Apr 19 20:09 CBMC_8K_13AB_10x.scot.h5ad
    -rw-rw-r--. 1 rcannood rcannood    0 Apr 19 20:09 .command.begin
    -rw-rw-r--. 1 rcannood rcannood  191 Apr 19 20:09 .command.err
    -rw-rw-r--. 1 rcannood rcannood  262 Apr 19 20:09 .command.log
    -rw-rw-r--. 1 rcannood rcannood   71 Apr 19 20:09 .command.out
    -rw-rw-r--. 1 rcannood rcannood 3224 Apr 19 20:09 .command.run
    -rw-rw-r--. 1 rcannood rcannood  463 Apr 19 20:09 .command.sh
    -rw-rw-r--. 1 rcannood rcannood    1 Apr 19 20:09 .exitcode

``` bash
cat work/3d/0d6afe9c27ab68d3f10551c3d3104c/.command.err 
```

    Traceback (most recent call last):
      File "/tmp/viash-run-mse-WausLu", line 39, in <module>
        adata.uns["metric_value"] = area_under_curve
    NameError: name 'area_under_curve' is not defined

It seems that some error occurred within the Python script. Luckiky, the
input file of this process is in this directory. We can manually run the
component by running:

``` bash
viash run src/modality_alignment/metrics/mse/config.vsh.yaml -- \
  -i work/3d/0d6afe9c27ab68d3f10551c3d3104c/CBMC_8K_13AB_10x.scot.h5ad -o test.h5ad
```

Alternatively, you can edit
`src/modality_alignment/metrics/mse/script.py` and replace the header
by:

``` python
## VIASH START
# The code between the the comments above and below gets stripped away before 
# execution. Here you can put anything that helps the prototyping of your script.
par = {
    "input": "work/3d/0d6afe9c27ab68d3f10551c3d3104c/CBMC_8K_13AB_10x.scot.h5ad",
    "output": "test.h5ad"
}
## VIASH END

## ... the rest of the script ...
```

Now you can work on the `script.py` file in your preferred editor
(vim?). For easy prototyping, viash will automatically strip away
anything between the `## VIASH START` and `## VIASH END` codeblock at
runtime.

## Benefits of using Nextflow + viash

### The pipeline is **language-agnostic**

This means that each component can be written in whatever scripting
language the user desires. Here are examples of a
[Python](src/modality_alignment/methods/scot/) and an
[R](src/modality_alignment/methods/mnn) component.

By default, viash supports wrapping the following scripting languages:
Bash, Python, R, JavaScript, and Scala. If viash doesn???t support your
preferred scripting language, feel free to ask the developers to [add
it](https://github.com/data-intuitive/viash/issues). Alternatively, you
can write a Bash script which calls your desired programming language.

### One Docker container per component

By running the `bin/viash_build` command, viash will build one Docker
container per component. While this results in some initial
computational overhead, this makes it a lot easier to add a new
component to the pipeline with dependencies which might conflict with
those of other components.

### Reproducible components

A component built by viash is meant to be reproducible. If you send the
`target/docker/modality_alignment/methods/foo/foo` file to someone, they
can run `./foo ---setup cachedbuild` and then will be able to use the
`foo` component however they like.

``` bash
# pretend to send the component to someone through 'cp'
cp target/docker/modality_alignment/methods/foo/foo foo_by_email

# build container
./foo_by_email ---setup cachedbuild
```

    [notice] Running 'docker build -t modality_alignment/methods_foo:0.0.1 /home/rcannood/workspace/viash_temp/viashsetupdocker-foo-TsNVGL'

``` bash
# view help
./foo_by_email -h
```

    foo 0.0.1
    Replace this with a (multiline) description of your component.

    Options:
        -i, --input
            type: file, required parameter
            Describe the input file.

        -o, --output
            type: file, required parameter, output
            Describe the output file.

        --option
            type: string
            default: default-
            Describe an optional parameter.

``` bash
# run component
./foo_by_email -i LICENSE -o foo_output.txt
```

    This is a skeleton component
    The arguments are:
     - input:  /viash_automount/home/rcannood/workspace/opsca/opsca-viash/LICENSE
     - output:  /viash_automount/home/rcannood/workspace/opsca/opsca-viash/foo_output.txt
     - option:  default-

### Reprodicible components on Docker Hub

You might notice that the `---setup cachedbuild` builds the docker
container from scratch, rather than pulling it from Docker hub.

With `bin/viash_build`, you can build a versioned release of all the
components in the repository and push it to Docker hub.

``` bash
bin/viash_build -m release -v '0.1.0' -r singlecellopenproblems
```

    In release mode...
    Exporting src/modality_alignment/methods/mnn/ (modality_alignment/methods) =docker=> target/docker/modality_alignment/methods/mnn
    Exporting src/modality_alignment/methods/scot/ (modality_alignment/methods) =docker=> target/docker/modality_alignment/methods/scot
    Exporting src/common/extract_scores/ (common) =docker=> target/docker/common/extract_scores
    Exporting src/modality_alignment/metrics/knn_auc/ (modality_alignment/metrics) =docker=> target/docker/modality_alignment/metrics/knn_auc
    Exporting src/modality_alignment/datasets/scprep_csv/ (modality_alignment/datasets) =docker=> target/docker/modality_alignment/datasets/scprep_csv
    Exporting src/modality_alignment/metrics/mse/ (modality_alignment/metrics) =docker=> target/docker/modality_alignment/metrics/mse
    > docker build -t singlecellopenproblems/common_extract_scores:0.1.0 --no-cache /home/rcannood/workspace/viash_temp/viashsetupdocker-extract_scores-FyHtgS
    > docker build -t singlecellopenproblems/modality_alignment/metrics_mse:0.1.0 --no-cache /home/rcannood/workspace/viash_temp/viashsetupdocker-mse-r2LSpO
    > docker build -t singlecellopenproblems/modality_alignment/metrics_knn_auc:0.1.0 --no-cache /home/rcannood/workspace/viash_temp/viashsetupdocker-knn_auc-S8dJP5
    > docker build -t singlecellopenproblems/modality_alignment/datasets_scprep_csv:0.1.0 --no-cache /home/rcannood/workspace/viash_temp/viashsetupdocker-scprep_csv-lItAG1
    > docker build -t singlecellopenproblems/modality_alignment/methods_scot:0.1.0 --no-cache /home/rcannood/workspace/viash_temp/viashsetupdocker-scot-xUKof3
    > docker build -t singlecellopenproblems/modality_alignment/methods_mnn:0.1.0 --no-cache /home/rcannood/workspace/viash_temp/viashsetupdocker-mnn-0rjhKc

The images themselves can be pushed to Docker Hub with the
`bin/viash_push` command. I???d have to make a small change to viash to
ensure that the component names don???t contain any slashes because the
images listed above can???t be pushed to Docker hub. However, the output
would look something like this:

``` bash
bin/viash_push -m release -v '0.1.0' -r openproblems
In release mode...
Using version 0.1.0 to tag containers
```

    > openproblems/modality_alignment_metrics_knn_auc:0.1.0 does not exist, try pushing ... OK!
    > openproblems/modality_alignment_methods_scot:0.1.0 does not exist, try pushing ... OK!
    > openproblems/modality_alignment_metrics_mse:0.1.0 does not exist, try pushing ... OK!
    > openproblems/modality_alignment_datasets_scprep_csv:0.1.0 does not exist, try pushing ... OK!
    > openproblems/common_extract_scores:0.1.0 does not exist, try pushing ... OK!
    > openproblems/modality_alignment_methods_mnn:0.1.0 does not exist, try pushing ... OK!

<!-- cleaning up temporary files -->
