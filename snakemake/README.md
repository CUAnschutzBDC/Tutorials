# How to build and run a snakemake pipeline

Snakemake is a "tool to create reproducible and scalable data analysis". Snakemake is written in python and can be easily used on different servers and in different envrionments. I like it because it allows me to be confident that I am doing the exact same analysis on all samples and I can run through all steps of the analysis with one single submission script. It is also very easy to update between projects and you often need to edit just one configuration file.

## Quick tips
* I always do a dry run first. This way if there are any errors you can quickly see those and fix them.
* Use `expand` to analyze multiple samples without needing to type out all the output files individually.
* Use `rule all` as a driver of all output files you want generated.
* Use `params` to pass arguments you don't want to hard code

## Installation

The recommendation is to install snakemake with `conda`.

```bash
$ conda install -c bioconda -c conda-forge snakemake
``` 

Personally, I like to put conda into it's own environment. I also use `mamba` because the packages are more up to date than `conda` so I run

```bash
$ conda create -n snakemake
$ conda activate snakemake
$ conda install -c conda-forge mamba
$ mamba install -c bioconda -c conda-forge snakemake
```

If you do it this way, you will need to activate snakemake before running snakemake or submitting a script that runs snakemake.

```bash
$ conda activate snakemake
```

## Basic usage
To use snakemake you write "rules". These rules are individual steps of your analysis pipeline. For example for my RNA-seq analysis, I have a rule to run fastqc, a rule to make a summary of fastqc output, a rule to run star, a rule to make a summary of star output, a rule to run featureCounts, and a rule to make a count table. But let's start with a simple rule.

All snakemake pipelines need a file called the "Snakefile" *This is not 100% necessary, but if you don't call it Snakefile you must put the name of your snakefile into the command when calling snakemake*. This Snakefile contains all of the code and rules necessary to run your pipeline. There are two options for how to write your rules.

1. You can write all of your rules directly into your Snakefile. This is the way I did it in grad school so you can see this option on my github.
2. You can make scripts for sets of rules (like run fastqc and write a summary) and link to them from the snakefile. I personally like this option because it makes your code much more readable. If someone wants to know how you aligned your reads, they can just look at the alignment file rather than searching through all of your rules. It also makes it easier to share rules between pipelines. For example, you are working with a sample that has both mouse and human sequences so you need to change your alignment step, you can keep the fastqc, adapter trimming, and counting rules and just write a new alignment rule. If you choose this method, you need to link to the files from your Snakefile. 

```python
include: "src/rules/rule_set.snake"
```

Snakemake recommends a very specific orginazitional strategy that I like, but you don't need to follow. They recommend the scripts go into a directory named `src` and rules go into a directory named `src/rules`. Analysis scripts can also go into the `src` directory. For reasons we can discuss later I put them here: `src/rules/scripts/`. If you chose to have rule files, they have the `.snake` extension.

### How to write a rule
*This corresponds to the rule `first_rule.snake`*

Let's make a rule that will write the name of a sample to a file. To write the first rule, you need the following structure

```python
rule make_file:
    output:
        "results/first_file.txt"
    shell:
        """
        echo "this is a file" > {output}
        """
```

A few things to notice. 

1. You always need to name your rules (and the names need to be unique). 
2. You always need to have an output argument
3. You can run shell commands with the shell argument. If you run a shell command you need the three """ at the start and end of the command
4. You can directly call the output file without needing to type the path again.
5. Any shell command (or set of commands) can work. As long as all of your commands are within the """, you can run as complex of a shell script as you want here.

We can now test what will happen by running this rule using the -np flag
* -n means "dry-run", don't execute
* -p prints out the shell commands that will be excuted

```bash
$ snakemake -np results/first_file.txt
```

```
Building DAG of jobs...
Job counts:
    count   jobs
    1   make_file
    1
[Mon Mar 15 14:39:29 2021]
rule make_file:
    output: results/first_file.txt
    jobid: 0
        echo "this is a file" > results/first_file.txt
Job counts:
    count   jobs
    1   make_file
    1
This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
```

Here we told snakemake what command to run by putting the path to the `output` file, here that was `results/first_file.txt`. Notice that it tells you both the name of the job it will run and the command it will run. At the end it also tells you "job counts" which is an overview of all the jobs that will be run. 

To actually run this, we can run
```bash
$ snakemake results/first_file.txt --cores 1
```

```
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 1 (use --cores to define parallelism)
Rules claiming more threads will be scaled down.
Job counts:
    count   jobs
    1   make_file
    1
Select jobs to execute...
[Mon Mar 15 14:48:16 2021]
rule make_file:
    output: results/first_file.txt
    jobid: 0
[Mon Mar 15 14:48:16 2021]
Finished job 0.
1 of 1 steps (100%) done
Complete log: /beevol/home/wellskri/Analysis/Tutorials/snakemake/.snakemake/log/2021-03-15T144815.335079.snakemake.log
```

This gives you the output of the job. If the job runs successfully, it will tell you the job completed successfully. It also tells you the path to the log file that saved this output. The log files are in hidden directories.

To see the output
```bash
$ cd results
$ ls
```

```
first_file.txt
```

```bash
$ cd ../
```

One more important point, I didn't have to create a results directory. Snakemake will create output directories for you as long as the directory is part of the path in your `output` argument.

### Generalizing a rule
*This corresponds to the rule `general_rule.snake`*

In addition to writing a rule for one specific file, we can also make rules more general by using something called `wildcards`. One of the most common ways I use wildcards is to make a rule general to any samples, but there are other uses too.

```python
rule make_file_general:
    output:
        "results/{sample}_first_file.txt"
    shell:
        """
        echo {wildcards.sample} > {output}
        """
```

Notice that you can use a wildcard in the name of the output file and you can call the wildcard as a variable in the shell command.

We can again run a dry test of this

```bash
$ snakemake -np results/sample_1_first_file.txt
```

```
Building DAG of jobs...
Job counts:
    count   jobs
    1   make_file_general
    1
[Mon Mar 15 14:47:04 2021]
rule make_file_general:
    output: results/sample_1_first_file.txt
    jobid: 0
    wildcards: sample=sample_1
        echo sample_1 > results/sample_1_first_file.txt
Job counts:
    count   jobs
    1   make_file_general
    1
This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
```

Again, it tells us what will be done. Here you can see that the job name is different than the pevious job, that was "make_file" while this is "make_file_general".

Again, we can run it
```{bash}
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 1 (use --cores to define parallelism)
Rules claiming more threads will be scaled down.
Job counts:
    count   jobs
    1   make_file_general
    1
Select jobs to execute...
[Mon Mar 15 14:52:18 2021]
rule make_file_general:
    output: results/sample_1_first_file.txt
    jobid: 0
    wildcards: sample=sample_1
[Mon Mar 15 14:52:18 2021]
Finished job 0.
1 of 1 steps (100%) done
Complete log: /beevol/home/wellskri/Analysis/Tutorials/snakemake/.snakemake/log/2021-03-15T145218.383705.snakemake.log
```

To see the output
```{bash}
$ cd results
$ ls
```

```
first_file.txt  sample_1_first_file.txt
```

Now we can see both files present.

If the file has already been created, snakemake will not run the rule again

```bash
$ snakemake results/sample_1_first_file.txt --cores 1
```

```
Building DAG of jobs...
Nothing to be done.
Complete log: /beevol/home/wellskri/Analysis/Tutorials/snakemake/.snakemake/log/2021-03-15T150541.044925.snakemake.log
```

### Creating multiple rules
*This corresponds to the rule `general_rule.snake`*

One of the powers of snakemake is that multiple rules can be strung together. Here, the output file for one rule becomes the input file for another rule. Before running the second rule, Snakemake will first check that the input file exists. If it doesn't, it will create the input file by running the first rule.

Let's look at an example using our first rule

```python
rule make_file_general:
    output:
        "results/{sample}_first_file.txt"
    shell:
        """
        echo {wildcards.sample} > {output}
        """
rule make_second_file:
    input:
        "results/{sample}_first_file.txt"
    output:
        "results/{sample}_second_file.txt"
    shell:
        """
        echo {wildcards.sample} > {output}
        """
```

Here we now have an input file. Just like with the output, we can include wildcards with the input file. Make sure the same wildcard is used in both the input and output files (there are some exceptions that we will get to later).

If snakemake can't locate the input file, it will look to see if any other rules will generate that input file. If it can't find any rules to generate that file, it will throw and error.

```
MissingInputException in line 18 of /beevol/home/wellskri/Analysis/Tutorials/snakemake/Snakefile:
Missing input files for rule make_second_file:
results/sample_1_first_file1.txt
```

If you see an error like this, make sure that you have a rule with an output file that is identical to the input file.

To test this rule we can run
```bash
$ snakemake -np results/sample_1_second_file.txt
```

```
Building DAG of jobs...
Job counts:
    count   jobs
    1   make_second_file
    1
[Mon Mar 15 15:06:45 2021]
rule make_second_file:
    input: results/sample_1_first_file.txt
    output: results/sample_1_second_file.txt
    jobid: 0
    wildcards: sample=sample_1
        echo sample_1 > results/sample_1_second_file.txt
Job counts:
    count   jobs
    1   make_second_file
    1
This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
```

Here you can see that only the second rule will be run. But a great thing about snakemake is it will run everything by checking that the input files exist. To show how this works, I will delete my results folder

```bash
$ rm -r results
$ snakemake -np results/sample_1_second_file.txt
```

```
Building DAG of jobs...
Job counts:
    count   jobs
    1   make_file_general
    1   make_second_file
    2
[Mon Mar 15 15:03:38 2021]
rule make_file_general:
    output: results/sample_1_first_file.txt
    jobid: 1
    wildcards: sample=sample_1
        echo sample_1 > results/sample_1_first_file.txt
[Mon Mar 15 15:03:38 2021]
rule make_second_file:
    input: results/sample_1_first_file.txt
    output: results/sample_1_second_file.txt
    jobid: 0
    wildcards: sample=sample_1
        echo sample_1 > results/sample_1_second_file.txt
Job counts:
    count   jobs
    1   make_file_general
    1   make_second_file
    2
This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
```

Now you can see that even though we only referred to the second file, snakemake will run both rules. 

Again to run it
```bash
$ snakemake results/sample_1_second_file.txt --cores 1
```

```
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 1 (use --cores to define parallelism)
Rules claiming more threads will be scaled down.
Job counts:
    count   jobs
    1   make_second_file
    1
Select jobs to execute...
[Mon Mar 15 15:08:46 2021]
rule make_second_file:
    input: results/sample_1_first_file.txt
    output: results/sample_1_second_file.txt
    jobid: 0
    wildcards: sample=sample_1
[Mon Mar 15 15:08:46 2021]
Finished job 0.
1 of 1 steps (100%) done
Complete log: /beevol/home/wellskri/Analysis/Tutorials/snakemake/.snakemake/log/2021-03-15T150846.110197.snakemake.log
```

We can check the output
```bash
$ cd results
$ ls
```

```
sample_1_first_file.txt  sample_1_second_file.txt
```

### Running multiple samples
*This corresponds to the rule `general_rule.snake`*

Snakemake can also help us combine multiple samples and run them at the same time. For example, you can be running an analysis and the last step is to create a combined summary of all the runs. For this, you could create a rule where the input file is the final output from each individual sample. To create this rule you could either list all the expected output samples

```python
rule combine_samples:
    input:
        "results/sample_1_second_file.txt",
        "results/sample_2_second_file.txt",
        "results/sample_3_second_file.txt"
```

But this isn't very flexible if you add a new sample or if you want to copy to a new analysis. It also becomes tiring if you want to do this for many samples. Another option is to use `expand`. To use `expand`, first include a list of samples in your snakefile (later we will put these into a `configfile`)

First we add a list to the snakefile
```python
samples=["sample_1", "sample_2", "sample_3", "sample_4"]
```

Then we can write our rule
```python
rule combine_samples:
    input:
        expand(
            "results/{sample}_second_file.txt",
            sample = samples
        )
    output:
        "results/combined_file.txt"
    shell:
        """
        cat {input} > {output}
        """
```

*Note: here we do not have the same wildcards for our input and output files*

In this rule we use something called `expand`. This basically will return a list replacing the wildcards with the samples in your list:

```
("results/sample_1_second_file", "results/sample_2_second_file",
"results/sample_3_second_file", "results/sample_4_second_file")
```

If you use two wildcards with expand, it will return all pairwise options:
```python
samples = ["sample_1", "sample_2"]
types = ["file", "matrix"]
expand("results/{sample}_second_{type}", sample = samples, type = types)
```

gives

```
("results/sample_1_second_file", "results/sample_1_second_matrix",
"results/sample_2_second_file", "results/sample_2_second_matrix")
```

And do a test run again
```bash
$ snakemake -np results/combined_file.txt
```

```
Building DAG of jobs...
Job counts:
    count   jobs
    1   combine_samples
    3   make_file_general
    3   make_second_file
    7
[Mon Mar 15 16:28:12 2021]
rule make_file_general:
    output: results/sample_3_first_file.txt
    jobid: 6
    wildcards: sample=sample_3
        echo sample_3 > results/sample_3_first_file.txt
[Mon Mar 15 16:28:12 2021]
rule make_file_general:
    output: results/sample_4_first_file.txt
    jobid: 8
    wildcards: sample=sample_4
        echo sample_4 > results/sample_4_first_file.txt
...
```


I'm not going to paste the entire output of this command, we can see it when we run it, but a few things I want to draw your attention to. First, look at the job counts table. This shows that 3 separate rules will be run: `combine_samples`, `make_file_general`, and `make_second_file`. Second, look at the count column, this shows that `combine_samples` will be run once, `make_file_general` and `make_second_file` will be run 3 times. This is because we need to run the first two rules for samples 2-4 (we already ran these rules for sample 1). Here you can also see that it is taking the sample names from our list and using them as the sample `wildcards`. 

We can again run this

```bash
$ snakemake --cores 1 results/combined_file.txt
$ cd results
$ more combined_file.txt
```

```
sample_1
sample_2
sample_3
sample_4
```

Here we see the output of each individual sample file combined into one file.

### Rule all
*This corresponds to the rule `rule_all_ex.snake`*

Another way to run many samples is by using the `rule all`. `rule all` only takes an input argument and no other arguments.

Here is an example

```python
rule four:
    input: "results/{sample}_second_file.txt"
    output: "results/{sample}_fourth_file.txt"
    shell:
        """
        echo {wildcards.sample} > {output}
        """
rule all:
    input:
        expand(
            "results/{sample}_fourth_file.txt",
            sample = samples
        )
```

And again we can run it with a dry run first
```bash
$ snakemake -np
```

```
Building DAG of jobs...
Job counts:
    count   jobs
    1   all
    4   four
    5
[Mon Mar 15 17:56:10 2021]
rule four:
    input: results/sample_3_second_file.txt
    output: results/sample_3_fourth_file.txt
    jobid: 7
    wildcards: sample=sample_3
        echo sample_3 > results/sample_3_fourth_file.txt
[Mon Mar 15 17:56:10 2021]
rule four:
    input: results/sample_4_second_file.txt
    output: results/sample_4_fourth_file.txt
    jobid: 10
    wildcards: sample=sample_4
        echo sample_4 > results/sample_4_fourth_file.txt
...
```

There are a couple of things I want to draw attention to here.

1. When you use `rule all` you don't need to specify an output file when calling snakemake. It will look at the output files in `rule all` and generate those
2. One of the jobs listed is `all`. This isn't really going to do anything, but the rule exists to tell snakemake what output files you want to generate.
3. Again here you can see that rule four is being run for all of our samples because we said that we wanted all of our samples run through that rule.
4. We don't generate an extra output file that combines all the outputs from `rule four`, instead we just say we want these files as an output.


`rule all` can be used to walk though your entire pipeline. Let's say here you want `rule four` run on individual samples (but don't plan on combining) and you want `rule combine` run as well. You can add the output from `rule combine` to your `rule all` as well. You can add as many inputs to `rule all` as you want. To add multiple inputs, just separate them by a comma

```python
rule all:
    input:
        # Create the fourth file for all samples
        expand(
            "results/{sample}_fourth_file.txt",
            sample = samples
        ),
        # Create the combined file from file two
        "results/combined_file.txt"
```

As you can see, here we are mixing input files with and without `expand`.

I like to put `rule all` in the main snakefile. I also add lots of comments explaining what files are generated. This just helps others work their way through your code. 

We can now run rule all

```bash
$ snakemake --cores 1
```

```
...
[Mon Mar 15 18:44:19 2021]
localrule all:
    input: results/sample_1_fourth_file.txt, results/sample_2_fourth_file.txt, results/sample_3_fourth_file.txt, results/sample_4_fourth_file.txt, results/combined_file.txt
    jobid: 0
...
```

You can see that now rule all is generating all the output files we want.

### Adding parameters
*This corresponds to the rule `rule_params.snake`*

In addition to having `input`, `output`, and `shell` in your rules you can add extra parameters with a `params` argument. You can then refer to specific paramaters in the shell command wih `{params.param_name}`

```{python}
rule params_example:
    input: "results/{sample}_fourth_file.txt"
    output: "results/{sample}_params_ex.txt"
    params:
        file_text = "this is file text",
        other_file = "/beevol/home/wellskri/Analysis/Tutorials/snakemake/results/combined_file.txt"
    shell:
        """
        cp {params.other_file} {output}
        echo {params.file_text} >> {output}
        """
```

Let's add this to rule all and run it

```python
rule all:
    input:
        # Create the fourth file for all samples
        expand(
            "results/{sample}_fourth_file.txt",
            sample = samples
        ),
        # Create the combined file from file two
        "results/combined_file.txt",
        # Create a file using params
        expand(
            "results/{sample}_params_ex.txt",
            sample = samples
            )
```

Testing it

```bash
$ snakemake -np
```

```
Building DAG of jobs...
Job counts:
    count   jobs
    1   all
    4   params_example
    5
[Mon Mar 15 18:49:05 2021]
rule params_example:
    input: results/sample_3_fourth_file.txt
    output: results/sample_3_params_ex.txt
    jobid: 16
    wildcards: sample=sample_3
        cp /beevol/home/wellskri/Analysis/Tutorials/snakemake/results/combined_file.txt results/sample_3_params_ex.txt
        echo this is file text >> results/sample_3_params_ex.txt
```

Here you can see that it has replaced the `{params.name}` with the data contained by the parameter data.

Now to run it and check the output
```{bash}
$ snakemake --cores 1
$ cd results
$ more sample_1_params_ex.txt
```

```
sample_1
sample_2
sample_3
sample_4
this is file text
```

While the params argument is a bit silly in this example, I often use it to pass arguments like the gtf file or other alignment parameters I want more control over and don't want to hard code into the alignment command.

### Adding a log file
We can also add a log file which will save all of the output and error information with the prefix that we chose once we submit to a cluster. It is a simple additional argument

```python
rule log_example:
    input: "results/{sample}_fourth_file.txt"
    output: "results/{sample}_log_ex.txt"
    log: "results/log/{sample}_log"
    params:
        file_text = "this is file text",
        other_file = "/beevol/home/wellskri/Analysis/Tutorials/snakemake/results/combined_file.txt"
    shell:
        """
        cp {params.other_file} {output}
        echo {params.file_text} >> {output}
        """
```

One thing to notice is that the log file path must contain the same wildcards as the rest of the rule.

Let's add this to rule all and run it

```python
rule all:
    input:
        # Create the fourth file for all samples
        expand(
            "results/{sample}_fourth_file.txt",
            sample = samples
        ),
        
        # Create the combined file from file two
        "results/combined_file.txt",
        
        # Create a file using params
        expand(
            "results/{sample}_params_ex.txt",
            sample = samples
            ),
        # Create a log file
        expand(
            "results/{sample}_log_ex.txt",
            sample = samples
            )
```

Let's test running it. I always test with a dry run

```bash
$ snakemake -np
```
```
Building DAG of jobs...
```

```bash
$ snakemake --cores 1
$ cd results/
$ ls
```
```
combined_file.txt         sample_1_params_ex.txt    sample_2_params_ex.txt    sample_3_params_ex.txt    sample_4_params_ex.txt
log                       sample_1_second_file.txt  sample_2_second_file.txt  sample_3_second_file.txt  sample_4_second_file.txt
sample_1_first_file.txt   sample_2_first_file.txt   sample_3_first_file.txt   sample_4_first_file.txt
sample_1_fourth_file.txt  sample_2_fourth_file.txt  sample_3_fourth_file.txt  sample_4_fourth_file.txt
sample_1_log_ex.txt       sample_2_log_ex.txt       sample_3_log_ex.txt       sample_4_log_ex.txt
```

A new directory has shown up called "log". Snakemake created this for us because this directory was specified under the `log:` argument. 

```bash
$ cd log/
$ ls
```

```
sample_1_log  sample_2_log  sample_3_log  sample_4_log
```

At the moment these files are empty, but this will become more helpful when we submit the job to the cluster.

## The config file
The way our pipeline is currently written, if you want to change the samples, you need to open the `Snakefile` and find where that is defined. This can be a bit tricky, and you don't necessarily want people directly modifying the `Snakefile`. A great alternative to this is the addition of a `configfile`. This is a file written in `yaml` format that can be modified by the user and then read into snakemake. The config file can contain anything you want. I most often use it for sample names and paths to files (like the genome reference for alignment). I sometimes also use it for parameters that may be changed (like the adaptor sequences for trimming).

You can either include the `configfile` in the `Snakefile` or you can include it as an argument to `snakemake` (ex `--configfile config.yaml`). The benefit of the second option is if you wanted to run the same snakemake pipeline with different options to see the output. You can keep both config files and just provide the different files in arguemnts to `snakemake`.

`YAML` syntax is relatively straight forward. The full document is read into `snakemake` as a dictionary. Each key in the dicitonary is defined in the following way:

```yaml
DATA_PATH:
```

Here, the key just needs to end with `:`.

I will write a quick config file here and show you how it is interpreted by `snakemake`

```yaml
DATA_PATH:
  /beevol/home/wellskri/Analysis/Tutorials/
SAMPLES:
  - sample1
  - sample2
  - sample3
  - sample4
GROUPS:
  group1:
    - sample1
    - sample2
  group2:
    - sample3
    - sample4
```

If we now update the `Snakefile` like so:
```python
print(config)
```

And run it:

```bash
$ snakemake -np --configfile config.yaml
```

Notice, we need to now add the `-configfile` flag

We get the following output from our print function.
```
{'DATA_PATH': '/beevol/home/wellskri/Analysis/Tutorials/', 'SAMPLES': ['sample_5', 'sample_6', 'sample_7', 'sample_8'], 'GROUPS': OrderedDict([('group1', ['sample_5', 'sample_6']), ('group2', ['sample_7', 'sample_8'])])}
```

Notice how the `config` variable is now a dictionary with a variety of value types. We can walk through the individually.

First, we have the `DATA_PATH` variable
```yaml
DATA_PATH:
  /beevol/home/wellskri/Analysis/Tutorials/
```

Is read in `snakemake` as a string.

```
'DATA_PATH': '/beevol/home/wellskri/Analysis/Tutorials/'
```

We can pull out this string by using the normal dictionary syntax of `dictionary[key]`

```python
DATA_PATH = config["DATA_PATH"]
```

Next, we have the `SAMPLES` variable
```yaml
SAMPLES:
  - sample_5
  - sample_6
  - sample_7
  - sample_8
```

Is read into `snakemake` as a list.

```
'SAMPLES': ['sample_5', 'sample_6', 'sample_7', 'sample_8']
```

We can pull out this string by using the normal dictionary syntax of `dictionary[key]`

```python
SAMPLE_LIST = config["SAMPLES"]
```

```yaml
GROUPS:
  group1:
    - sample_5
    - sample_6
  group2:
    - sample_7
    - sample_8
```

Is read into `snakemake` as a dictionary.
```
'GROUPS': OrderedDict([('group1', ['sample_5', 'sample_6']), ('group2', ['sample_7', 'sample_8'])])
```

We can pull out this string by using the normal dictionary syntax of `dictionary[key]`

```python
GROUPS = config["GROUPS"]
```

Let's update our Snakefile to read in all of these variables.

```python
SAMPLES_TWO = config["SAMPLES"]
DATA_PATH = config["DATA_PATH"]
GROUPS = config["GROUPS"]
```

And let's update our `rule_all` to include these new samples

```python
rule all:
    input:
        # Create the fourth file for all samples
        expand(
            "results/{sample}_fourth_file.txt",
            sample = samples
        ),
        
        # Create the combined file from file two
        "results/combined_file.txt",
        
        # Create a file using params
        expand(
            "results/{sample}_params_ex.txt",
            sample = samples
            ),
        # Create a log file
        expand(
            "results/{sample}_log_ex.txt",
            sample = samples
            ),
        # Trying with samples from the config file
        expand(
            "results/{sample}_params_ex.txt",
            sample = SAMPLES_TWO
        )
```

```bash
$ snakemake -np --configfile config.yaml
```

```
Building DAG of jobs...
Job counts:
    count   jobs
    1   all
    4   four
    4   make_file_general
    4   make_second_file
    4   params_example
    17
[Tue May  4 17:04:57 2021]
rule make_file_general:
    output: results/sample_5_first_file.txt
    jobid: 24
    wildcards: sample=sample_5
...
```

Again, I won't show you all here, but notice that the wildcard is now `sample_5` which came from the `config.yaml` document.  Also, because I just specified `results/{sample}_params_ex.txt` Snakemake checks for each of the necessary input files and makes all files necessasry (all of them for samples 5-8)

Quickly running

```bash
$ snakemake --configfile config.yaml --cores 1
$ cd results/
$ ls
```

```
combined_file.txt         sample_2_fourth_file.txt  sample_3_second_file.txt  sample_5_params_ex.txt    sample_7_params_ex.txt
log                       sample_2_log_ex.txt       sample_4_first_file.txt   sample_5_second_file.txt  sample_7_second_file.txt
sample_1_first_file.txt   sample_2_params_ex.txt    sample_4_fourth_file.txt  sample_6_first_file.txt   sample_8_first_file.txt
sample_1_fourth_file.txt  sample_2_second_file.txt  sample_4_log_ex.txt       sample_6_fourth_file.txt  sample_8_fourth_file.txt
sample_1_log_ex.txt       sample_3_first_file.txt   sample_4_params_ex.txt    sample_6_params_ex.txt    sample_8_params_ex.txt
sample_1_params_ex.txt    sample_3_fourth_file.txt  sample_4_second_file.txt  sample_6_second_file.txt  sample_8_second_file.txt
sample_1_second_file.txt  sample_3_log_ex.txt       sample_5_first_file.txt   sample_7_first_file.txt
sample_2_first_file.txt   sample_3_params_ex.txt    sample_5_fourth_file.txt  sample_7_fourth_file.txt
```

We successfully made output files for all of our samples now.

Now if you wanted to change the samples, instead of changing the `Snakefile` you can just change the `config.yaml`. It is good practice to comment your `configfile` so that it is clear to someone changing it what each argument is for.

## More complicated rules
We can add even more complexity to the rules by including more `python` (and even `R`!) 

### Python code as the rule
*This corresponds to the rule `rule_python.snake`*

Up until this point, we have included only shell scripts as the rules, but we can also write a rule that runs some python code. To do this, you just change the `shell:` argument to `run:`

```python
rule python_example:
    input: "results/{sample}_params_ex.txt"
    output:
        python_output = "results/{sample}_python_code.txt",
        shell_output = "results/{sample}_python_shell.txt"
    params:
        sample_list = SAMPLES_TWO
    run:
        # Some code
```

Now you can run any python code you want using any of the variables defined. You can easily access the variables in the python script.

* To access input files
    * By index: `input[0]` or any index that exists
    * By name: `input.name`
* To access output files
    * By index: `output[0]`
    * By name: `output.python_output`
* To access parameters
    * By index: `params[0]`
    * By name: `params.sample_list`

We can now use this to write a simple python script

```python
rule python_example:
    input: "results/{sample}_params_ex.txt"
    output:
        python_output = "results/{sample}_python_code.txt",
        shell_output = "results/{sample}_python_shell.txt"
    params:
        sample_list = SAMPLES_TWO
    run:
        # We can pull any data from the above, just like with shell
        with open(output.python_output, "w") as out_file:
            with open(input[0], "r") as in_file:
                for line in in_file:
                    out_file.write(line)
            for sample in params.sample_list:
                out_file.write(sample + "\n")
```

You can also run a shell command from your python script.

```python
rule python_example:
    input: "results/{sample}_params_ex.txt"
    output:
        python_output = "results/{sample}_python_code.txt",
        shell_output = "results/{sample}_python_shell.txt"
    params:
        sample_list = SAMPLES_TWO
    run:
        # We can pull any data from the above, just like with shell
        with open(output.python_output, "w") as out_file:
            with open(input[0], "r") as in_file:
                for line in in_file:
                    out_file.write(line)
            for sample in params.sample_list:
                out_file.write(sample + "\n")
        # We can also run the shell here
        shell(
            """
            cp {output.python_output} {output.shell_output}
            """
            )
```

This example may seem a little silly, but let's say you want the shell command to be slightly different for different circumstances (like Paired End sequencing vs Single End). You can include a T/F test as the python script and then run the shell command normally. Here's an example from one of my Snakefiles for running RNA-seq analysis.

```python
rule cutadapt:
    input:
        input_list   = _get_input,
        fastqc_output = "{results}/fastqc_pre_trim/fastqc_{sample}_summary_untrimmed.txt"
    output:
        "{results}/cutadapt_trim/{sample}.txt"
    params:
        settings  = CMD_PARAMS["cutadapt"],
        fastq1    = "{results}/cutadapt_trim/{sample}_R1_trimmed.fastq.gz",
        fastq2    = "{results}/cutadapt_trim/{sample}_R2_trimmed.fastq.gz"
    log:
        "{results}/logs/trimming/cutadapt_{sample}"
    message:
        "Trimming illumina adapters for {wildcards.sample}" 
    threads:
        1
    run:
        if IS_PAIRED:
            shell(
                """
                cutadapt \
                    {params.settings} \
                    -o {params.fastq1} \
                    -p {params.fastq2} \
                    {input.input_list[0]}  \
                    {input.input_list[1]}
                """
            )
        else:
            shell(
                """
                cutadapt \
                    {params.settings} \
                    -o {output} \
                    {input.input_list[0]} 
                """
            )
        with open(output[0], "w") as out:
            out.write("done with cutadapt trimming/n")
```

Here, the shell command is slightly different for paired end and single end. `IS_PAIRED` is inferred based on if two fastq files were provided for the sample so this will automatically run the correct command based on the input provided by the user.

There are also a few other additions in this attached file. I've included a `message`, which is simply what will be printed to the console or the log file about what step is being completed. One more addition is "settings" under `params` which links to a line in the config where `cutadapt` settings can be specified. Finally, I've included an input function which I will expalin shortly.

### Running a script as the rule
Instead of including a `script` written under the `run` variable, we can include a full script as part of the pipeline. These scripts allow us to write a much more complicated file and include is in our pipeline.

When writing a python script, you can call any snakemake parameters the same way you did in the `run` section but also including `snakemake`. For example, access input files with `snakemake.input[0]`. I like to name all variables at the top of the script. So we can write a rule as follows. Notice that we have replaced `run` with `script` and included a path to the script.

```python
rule python_script:
    input: "results/{sample}_python_code.txt"
    output:
        python_output = "results/{sample}_python_script.txt"
    params:
        sample_list = SAMPLES_TWO
    script:
        "../scripts/python_example.py"
```

Snakemake will look for scripts based on the location of the rule file so we need to go up one directory to get to the scripts directory.

Below is the contents of `python_example.py`

```python
# I first set all variables
output_file = snakemake.output.python_output
input_file = snakemake.input[0]
samples = snakemake.params.sample_list
with open(output_file, "w") as out_file:
    with open(input_file, "r") as in_file:
        for line in in_file:
           out_file.write(line)
        for sample in samples:
            out_file.write(sample + "\n")
```

This doesn't only work for only python scripts, it also works for R. To access the snakemake variables in R, the syntax is slightly different.

```R
input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
interesting_parameter <- snakemake@params[[1]]
# The parameters can also be named
input_matrix <- snakemake@input[["matrix"]]
input_gtf <- snakemake@input[["gtf"]]
```

One important thing to note is that the variables are now 1 based and not 0 based.

You can see an R script that I wrote to work with `Snakemake` here: https://github.com/kwells4/mtec_analysis/blob/master/scripts/figures.R

The whole Snakefile is here: https://github.com/kwells4/mtec_analysis/blob/master/Snakefile

*This was the first Snakefile I wrote. Feel free to go through it. It is a little more simple than the other pipelines I will share with you. It does not have separate rule files (which I've found to be helpful) and it includes a conda environment for each rule. There were many ways that the conda environments were very useful, but I won't go into it here. One other important note, in the final_analysis rule I have the path to the script as a parameter, this is no longer necessary in newer versions of Snakemake, you can now include wildcards in the path to scripts*

### Using a function for input (or paramaters)
*This corresponds to the file `functions.snake`*

Instead of writing out a list of input files, we may want a more complicated input (or parameter). To do this, we can include a function. Let's say that we want to run one rule twice, but the input files will be slightly different. For example, you may want to do something with the output files from both `python_example` and `python_script`. You could write a rule that looks like this

```python
rule function_example:
    input: _get_input
    output: "results/{run_number}/output_{sample}.txt"
    shell:
        """
        echo {input} > {output}
        """
```

Notice that we now have two wildcards: `run_number` and `sample`. We also have a function `_get_input` instead of a file. Here is what that function could look like.

```python
def _get_input(wildcards):
    if wildcards.run_number == "first_run":
        return("results/" + wildcards.sample + "_python_code.txt")
    elif wildcards.run_number == "second_run":
        return("results/" + wildcards.sample + "_python_script.txt")
```

Here, our arguments to the function are `wildcards`. This contains a `wildcards` object from which you can extract any of the `wildcards` variables. For example, we can extract the `run_number` by using `wildcards.run_number`. We can also extract the `sample` with `wildcards.sample`. This will work for any number of `wildcards`.

So the whole `functions.snake` Now looks like this
```python
def _get_input(wildcards):
    if wildcards.run_number == "first_run":
        return("results/" + wildcards.sample + "_python_code.txt")
    elif wildcards.run_number == "second_run":
        return("results/" + wildcards.sample + "_python_script.txt")
rule function_example:
    input: _get_input
    output: "results/{run_number}/output_{sample}.txt"
    shell:
        """
        echo {input} > {output}
        """
```

We can now run this rule with either set of input files by updating the `rule_all`. For example, to run the rule on the output of `python_example` we can all this to rule all

```python
expand(
    "results/{run_number}/output_{sample}.txt",
    run_number = "first_run", sample = SAMPLES_TWO
    )
```

So `rule_all` now is

```python
rule all:
    input:
        # Create the fourth file for all samples
        expand(
            "results/{sample}_fourth_file.txt",
            sample = samples
        ),
        
        # Create the combined file from file two
        "results/combined_file.txt",
        
        # Create a file using params
        expand(
            "results/{sample}_params_ex.txt",
            sample = samples
            ),
        # Create a log file
        expand(
            "results/{sample}_log_ex.txt",
            sample = samples
            ),
        # Trying with samples from the config file
        expand(
            "results/{sample}_params_ex.txt",
            sample = SAMPLES_TWO
        ),
        # Running python code
        expand(
            "results/{sample}_python_code.txt",
            sample = SAMPLES_TWO
            ),
        # Running python script
        expand(
            "results/{sample}_python_script.txt",
            sample = SAMPLES_TWO
            ),
        expand(
            "results/{run_number}/output_{sample}.txt",
            run_number = "first_run", sample = SAMPLES_TWO
            )
```

Now we can try running it
```bash
$ snakemake -np --configfile config.yaml
```

```
Building DAG of jobs...
Job counts:
    count   jobs
    1   all
    4   function_example
    5
[Thu May 27 15:11:31 2021]
rule function_example:
    input: results/sample_7_python_code.txt
    output: results/first_run/output_sample_7.txt
    jobid: 48
    wildcards: run_number=first_run, sample=sample_7
        echo results/sample_7_python_code.txt > results/first_run/output_sample_7.txt
```

Again, I'll only show one, but there's a couple of things I want you to notice. For the rule `function_example`, it shows two `wildcards` now. This is because there are two `wildcards` in the output file. Also, while this rule did not have an `input` file (remember, `input` just had a function `get_input`), we now see an `input` file listed "results/sample_7_python_code.txt". This was returned by the function. Now, even though we did not explictely write any input, `Snakemake` will still make sure the files returned by the function exist before running the rule.

One other thing to note here, I've change the output to be its own directory. This is something that I often do to make my results cleaner. But as long as this directory is in the `output` argument, `Snakemake` will make the directory for you. If the direcotory is only in the parameters, you will need to make the direcotry yourself.

Let's also add the second option to rule all.

```python
rule all:
    input:
        # Create the fourth file for all samples
        expand(
            "results/{sample}_fourth_file.txt",
            sample = samples
        ),
        
        # Create the combined file from file two
        "results/combined_file.txt",
        
        # Create a file using params
        expand(
            "results/{sample}_params_ex.txt",
            sample = samples
            ),
        # Create a log file
        expand(
            "results/{sample}_log_ex.txt",
            sample = samples
            ),
        # Trying with samples from the config file
        expand(
            "results/{sample}_params_ex.txt",
            sample = SAMPLES_TWO
        ),
        # Running python code
        expand(
            "results/{sample}_python_code.txt",
            sample = SAMPLES_TWO
            ),
        # Running python script
        expand(
            "results/{sample}_python_script.txt",
            sample = SAMPLES_TWO
            ),
        # Test input functions
        expand(
            "results/{run_number}/output_{sample}.txt",
            run_number = "first_run", sample = SAMPLES_TWO
            ),
        expand(
            "results/{run_number}/output_{sample}.txt",
            run_number = "second_run", sample = SAMPLES_TWO
            )
```

And we can test again

```bash
$ snakemake -np --configfile config.yaml
```

```
Building DAG of jobs...
Job counts:
    count   jobs
    1   all
    8   function_example
    9
[Thu May 27 15:25:46 2021]
rule function_example:
    input: results/sample_5_python_script.txt
    output: results/second_run/output_sample_5.txt
    jobid: 50
    wildcards: run_number=second_run, sample=sample_5
        echo results/sample_5_python_script.txt > results/second_run/output_sample_5.txt
...
[Thu May 27 15:25:46 2021]
rule function_example:
    input: results/sample_5_python_code.txt
    output: results/first_run/output_sample_5.txt
    jobid: 46
    wildcards: run_number=first_run, sample=sample_5
        echo results/sample_5_python_code.txt > results/first_run/output_sample_5.txt
```

I've pulled out two examples of the output here. Both show that the rule being run is `function_example`, but their input files are different because the input file is determined from the function we wrote.

One important note here, there is another, more clear, way to write this particular rule (although we wrote it with a function for example purposes). Instead of writing the rule to need an input function, we can easily write this rule to just need a second `wildcard`.

```python
rule multi_input_example:
    input: "results/{sample}_python_{run_type}.txt"
    output: "results/{run_type}/output_{sample}.txt"
    shell:
        """
        echo {input} > {output}
        """
```

Now in `rule_all`
```python
rule_all:
        # Test input functions
        expand(
            "results/{run_type}/output_{sample}.txt",
            run_type = "code", sample = SAMPLES_TWO
            ),
        expand(
            "results/{run_type}/output_{sample}.txt",
            run_type = "script", sample = SAMPLES_TWO
            )
```

This will get us the same result, but the output files will be named in a more consistent way to the input files. *Note: we can't have both `function_example` and `multi_input_example` in one pipeline because they have identical output files, `Snakemake` won't know which one to run* While we could easily rewrite this rule to not need an input function, there are many times input functions are very helpful. Here's my example from `cutadapt` again.

```python
# Function to return paths of input files
def _get_input(wildcards):
    # Grab path of the fastq file
    fastq1 = SAMPLE_LIST.loc[wildcards.sample, "fastq1"]
    data_dir = SAMPLE_LIST.loc[wildcards.sample, "data_dir"]
    fastq1 = data_dir + "/" + fastq1
    # Make sure file exists
    fastq1  = _check_path(fastq1)
    if IS_PAIRED:
        # Grab path of second read
        fastq2 = SAMPLE_LIST.loc[wildcards.sample, "fastq2"]
        fastq2 = data_dir + "/" + fastq2
        # Make sure file exists
        fastq2 = _check_path(fastq2)
        return(fastq1, fastq2)
    else:
        return(fastq1)
""" Rules for trimming reads with cutadapt """
# Run cutadapt for single-end reads
rule cutadapt:
    input:
        input_list   = _get_input,
        fastqc_output = "{results}/fastqc_pre_trim/fastqc_{sample}_summary_untrimmed.txt"
    output:
        "{results}/cutadapt_trim/{sample}.txt"
    params:
        settings  = CMD_PARAMS["cutadapt"],
        fastq1    = "{results}/cutadapt_trim/{sample}_R1_trimmed.fastq.gz",
        fastq2    = "{results}/cutadapt_trim/{sample}_R2_trimmed.fastq.gz"
    log:
        "{results}/logs/trimming/cutadapt_{sample}"
    message:
        "Trimming illumina adapters for {wildcards.sample}" 
    threads:
        1
    run:
        if IS_PAIRED:
            shell(
                """
                cutadapt \
                    {params.settings} \
                    -o {params.fastq1} \
                    -p {params.fastq2} \
                    {input.input_list[0]}  \
                    {input.input_list[1]}
                """
            )
        else:
            shell(
                """
                cutadapt \
                    {params.settings} \
                    -o {output} \
                    {input.input_list[0]} 
                """
            )
        with open(output[0], "w") as out:
            out.write("done with cutadapt trimming/n")
```

Here, the input files are different if the data is paired end or not. We don't want to try to run the rule if just one of the two fastq files is present in paired end data, but if we don't have paired end data, providing two files as input will crash the pipeline. I also include the output of `fastqc` because I don't want both programs to have the same fastq file open at the same time, so now `cutadapt` will wait to run until `fastqc` has finished.

## Running jobs on the server
*This corresponds to the rule `functions.snake`*

Most importantly, `Snakemake` plays well with servers and can run all of your jobs in parallel. Right now, we are just running all of our jobs one at a time, but if we submit to the server, all samples can be run at the same time. It's pretty simple to do.

First, you need to install drmaa
```bash
$ mamba install drmaa
```

Now we need to add two paramaters to all of our rules.

1. `job_name` The name of the job that apears in the queue
2. `memory` The amount of memory to give the job.

We also need to add a new argument: `threads`

Below is an example of a rule with the parameters
```python
rule submit_example:
    output:
        "results/submit_ouput/{sample}_submit.txt"
    params:
        sample_list = SAMPLES_TWO,
        job_name    = "test_jobs_{sample}",
        memory      = "select[mem>4] rusage[mem=4]"
    threads:
        1
    shell:
        """
        echo "completed" > output
        """
```

Now, I also have a driver script called `snakecharmer.sh`. Here is what that looks like
```bash
#!/usr/bin/env bash
#BSUB -J RNAseq
#BSUB -o logs/snakemake_%J.out
#BSUB -e logs/snakemake_%J.err
#BSUB -R "select[mem>4] rusage[mem=4] " 
#BSUB -q rna
set -o nounset -o pipefail -o errexit -x
# Load modules
module load fastqc/0.11.7
module load samtools/1.5
module load STAR/2.5.2a
module load subread
# LSF arguments
args=' 
  -q rna 
  -o {log}.out 
  -e {log}.err 
  -J {params.job_name} 
  -R "{params.memory} span[hosts=1] " 
  -n {threads} ' 
# Run snakemake pipeline
snakemake \
    --drmaa "$args" \
    --snakefile Snakefile \
    --configfile config.yaml \
    --jobs 60 \
    --latency-wait 60 \
    --rerun-incomplete
```

To walk through this, there are a few things going on.

First, this is a normal submit script, the `#BSUB` arguments probably look familiar to you

Next, we set shell options, these are especially helpful for exiting if something fails and for printing errors.

* -o nounset
    * Treat unset variables and parameters other than the special parameters `@` or `*` as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit.
* -o pipefail
    * If set, the return value of a pipeline is the value of the last (rightmost) command to exit with a non-zero status, or zero if all commands in the pipeline exit successfully. This option is disabled by default.
* -o errexit
    * exit if any command returns non zero exit status
* -x
    * print trace of commands

Then, we load all of the modules we will need (this is pretty unique to each of your pipelines).

Now is the interesting part, the arguments. Here we are tellling the system to use the arguments from snakemake.

1. Logs (`-o` and `-e`)
    * It will look under the log argument for the base name of the log file. It will then output that name `.err` and `.out` as the output files. I like to organize it so all of my logs are in one `logs` directory and then each rule gets its own directory. Each sample then gets it's own file. If the log does not contain all the same wildcards as the output files, you would not get unique log files for each rule. Because of this, snakemake won't run if the naming isn't unique. The log for my trimming example is `{results}/logs/trimming/cutadapt_{sample}`. Here you see I have multiple subdirectories before finally getting to the file name (based on the sample).
2. job name (`-J`)
    * This tells `Snakemake` to find the `job_name` under `params`.
3. memory (`-R`)
    * This tells `Snakemake` to find the `memory` under `params`
4. Threads (`-n`)
    * This tells `Snakemake` to find the number of threads under the `threads` argument

Finally, the call to snakemake is made. Here `drmaa` provides the arguments listed above to snakemake based on each specific rule. We also provide `Snakemake` with the name of the `snakefile` and `configfile`. These are easy to change if you want to test out different parameters in the `configfile`. The last three arguments are telling snakemake how many jobs to submit (60 in this case), how long to wait to see completed files (60 seconds, this is just to account for any lag associated with submitting jobs), and to rerun any jobs that were not completed. Other arguments that may be helpful include the `--keep-going` argument. This is helpful because if one job fails, by default `Snakemake` won't submit any new jobs. Sometimes this is nice (remember, we can't restart a `Snakemake` run until all jobs have finished), but this is sometimes annoying (like if you left it to run overnight, but one sample failed on step one). If you'd rather `Snakemake` keep submitting jobs even if one failed, you can include the `--keep-going` flag.

Once this file is written, you can easily submit with
```bash
$ bsub < snakecharmer.sh
```

When this happens, there are a few things `Snakemake` will do. It will first check all of the output files listed under `rule_all`. It will then back track through your rules until it finds a rule (upsteam of `rule_all`) where the input files exist. It will now submit that rule for each sample as an individual job. When that job completes for each sample, it will submit the next job. It will keep doing this until it makes all of the output files requested by `rule_all`. Because each sample is now idependent, as each task finishes for the sample, a new one will start. 

## Some headaches and errors
1. Can't determine wildcards
    * This often happens if you have a case like this for your output `some/path/to/{sample}_{run}`. If your sample is `WT_1` and your run is `run_1`, `Snakemake` can't determine what the `wildcards` are. For example is it sample = `WT_1_run` and run is `1`? or `WT_1` and `run_1` or `WT` and `1_run_1`? I often get around this by just adding some more specificity. For example, I'd change the output to `some/path/to/{sample}_run_{run}` and then have my sample as `WT_1` and the run as `1`. Notice here, there is no way to mix up the two `wildcards`.
2. Can't find input files
    * This almost always means I have a typo somewhere between `rule_all`, the rule that makes the file, or the rule(s) that use the file
3. Output files don't exist
    * This happens most often when I am creating a file specified in the `params` not `output` variables. If there is a new directory, it will only be created by `Snakemake` if it's specified in the `output`. If this happens and I want to keep the one output in the `params`, I can instead add a directory to the `output` params (which is easy, just do `directory(path/to/directory`))
4. The snakefile is locked.
    * Unfortunately (or fortunately), you can only run instance of `Snakemake` in a directory at one time. This is a good thing because you can't overwrite files that are in progress. If you try to do this, `Snakemake` will lock the directory. Now when you try to run  again, you will get an error that the directory is locked (even if `Snakemake` is no longer actively running). It is easy to fix though. First make sure there is no instance of `Snakemake` running (check all of your jobs and windows, my biggest problem is if one job didn't finish running and I didn't notice). Once you've checked run `Snakemake --configifile config.yam --unlock`. Now everything should work again.

There are many more errors that can happen, these are just the most common (and honestly, they happen to me a lot, it's why I always test first with `-np`). Let me know if you run into some problem you don't know how to fix. Either I've seen it or I may have some idea how to fix it.
