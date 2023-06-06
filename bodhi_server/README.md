# How to run jobs on a server

## Key points

* Connect using `ssh`
* Always immediately start an interactive session using `qlogin -R "rusage[mem=10]"`
* Submit jobs with `bsub`
* Submit bash scripts with `bsub < test.sh`
* The server does not have a great backup system and sometimes goes down, always have another backup of your data and scripts to make your final (github is a great place to store scripts)
* Edit files using `sshfs` to edit locally or `vim`/`emacs` to edit on the server.

## Server basics

* A server is composed of many compute nodes. You can think of these as individual computers that make up the cluster (or server). Each of these nodes has different resources that range in the number of cpus and memory.
* Servers are great because you can run multiple jobs at once on a single compute node or multiple compute nodes.
* LSF scheduler
    * To run jobs, Bodhi uses an lsf scheduler (although other schedulers exist like slurm). These job schedulers attempt to fairly distribute resources between individuals rather than allowing one person to hog all of the resources.
    * The scheduler will decide a compute node to run your job (you don't have to do this!) and will start running the job on that node.
    * To submit a job, you must estimate the resources you need so the scheduler knows what compute node to send your job to (based on what resources are currently available and where). If you end up needing more resources than you requesed, your job will be terminated and you will need to request resources again.
    * When you submit a job to the scheduler, depending on how many other jobs you and others have running, your job will either start running immediately or be "pending" until there is enough space to run a new job. How quickly your job starts will also depend on how many resources you think you will need to run your job.
    * When the job is completed output files will appear that show you any output from the commands run. 

## Connecting to the server

* If you are off campus, you will first need to connect to the VPN
* Once on the VPN, open up a terminal window
* Connect to the server using your username

```{bash}
ssh your_username@amc-bodhi.ucdenver.pvt
```

There will likely be a message about the host. Type yes

```
The authenticity of host '[hostname] ([IP address])' can't be established.
RSA key fingerprint is [key fingerprint].
Are you sure you want to continue connecting (yes/no)? yes
```

You computer may also ask about adding to known hosts. Type y and you won't be asked again.

You will also be asked to give the password that David provided. You have three attempts before you have to rerun the ssh command.

## The interactive session

Once connected, you will be on the head node or the login node. This is one of the many compute nodes on bodhi, but it is not intended for computing. There are very few resources available on this node, so you can quickly use up all of the resources. To be able to interactively compute, you can use an interactive session. To start an interactive session use `qlogin`

```{bash}
$ qlogin
```

```
Interactive
******************************************
* You did not specify a memory resource  *
* A default value of 4G has been applied *
* If your application exceeds this value *
* Your job will be killed                *
******************************************
Job <641243> is submitted to queue <interactive>.
<<Waiting for dispatch ...>>
<<Starting on compute00>>
```

This tells you many things. First, just typing `qlogin` will start an interactive session, but it gives a warning. This warning says that your interactive job will be given 4G of memory. It also says if you need more than this your job will be killed, so no matter what you are doing, the job will be immediately terminated and anything that wasn't saved will be lost. Finally, it tells you the name of the job `<641243>`, the queue `<interactive>`, and the compute node `compute00`.

Notice that the path before the `$` has now changed from `[wellskri@amc-bodhi ~]` to `[wellskri@compute00 ~]`

Now you can look to see if this job is running by typing `bjobs`

```{bash}
$ bjobs
```

```
JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
641243  wellskr RUN   interactiv amc-bodhi.c compute00   bash -l    Mar  4 13:55
```

This also tells you the job name, queue, and compute node. In addition, it tells you "STAT" the status (RUN, PEND), the name, and the submit time.

What if you want more resources? You can ask for more with your qlogin command.

```{bash}
# First end your current interactive session
$ exit
# Then start a new session with more resources
$ qlogin -R "rusage[mem=10]"
```

```
Interactive -R rusage[mem=10]
rusage[mem=10]
MEMORY LIMIT: 10
Job <641188> is submitted to queue <interactive>.
<<Waiting for dispatch ...>>
<<Starting on compute24>>
```

This now tells you that you asked for more memory. You can change this command to be 20, 30, or even 40 G of memory. One note is that as you ask for more resources, it will take longer for your job to start.


## Submitting jobs

For submitting a job and learning other fun cluster resources, I suggest you use the tutorial created by Dave https://github.com/dpastling/lsf_tutorial

General pros of job submission and some tips

* I recommend writing bash scripts and submitting your jobs that way. When you do that, you always have a record of what you did
* When you submit a job, always include paths for the log files, the `-o` and `-e` below. Just make sure you make a log directory or you won't be able to see the output.
```{bash}
#!/usr/bin/env bash
#BSUB -J ShortName
#BSUB -n 2
#BSUB -R "select[mem>1] rusage[mem=1] span[hosts=1]"
#BSUB -o logs/stdout_%J.out
#BSUB -e logs/stderr_%J.err
perl some_script.pl
```
* If you submit a job, you can now shut down your computer, and as long as there were no bugs you can come back tomorrow and it should be done.
* Submitting to a cluster means you can ask for more resources. For example, when using `star`, you can include a `threads` argument to utilize more resources. Also for programs like `star`, setting `-R "select[mem>40] rusage[mem=40]` provides much more memory than you could get on your local computer.
* An additional advantage is that you can submit all of your jobs at once. So if you have 8 samples to align, you can align them at the same time. You can do this by submitting a job array (as described in the tutorial) or by using a program like Snakemake. Snakemake is a very cool tool that allows you to write rules that can be applied to any samples. For example, you can have a rule to run `fastqc`, another to run `star`, another to run `featureCounts` and finally one to create a count table. You only have to start the driver and Snakemake will take care of all the steps for all the samples you specify. If you have any interest in learning how to use Snakemake, let me know and we can do another of these.

## Connecting between your computer and server

One tricky part about working on a server is connecting between your local computer and the server. I have two ways that I recommend to do this.

### Rsync
To copy files between the server and your local computer, I recommend using `rsync`. This will check to see if files with identical sizes already exist and will copy missing files. A benefit is that if a file did not correctly transfer, `rsync` will recopy that file if you try again. To run rsync to the server:

From the server to your local:
```{bash}
# Run this on your local computer running bash
rsync -avz --progress your_username@amc-bodhi.ucdenver.pvt:/beevol/home/your_username/path/to/data/* path/to/local/directory
```

```{zsh}
# Run this on your local computer running zsh, notice zsh requires ""
rsync -avz --progress "your_username@amc-bodhi.ucdenver.pvt:/beevol/home/your_username/path/to/data/*"" path/to/local/directory
```

Here, -a is a shortcut to many other arguments (rlptgoD no -H, -A, -X)
* -r recursive (go into all directories recursively)
* -l copy symlinks as symlinks
* -p preserve permissions
* -t preserve modification times
* -g preserve group
* -o preserve owner
* -D preserve device and special files
* Exclude: -H preserve hard links
* Exclude: -A preserve ACLs (access control list, who gets permission)
* Exclude: -X preserve extended attributes

* -v verbose (print warnings)
* -z compress (not necessary if what you are moving is already compressed)
* --progress  print progress, really nice to have if you are copying many files or large files.

From your local to the server:
```{bash}
# Run this on your local computer running bash
rsync -avz --progress path/to/local/directory/* your_username@amc-bodhi.ucdenver.pvt:/beevol/home/your_username/path/to/directory/on/server/
```

```{zsh}
# Run this on your local computer running zsh, notice zsh requires ""
rsync -avz --progress path/to/local/directory/* "your_username@amc-bodhi.ucdenver.pvt:/beevol/home/your_username/path/to/directory/on/server/"
```

### Sshfs
There are many ways to modify files on the server.
1. You may modify them by using a text editor like `vim` or `emacs` (I use `vim`)
2. You may modify them by using your favorite local text editor and then copying that back to the server (I don't recommend this)
3. You can install macFUSE and sshfs and use sshfs (my prefered method)

#### What is sshfs?
Sshfs allows you to directly connect to the server and mount the file system onto your local computer. This means that you can open files from the server on any of your favority local text editors and it will automatically update the file on the server as soon as you save.

#### How do you set up sshfs?
1. Go to https://osxfuse.github.io/
2. Download both sshfs and macFUSE and install using the instructions

That's it! You should now be able to use sshfs

#### How do you use sshfs?
1. Set up an empty directory to use for sshfs (full disclosure, I started doing this my first year of grad school when I was just learning and haven't looked into it seriously again because it's worked well. There may be better ways to do it). This only needs to be done once.
```{bash}
mkdir Documents/sshfs
```

2. Connect to the server
```{bash}
# You can pass any folder on the server that you want
sshfs your_username@amc-bodhi.ucdenver.pvt:/home/beevol/your_username Documents/sshfs
```

3. Now you can view the files on the server from your local computer either using the terminal or you can even view the files in finder. This is very convienient if you have html output (like from `fastQC`) because you can now view these files without transferring to your local.

4. You can edit files on the server any way you want. I personally like to use `sublime` because it has language specific highlighting
```{bash}
sublime Documents/sshfs/my_server_file.txt
```

5. Before you put your computer to sleep or log off, make sure to eject the connection. I personally go to my documents folder in finder and eject (in case there is something still running), but you can also eject from the terminal
```{bash}
umount Documents/sshfs
```

6. Sometimes if your connection is interrupted or your computer goes to sleep while still connected, the created `sshfs` folder may be broken. If so, first try to umount it. This almost always works. If that fails, you may want to restart your computer.
```{bash}
umount Documents/sshfs
```

## Server backup

The server does not have a great backup system, so always have your raw data backed up somewhere else (like your hard drive). Also make sure your processing scripts are saved somewhere (you don't need to backup intermediate files as long as you have the scripts to recreate them). One great way to back up your scripts is with github (this also make life easy when you publish, all you need to do it make your git repo public). If you are interested in learning github we can do another one of these for that as well.
