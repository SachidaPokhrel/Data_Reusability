# PBS Pro allocation submission

../submit.py validates qsub as PBS, inspects the requested queue using qstat -Qf, and prepares
qsub resource arguments for one select chunk with CPUs/memory and walltime. Account/queue must
be supplied. Native depend=afterok is supported. Generated scripts/logs are under ignored logs/hpc.
Run `./run.sh pilot --executor pbs --account YOUR_ACCOUNT --queue YOUR_QUEUE --cores 8` from the
root after setup and input validation. This has not been tested on a real PBS installation.
