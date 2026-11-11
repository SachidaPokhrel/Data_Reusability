# SLURM allocation submission

../submit.py checks sbatch and sinfo for the supplied partition, then prepares one allocation
with account, partition, CPUs, memory, time and separate %j logs. Native afterok dependencies are
supported. Run `./run.sh pilot --executor slurm --account YOUR_ACCOUNT --queue YOUR_PARTITION
--cores 8` from the root after setup. Generated scripts/logs are ignored. No real SLURM test has run.
