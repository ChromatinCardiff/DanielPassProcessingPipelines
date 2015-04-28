<b>Job submission to slurm cluster:</b>
```
sbatch slurm_submission.sh
```
<b>Script looks like this:</b>
```
#!/bin/bash
#
#SBATCH --job-name=parsing
#SBATCH --output=slout.txt
#
#SBATCH --ntasks=24
#SBATCH --time=100:00
#
#SBATCH --error=slout.err
#SBATCH --output=slout.out

srun hostname
blast -in abc.fas -out out.txt -a 24
```

<b>Useful commands:</b>
