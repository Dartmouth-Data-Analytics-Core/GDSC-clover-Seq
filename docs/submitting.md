# Submitting the GDSC-Clover-Seq Pipeline

Before submitting the GDSC-Clover-Seq Pipeline ensure you have done the following:

1. In the `prebuilt_config` corresponding to your organism, ensure you have toggled on/off the database build module (default is to NOT build.)

2. In the `prebuilt_config` corresponding to your organism, ensure you have set the `refLevel` to the approproate reference level group for your sample. This needs to match a valid entry in your `Sample_list_SE.txt` file (namely, the Group column).

3. Ensure that in `job.script.sh` you have edited the SBATCH email directive to point to your email. 

4. Ensure that in `job.script.sh` that you have modified line 16 (`CONFIG`) to be one of `hg38`, `mm10`, or `dm6`. This will automatically set the pipeline to run using the appropriate `prebuilt_config`. 

5. Ensure your `Sample_list_SE.txt` has the following 3 column names, comma-separated and case-sensitive: `Sample_ID,fastq_1,Group`.

6. The paths in your `Sample_list_SE.txt` file are accessibile. Your data can either live in the working directory, or elsewhere on discovery, but full file paths are needed if the later option is specified. 

Once this is all finished, you can submit the pipeline using the following command:

```shell
sbatch job.script.sh
```

The job will spawn a log file called `clvSeq_xxxxxxx.log` and create a folder called `slurm_logs` where individual rule logs will populate as the pipeline runs. Error logs in this folder are often very informative and should be utilized if the pipeline fails. 

To check the status of your job, run the `squeue` command and grep for your Dartmouth NetID.

```shell

# Replace f007qps with your netID
squeue | grep "f007qps"

```