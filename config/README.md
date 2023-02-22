### Please note:

These are Nextflow configuration files to run workflows across the multiple compute infrastructures that the repo author is affiliated with. These include the high-throughput computing cluster (HTC), maintained by the [University of Wisconsin-Madison's Center for High Throughput Computing](https://chtc.cs.wisc.edu/), the Beartooth cluster maintained by [the University of Wyoming Advanced Research Computing Center (ARCC)](https://www.uwyo.edu/arcc/), and the dhogal and dhogal2 clusters maintained by [David O'Connor's research group at University of Wisconsin-Madison](https://dho.pathology.wisc.edu/). Also included are two config files that will enable workflows to run on either x86 (Intel) or ARM (Apple Silicon) MacOS machines. All config files are specified with the nextflow command line flag `-c`.

These files may or may not be useful to other researchers and should be used and modified with caution.

Additionally, this directory contains a Dockerfile and Singularity recipe file to containerize all software dependencies. Images will be regularly updated and versioned on Docker Hub.