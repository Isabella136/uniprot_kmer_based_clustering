Requires diamond to be installed; if conda or mamba is installed, this can be done via the command
```
[conda or mamba] env create -f diamond.yaml
```
Requires rust to be installed; click [here](https://www.rust-lang.org/tools/install) for download and installation instructions.

Also requires rust nightly to be enabled via the following command:
```
rustup default nightly
```
The actual program can be run as an sbatch job via the command below, but users may need to update some of the job specifications (notably the account and partition if outside CBCB)
```
sbatch run.sh
```
Alternatively, the program can also be run locally via the command
```
cargo run --release -- uniprot_arg.fasta [threads]
```
