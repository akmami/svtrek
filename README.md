# SVTrek

SVTrek is a bioinformatics tool for auditing structural variation (SV) discoveries and can also perform SV discovery independently using long-read analysis.

## Dependencies

- A C compiler (`gcc` or `clang`) and `make`
- `autoconf`, `automake`, and `libtool` (only needed when htslib is built locally, for `autoreconf`)
- `zlib` development headers
- `pkg-config` (used to detect a system htslib)

## Installation

```
# clone the repository
git clone --recursive https://github.com/akmami/SVTrek.git
cd SVTrek

# build dependencies
make install

# compile the program
make
```

If you already cloned without `--recursive`, initialise the submodule first:

```
git submodule update --init --recursive
```

### Choosing the htslib source

`make install` (and `make`) accept a `HTSLIB` variable:

- `make install HTSLIB=auto` — use system htslib if `pkg-config` finds it, otherwise build locally (this is the default).
- `make install HTSLIB=system` — require the system htslib reported by `pkg-config`.
- `make install HTSLIB=local` — always clone + build htslib under `deps/htslib`, ignoring any system copy.

Use the same `HTSLIB=...` value for the subsequent `make` so the compile and link flags match.

### Other make targets

- `make clean` — remove the `svtrek` binary and object files.
- `make distclean` — also remove the locally built `deps/htslib` (a system htslib and the abPOA submodule are left untouched).

## Usage
```
./svtrek [MODE] [OPTIONS]
```

## Program Modes

### `disc`

SV discovery mode on a graph-alignment result. It takes a GFA reference graph, a GAF alignment, and the FASTQ reads.

#### Usage
```
./svtrek disc [-r|--gfa GFA] [-a|--gaf GAF] [-q|--fq FASTQ] [OPTIONS]
```

##### Required Parameters
- `-r, --gfa <GFA>`
  - Reference graph in GFA format.
- `-a, --gaf <GAF>`
  - Graph alignment in GAF format.
- `-q, --fq <FASTQ>`
  - Long reads in FASTQ format.

##### Options
- `-o, --output <filename>`
  - Output filename.
  - **Default:** `svtrek.out`
- `-t <num>`
  - Number of threads to use for processing.
  - **Default:** `4`
- `--verbose`
  - Enables verbose output.
  - **Default:** `false`
- `--consensus-interval-range <num>`
  - Interval that limits the refinement range.
  - **Default:** `500`
- `--consensus-interval <num>`
  - Interval within which locations are treated as the same position.
  - **Default:** `5`
- `--consensus-min-count <num>`
  - Minimum number of elements required for consensus determination.
  - **Default:** `3`

### `audt`

Variation auditing mode: parses a BAM file and validates the variations reported in a VCF file.

#### Usage
```
./svtrek audt [-b|--bam BAM] [-v|--vcf VCF] [OPTIONS]
```

##### Required Parameters
- `-b, --bam <BAM>`
  - BAM file to be processed.
- `-v, --vcf <VCF>`
  - VCF file to be audited.

##### Options
- `-o, --output <filename>`
  - Output filename.
  - **Default:** `svtrek.out`
- `-t <num>`
  - Number of threads to use for processing.
  - **Default:** `4`
- `--verbose`
  - Enables verbose output.
  - **Default:** `false`
- `--wider-interval <num>`
  - Offset interval for the start of the reads (DEL-START).
  - **Default:** `20000`
- `--median-interval <num>`
  - Offset interval for point positions (INS).
  - **Default:** `10000`
- `--narrow-interval <num>`
  - Offset interval for the end of the reads (DEL-END).
  - **Default:** `2000`
- `--consensus-interval-range <num>`
  - Interval that limits the refinement range.
  - **Default:** `500`
- `--consensus-interval <num>`
  - Interval within which reads are treated as being at the same position.
  - **Default:** `5`
- `--consensus-min-count <num>`
  - Minimum number of elements required for consensus determination.
  - **Default:** `3`

## Simulation & benchmarking

Helper scripts for generating a simulated benchmark live under `misc/`:

- `misc/simulate.sh` — simulates SVs/SNPs on a reference (SURVIVOR), generates long reads (Badread), and maps them (minimap2), producing a mutated genome, a gold-standard VCF, and a sorted BAM. Species profiles are in `misc/params/`.
- `misc/fix_survivor_vcf.awk` — normalises SURVIVOR's VCF output into a valid, tab-delimited VCF that `bcftools`/`truvari` accept.
