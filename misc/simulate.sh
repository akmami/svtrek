#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# defaults
REF_FASTA=""; PREFIX=""; SPECIES="default"; PARAMS_FILE=""
DEPTH_OVERRIDE=""; OUTDIR="sim_out"; THREADS=8; MIN_CONTIG_OVERRIDE=""
PARAMS_DIR="$SCRIPT_DIR/params"
SEED="${SEED:-42}"
COMPRESS_VCF="${COMPRESS_VCF:-1}"
SURVIVOR_TIMEOUT="${SURVIVOR_TIMEOUT:-600}"   # seconds; guards against endless retry loops

msg() {
    printf "\033[1;34m[%s]\033[0m %s\n" "$(date +%H:%M:%S)" "$*" >&2;
}

die() {
    printf "\033[1;31m[ERROR]\033[0m %s\n" "$*" >&2;
    return 1;
}

usage() {
    sed -n '3,40p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//';
    return "${1:-0}";
}

# parse args
while getopts ":r:p:s:c:d:o:t:m:P:h" opt; do
    case "$opt" in
        r) REF_FASTA="$OPTARG" ;;
        p) PREFIX="$OPTARG" ;;
        s) SPECIES="$OPTARG" ;;
        c) PARAMS_FILE="$OPTARG" ;;
        d) DEPTH_OVERRIDE="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        m) MIN_CONTIG_OVERRIDE="$OPTARG" ;;
        P) PARAMS_DIR="$OPTARG" ;;
        h) usage 0; exit 0 ;;
        \?) die "unknown option -$OPTARG (use -h)" ;;
        :)  die "option -$OPTARG needs an argument" ;;
    esac
done

[[ -n "$REF_FASTA" ]] || { usage 1; exit 1; }
[[ -n "$PREFIX"    ]] || { die "prefix (-p) is required"; exit 1; }
[[ -f "$REF_FASTA" ]] || { die "reference FASTA not found: $REF_FASTA"; exit 1; }

# resolve + load params file
[[ -z "$PARAMS_FILE" ]] && PARAMS_FILE="$PARAMS_DIR/${SPECIES}.params"
[[ -f "$PARAMS_FILE" ]] || { die "params file not found: $PARAMS_FILE"; exit 1; }
msg "Loading profile: $PARAMS_FILE"
# shellcheck disable=SC1090
source "$PARAMS_FILE"
[[ -n "$DEPTH_OVERRIDE" ]] && DEPTH="$DEPTH_OVERRIDE"
MIN_CONTIG="${MIN_CONTIG:-1000000}"
[[ -n "$MIN_CONTIG_OVERRIDE" ]] && MIN_CONTIG="$MIN_CONTIG_OVERRIDE"

# preflight tools
for t in SURVIVOR badread minimap2 samtools bgzip tabix awk; do
    command -v "$t" >/dev/null 2>&1 || { die "'$t' not found on PATH"; exit 1; }
done

mkdir -p "$OUTDIR"
OUT_PREFIX="$OUTDIR/$PREFIX"
PARAM_SV="$OUTDIR/${PREFIX}.survivor.params"
MUT_GENOME="${OUT_PREFIX}.fasta"
READS="${OUT_PREFIX}.reads.fastq.gz"
BAM="${OUT_PREFIX}.sorted.bam"


# index reference
[[ -f "${REF_FASTA}.fai" ]] || { msg "Indexing reference"; samtools faidx "$REF_FASTA"; }

WORK_REF="$REF_FASTA"
n_total=$(wc -l < "${REF_FASTA}.fai")
n_keep=$(awk -v m="$MIN_CONTIG" '$2>=m' "${REF_FASTA}.fai" | wc -l)
[[ "$n_keep" -gt 0 ]] || { die "no contigs >= MIN_CONTIG ($MIN_CONTIG bp) in $REF_FASTA; lower -m"; exit 1; }
if [[ "$n_keep" -lt "$n_total" ]]; then
    WORK_REF="$OUTDIR/${PREFIX}.ref.fa"
    msg "Filtering reference: keeping $n_keep/$n_total contigs >= ${MIN_CONTIG} bp -> $WORK_REF"
    awk -v m="$MIN_CONTIG" '$2>=m{print $1}' "${REF_FASTA}.fai" \
        | xargs samtools faidx "$REF_FASTA" > "$WORK_REF"
    samtools faidx "$WORK_REF"
else
    msg "All $n_total contigs >= ${MIN_CONTIG} bp; using reference as-is"
fi

# measure genome size (on the working reference)
GENOME_BP=$(awk '{s+=$2} END{print s}' "${WORK_REF}.fai")
[[ "$GENOME_BP" -gt 0 ]] || { die "could not determine genome size from ${WORK_REF}.fai"; exit 1; }
GENOME_MB=$(awk -v b="$GENOME_BP" 'BEGIN{printf "%.3f", b/1000000}')

# rate * genome_Mb -> integer count (>=0)
countf() { awk -v r="$1" -v mb="$GENOME_MB" 'BEGIN{ c=r*mb; if(c<0)c=0; printf "%d", c+0.5 }'; }

INDEL_N=$(countf "$INDEL_PER_MB")
DUP_N=$(countf "$DUP_PER_MB")
INV_N=$(countf "$INV_PER_MB")
TRA_N=$(countf "$TRA_PER_MB")
INVDEL_N=$(countf "$INVDEL_PER_MB")
INVDUP_N=$(countf "$INVDUP_PER_MB")
SNP_FREQ=$(awk -v s="$SNP_PER_KB" 'BEGIN{printf "%.8f", s/1000}')

msg "Genome: ${GENOME_BP} bp (${GENOME_MB} Mb) -> counts: INDEL=$INDEL_N DUP=$DUP_N INV=$INV_N TRA=$TRA_N INVdel=$INVDEL_N INVdup=$INVDUP_N ; SNP_freq=$SNP_FREQ"

# SURVIVOR parameter file + run
cat > "$PARAM_SV" <<EOF
PARAMETER FILE: DO JUST MODIFY THE VALUES AND KEEP THE SPACES!
DUPLICATION_minimum_length: $DUP_MIN
DUPLICATION_maximum_length: $DUP_MAX
DUPLICATION_maximum_num: $DUP_MAX_COPIES
DUPLICATION_number: $DUP_N
INDEL_minimum_length: $INDEL_MIN
INDEL_maximum_length: $INDEL_MAX
INDEL_number: $INDEL_N
TRANSLOCATION_minimum_length: $TRA_MIN
TRANSLOCATION_maximum_length: $TRA_MAX
TRANSLOCATION_number: $TRA_N
INVERSION_minimum_length: $INV_MIN
INVERSION_maximum_length: $INV_MAX
INVERSION_number: $INV_N
INV_del_minimum_length: $INVDEL_MIN
INV_del_maximum_length: $INVDEL_MAX
INV_del_number: $INVDEL_N
INV_dup_minimum_length: $INVDUP_MIN
INV_dup_maximum_length: $INVDUP_MAX
INV_dup_number: $INVDUP_N
Number_haploid: $N_HAPLOID
homozygous_ratio: $HOM_RATIO
EOF

msg "SURVIVOR simSV: mutated genome + gold VCF"
rm -f "$MUT_GENOME" "${OUT_PREFIX}.vcf"
TO=""; command -v timeout >/dev/null 2>&1 && TO="timeout ${SURVIVOR_TIMEOUT}"
SV_TRIES="${SV_TRIES:-12}"
rc=1
for try in $(seq 1 "$SV_TRIES"); do
    SV_SEED=$(( (SEED + try*104729) % 2000000000 )) \
        $TO SURVIVOR simSV "$WORK_REF" "$PARAM_SV" "$SNP_FREQ" 0 "$OUT_PREFIX" >/dev/null 2>&1
    if [[ -s "$MUT_GENOME" && -s "${OUT_PREFIX}.vcf" ]]; then rc=0; break; fi
    rm -f "$MUT_GENOME" "${OUT_PREFIX}.vcf"
    msg "  SURVIVOR attempt $try/$SV_TRIES failed to place SVs; retrying with a new seed"
    sleep 1.1   # ensures time()-based seed advances on the stock binary
done
[[ $rc -eq 0 ]] || { die "SURVIVOR could not place SVs in $SV_TRIES tries. The usable genome is too small for the requested SV load: lower *_MAX / *_PER_MB, or raise -m so fewer contigs are dropped."; exit 1; }
[[ -s "$MUT_GENOME" && -s "${OUT_PREFIX}.vcf" ]] || {
    die "SURVIVOR produced no output. Raise -m (currently $MIN_CONTIG) or reduce SV sizes/counts in $PARAMS_FILE.";
    exit 1;
}

GOLD_VCF="${OUT_PREFIX}.vcf"
if [[ "$COMPRESS_VCF" == "1" ]]; then
    msg "    sort + bgzip + tabix gold VCF"
    if command -v bcftools >/dev/null 2>&1; then
        # bcftools sort "${OUT_PREFIX}.vcf" -Oz -o "${OUT_PREFIX}.gold.vcf.gz" 2>/dev/null
        awk -f "$SCRIPT_DIR/fix_survivor_vcf.awk" "${OUT_PREFIX}.vcf" | bcftools sort -Oz -o "${OUT_PREFIX}.gold.vcf.gz" 2>/dev/null
    else
        { grep '^#' "${OUT_PREFIX}.vcf"; grep -v '^#' "${OUT_PREFIX}.vcf" | LC_ALL=C sort -k1,1 -k2,2n; } \
            | bgzip > "${OUT_PREFIX}.gold.vcf.gz"
    fi
    tabix -p vcf "${OUT_PREFIX}.gold.vcf.gz"
    GOLD_VCF="${OUT_PREFIX}.gold.vcf.gz"
fi

# Badread reads from mutated genome
msg "Badread: simulating ${DEPTH}x reads"
badread simulate --seed "$SEED" \
    --reference "$MUT_GENOME" --quantity "${DEPTH}x" \
    --length "${READ_MEAN_LEN},${READ_SD_LEN}" \
    --identity "$READ_IDENTITY" 2>"$OUTDIR/${PREFIX}.badread.log" \
    | gzip > "$READS"

# map, sort, index
msg "minimap2 (-x $MINIMAP_PRESET) -> (4) sort + index"
minimap2 -ax "$MINIMAP_PRESET" -t "$THREADS" "$WORK_REF" "$READS" 2>"$OUTDIR/${PREFIX}.minimap2.log" \
    | samtools sort -@ "$THREADS" -o "$BAM" - 2>/dev/null
samtools index -@ "$THREADS" "$BAM"

# summary
msg "Done. Outputs in $OUTDIR/"
{
    echo "  profile         : $PARAMS_FILE"
    echo "  reference       : $REF_FASTA"
    echo "  working ref     : $WORK_REF  (${GENOME_MB} Mb, $n_keep/$n_total contigs >= ${MIN_CONTIG} bp)"
    echo "  mutated genome  : $MUT_GENOME"
    echo "  gold VCF        : $GOLD_VCF"
    echo "  reads           : $READS  (${DEPTH}x)"
    echo "  BAM             : $BAM  (+ .bai)"
    echo
    echo "  gold variant counts:"
    grep -v '^#' "${OUT_PREFIX}.vcf" | grep -oE 'SVTYPE=[A-Z]+' | sort | uniq -c | sed 's/^/    SV /'
    snp=$(grep -v '^#' "${OUT_PREFIX}.vcf" | grep -cE 'SNP[0-9]+SURVIVOR' || true)
    echo "    SNP  $snp"
    echo
    echo "  mean depth:"
    samtools depth -a "$BAM" 2>/dev/null | awk '{s+=$3; n++} END{ if(n) printf "    %.1fx over %d bp\n", s/n, n }'
} >&2

# ../svtrek audt -b $BAM -v <(zcat $GOLD_VCF) -t $THREADS
# truvari bench -b $GOLD_VCF -c calls.vcf.gz -f $REF_FASTA -o $OUTDIR/truvari
