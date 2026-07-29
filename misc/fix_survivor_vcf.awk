#!/usr/bin/awk -f
# Normalize a SURVIVOR simSV VCF into a valid, tab-delimited VCF that
# bcftools / truvari accept. Fixes SURVIVOR's quirks:
#   1. space-delimited #CHROM line            -> force tabs
#   2. missing SAMPLE column in #CHROM header -> add it
#   3. records missing QUAL/FILTER            -> insert  .  / PASS
#   4. FORMAT of 9 keys with 1 value (GT:...:.. -> 1/1) -> collapse to GT
#   5. undeclared INFO/FILTER ids (dup_num, PRECISE, PASS, ...) -> inject defs
BEGIN{
    OFS="\t"
    # standard defs to guarantee present (id -> full header line)
    def["INFO:PRECISE"]  ="##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">"
    def["INFO:IMPRECISE"]="##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">"
    def["INFO:SVTYPE"]   ="##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"
    def["INFO:SVLEN"]    ="##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">"
    def["INFO:END"]      ="##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">"
    def["INFO:CHR2"]     ="##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate\">"
    def["INFO:SVMETHOD"] ="##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect SV\">"
    def["INFO:dup_num"]  ="##INFO=<ID=dup_num,Number=1,Type=Integer,Description=\"Number of duplicated copies\">"
    def["FILTER:PASS"]   ="##FILTER=<ID=PASS,Description=\"All filters passed\">"
    order="INFO:PRECISE INFO:IMPRECISE INFO:SVTYPE INFO:SVLEN INFO:END INFO:CHR2 INFO:SVMETHOD INFO:dup_num FILTER:PASS"
}
/^##/{
    if(match($0,/^##INFO=<ID=[^,>]+/))   seen["INFO:"   substr($0,RSTART+10,RLENGTH-10)]=1
    if(match($0,/^##FILTER=<ID=[^,>]+/)) seen["FILTER:" substr($0,RSTART+12,RLENGTH-12)]=1
    meta[++m]=$0; next
}
/^#CHROM/ && !done{
    for(i=1;i<=m;i++) print meta[i]
    k=split(order,ord," ")
    for(i=1;i<=k;i++) if(!(ord[i] in seen)) print def[ord[i]]
    print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE"
    done=1; next
}
!/^#/{
    n=split($0,a,/[ \t]+/)
    infoi=0
    for(i=6;i<=n;i++) if(a[i] ~ /=/ || a[i] ~ /(^|;)(PRECISE|IMPRECISE)($|;)/){ infoi=i; break }
    if(infoi==0){ print "WARN: no INFO field, skipped: " $0 > "/dev/stderr"; next }
    gt=(infoi+2<=n)?a[infoi+2]:"1/1"; sub(/:.*/,"",gt); if(gt=="")gt="1/1"
    print a[1],a[2],a[3],a[4],a[5],".","PASS",a[infoi],"GT",gt
}
