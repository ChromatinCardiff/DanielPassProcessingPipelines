
awk '{
    for (i = 1; i <= NF; i++) {
        if ($i ~ /gene_id|gene_name/) {
            printf "%s ", $(i+1)
        }
    }
    print ""
}' Homo_sapiens.GRCh37.70.gtf | sed -e 's/"//g' -e 's/;//g' -e 's/ /\t/' | sort -k1,1 | uniq > Homo_sapiens.GRCh37.70.txt
