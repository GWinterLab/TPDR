#!/bin/bash



split_gene_annotation_in_vcf () {
    local vcf_input_filename="${1}";
    local vcf_output_filename="${2}";
    
    if [ ${#@} -ne 2 ] ; then
        printf 'Usage: split_gene_annotation_in_vcf.sh input.vcf output.vcf\n';
        return 1;
    fi
    
    awk \
        -F '\t' \
        -v OFS='\t' '
        {
            if ( $1 ~ /^#/ ) {
                # Print header lines.
                print $0;
            } else {
                # Find start and length of CSQ field in column 8.
                match($8, /CSQ=[^;]+/);
                
                if ( RSTART != -1 ) {
                    nbr_splitted_CSQ_fields = split( substr($8, RSTART + 4, RLENGTH - 4), splitted_CSQ_field, ",");
                    
                    if (nbr_splitted_CSQ_fields >= 1) {
                        for ( nbr_splitted_CSQ_field_idx=1; nbr_splitted_CSQ_field_idx <= nbr_splitted_CSQ_fields; nbr_splitted_CSQ_field_idx++ ) {
                            # Print mutation line with info for each gene on a seperate line.
                            print $1, $2, $3, $4, $5, $6, $7, substr($8, 1, RSTART - 1) "CSQ=" splitted_CSQ_field[nbr_splitted_CSQ_field_idx] substr($8, RSTART + RLENGTH), $9, $10, $11;
                        }
                    } else {
                        # CSQ field only contains one gene, so print everything.
                        print $0;
                    }
                } else {
                    # No CSQ field found in column 8, print whole line unmodified.
                    print $0;
                }
            }
        }' \
        "${vcf_input_filename}" \
        > "${vcf_output_filename}"
    return $?
}



split_gene_annotation_in_vcf "${@}"


