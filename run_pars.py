from parsing_scripts import *


pr_ortho_tsv = r'C:\Users\Professional\PycharmProjects\Protein-crossing\Prot_ortho_results\po_no_sensitive.proteinortho.tsv'
#convert_to_family_file(pr_ortho_tsv).to_csv("micro_families.txt", index=False, sep='\t')
#mafft_prep(pr_ortho_tsv, 15).to_csv("bio_mart_input.txt", index=False, sep='\t')
maffr_res = r'C:\Users\Professional\Downloads\224_maffted.fa'
sign = r'C:\Users\Professional\PycharmProjects\Cafe_parser\Significant_families.txt'
base_count = r'C:\Users\Professional\PycharmProjects\Cafe_parser\Base_count.tab'
base_change = r'C:\Users\Professional\PycharmProjects\Cafe_parser\Base_change.tab'
#a = input_mafft(mafft_prep(pr_ortho_tsv), 1)
#print(mafft_prep(pr_ortho_tsv).head(10))
#брала 0, 20, 90
#replace_names(maffr_res, a)
significant(sign, pr_ortho_tsv)
print(incr_decr(base_count, base_change, sign))

