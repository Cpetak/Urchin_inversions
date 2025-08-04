
#chr=NW_022145594.1 start=11884704 stop=17611976 actustart=12702886 actustop=16793794
#chr=NW_022145594.1 start=38826129 stop=43049305 actustart=39429440 actustop=42445994
#chr=NW_022145597.1 start=14203516 stop=14314359 actustart=14219351 actustop=14298524
#chr=NW_022145600.1 start=33694860 stop=34730460 actustart=33842803 actustop=34582517
#chr=NW_022145601.1 start=28827224 stop=29689418 actustart=28950395 actustop=29566247
#chr=NW_022145603.1 start=8130688 stop=9254692 actustart=8291260 actustop=9094120
#chr=NW_022145606.1 start=16585555 stop=16762970 actustart=16610900 actustop=16737625
#chr=NW_022145609.1 start=29358050 stop=29845887 actustart=29427741 actustop=29776196
chr=NW_022145610.1 start=30642801 stop=31597195 actustart=30779143 actustop=31460853

cname=610
java -Xmx8g -jar snpEff.jar -v Strongylocentrotus_purpuratus ${chr}_filtered_renamed.vcf > ${cname}_snpeff.vcf
mv snpEff_summary.html ${cname}_snpEff_summary.html
bash get_html_vals.sh ${cname}_snpEff_summary.html ${cname}.csv

bash subset_vcf.sh $chr $actustart $actustop
cname=610_inv
java -Xmx8g -jar snpEff.jar -v Strongylocentrotus_purpuratus ${chr}_${actustart}_${actustop}.vcf > ${cname}_snpeff.vcf
mv snpEff_summary.html ${cname}_snpEff_summary.html
bash get_html_vals.sh ${cname}_snpEff_summary.html ${cname}.csv

