
#1 is chromosome
#2 is original region start
#3 is original region end
#4 is group type, eg homop
#5 is inversion label, eg inv1
#6 is window size

#padding=($3-$2)/5
padding=$(( ($3 - $2) / 5 ))
padding=$(printf "%.0f" "$padding")
#new_start=$2-$padding
new_start=$(( $2 - $padding ))
#new_end=$3+$padding
new_end=$(( $3 + $padding ))

#conda deactivate
#bash calc_LD_subgroups.sh $1 ${2}_${3} $4 $new_start $new_end $5 $6
#conda activate grn

if [ "$7" = "TRUE" ]; then
  bash calc_LD_subgroups2.sh $1 ${2}_${3} $4 $new_start $new_end $5 $6
  bash calc_LD_subgroups3.sh $5 $4 $1 $new_start $new_end $6 
else
  bash calc_LD_subgroups.sh $1 ${2}_${3} $4 $new_start $new_end $5 $6
fi
