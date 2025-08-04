
grep '<td class="numeric" bgcolor="[^"]*">[^<]*</td>' $1 | sed -n '1p' | sed -n 's/.*bgcolor="[^"]*">\s*\([^<]*\)\s*<\/td>/\1/p'| tr -d ','| sed 's/^[ \t]*//' >> $2
grep '<td class="numeric" bgcolor="[^"]*">[^<]*</td>' $1 | sed -n '11p' | sed -n 's/.*bgcolor="[^"]*">\s*\([^<]*\)\s*<\/td>/\1/p'| tr -d ','| sed 's/^[ \t]*//' >> $2  
grep '<td class="numeric" bgcolor="[^"]*">[^<]*</td>' $1 | sed -n '13p' | sed -n 's/.*bgcolor="[^"]*">\s*\([^<]*\)\s*<\/td>/\1/p'| tr -d ','| sed 's/^[ \t]*//' >> $2
grep '<td class="numeric" bgcolor="[^"]*">[^<]*</td>' $1 | sed -n '15p' | sed -n 's/.*bgcolor="[^"]*">\s*\([^<]*\)\s*<\/td>/\1/p'| tr -d ','| sed 's/^[ \t]*//' >> $2
grep '<td class="numeric" bgcolor="[^"]*">[^<]*</td>' $1 | sed -n '19p' | sed -n 's/.*bgcolor="[^"]*">\s*\([^<]*\)\s*<\/td>/\1/p'| tr -d ','| sed 's/^[ \t]*//' >> $2
grep '<td class="numeric" bgcolor="[^"]*">[^<]*</td>' $1 | sed -n '21p' | sed -n 's/.*bgcolor="[^"]*">\s*\([^<]*\)\s*<\/td>/\1/p'| tr -d ','| sed 's/^[ \t]*//' >> $2
grep '<td class="numeric" bgcolor="[^"]*">[^<]*</td>' $1 | sed -n '23p' | sed -n 's/.*bgcolor="[^"]*">\s*\([^<]*\)\s*<\/td>/\1/p'| tr -d ','| sed 's/^[ \t]*//' >> $2
grep '<th class="numeric">' $1 | sed -n '3p' | sed -n 's/<th class="numeric">\s*\([^<]*\)\s*<\/th>/\1/p'| tr -d ','| sed 's/^[[:space:]]*//' >> $2
grep 'Missense / Silent ratio:' $1 | sed 's/.*: //'| tr -d ','| sed 's/^[ \t]*//' >> $2
#grep '<td class="numeric" bgcolor="[^"]*">[^<]*</td>' $1 | sed -n '35p' | sed -n 's/.*bgcolor="[^"]*">\s*\([^<]*\)\s*<\/td>/\1/p'| tr -d ','| sed 's/^[ \t]*//' >> $2

python extract_intergenic.py $1 $2
