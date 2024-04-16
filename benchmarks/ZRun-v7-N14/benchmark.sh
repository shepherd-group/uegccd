../uegccd > shelloutput
for file in fort.59 fort.60 Output; do
    echo $file
    diff $file */$file | wc -l
done
