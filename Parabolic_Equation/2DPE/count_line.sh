file="./test.dat"
n=$(wc -l < "$file")
sed -i "1i $n" "$file"
