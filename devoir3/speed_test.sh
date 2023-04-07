# based on my project for LINFO1252 because bash is so boring
LC_ALL=C
make
# files=(square.geo tuningFork.geo)
make test_inversion -s
# echo "file,meshSize,time" >&1
echo "meshSize,time" >&1
# for file in "${files[@]}"
# do
for i in {3..9}
do
    for j in {1..5}
    do
        value=$( { time ./test_inversion $1 0.$i $2; } 2>&1 | grep real | awk '{print $NF}')
        echo "0.$i,$value" >&1
    done
done
# done
make clean -s