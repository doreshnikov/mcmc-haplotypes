python3 ../src/main/kotlin/scoring/graphs/stats.py
python3 ../src/main/kotlin/scoring/graphs/likelihood.py $1 -o history$2
python3 ../src/main/kotlin/scoring/scorer_precision_recall.py $1
python3 ../src/main/kotlin/scoring/scorer_ortools.py $1
