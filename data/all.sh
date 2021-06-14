for i in `seq 80 -5 40`;
	do ./do_all.sh test.$1.$i last;
done;
