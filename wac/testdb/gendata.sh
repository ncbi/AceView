
# the test program assumes that the class "lots" contains
# exactly 2000 objects

for a in 0 1 
do
for b in 0 1 2 3 4 5 6 7 8 9 
do
for c in 0 1 2 3 4 5 6 7 8 9 
do
for d in 0 1 2 3 4 5 6 7 8 9 
do
	echo 'lots : '$a$b$c$d
	echo 't '$a$b$c$d
	echo ''
done
done
done
done

for x in 0 1 2 3 4 5 6 7 8 9 
do
	echo 'lots : '000$x
	echo 't1 lots-000'$x
	echo ''
done


for x in 0 1 2 3 4 5 6 7 8 9 
do
	echo 'lots : '001$x
	echo 't2 '$x
	echo ''
done
