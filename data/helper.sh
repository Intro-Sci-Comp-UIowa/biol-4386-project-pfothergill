#/bin/bash

echo "Enter text"
read text
i=1
for word in $text
do
	echo "Word No-$i = $word"
	((i=$i+1))
done
