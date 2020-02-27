#/bin/bash

for name in Watson Micheal Sinha Amir Lily
do
	if [ $name == 'Amir' ]
	then
	for color in Red Green Blue White
	do
		if [ $color == "Blue" ]
		then
		echo "The favorite color of $name is $color"
		break 
		else
		echo "The current color is $color $name"
		fi
		done 
	else
		echo "$name"
	fi
	done
