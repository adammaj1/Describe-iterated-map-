#!/bin/bash 
 
# script file for BASH 
# which bash
# save this file as m.sh
# chmod +x m.sh
# ./m.sh
# checked in https://www.shellcheck.net/






if ! g++ main.cpp -lm -Wall
	then
    		printf" ERROR: compilation failed !!!!!!\n"
    		exit 1
    	else printf "program compiled without errors  \n"
fi



# for all input txt files in this directory run the program \n"
printf "run the program for input txt file:\n" 
for file in *.txt ; do
	# check if file exist 
  	if [ -f "$file" ]
  		then
  			
  			# check if txt file is input not output
        		## https://linuxize.com/post/how-to-check-if-string-contains-substring-in-bash/
        		if [[ "$file" == *"_out.txt" ]]
        			then
            				printf "\t$file is output file = skipped \n"
            			
        			else
  			
  					# b is name of file without extension
  					b=$(basename "$file" .txt)
  					# run the program using input txt file 
  					./a.out "$file" > "$b""_out.txt"
  					printf "\t$file = input file \t "$b""_out.txt" = output file \n"
  				fi
  		else
  			printf " but there is no input txt files in this directory \n"
  	fi
done

 
echo end
# end
