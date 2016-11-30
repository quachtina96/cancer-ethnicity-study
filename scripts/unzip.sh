# CD into the parent directory containing the folders 
# organizing the data

# Assumes that the data folder 
cd ..
cd gdc_data

for sampleDir in *; do
	echo $sampleDir
 	cd $sampleDir
 	for subdir in *; do
 		if [ -d "$subdir" ]; then
	 	 	echo $subdir
	 	 	cd $subdir
	 		for file in *; do
	 			if [[ $file == *.gz ]]; then
	 				echo "UNZIPPING FILE"
	 				echo $file
	 				gzip -d $file
	 			fi
	 			if [[ $file == *.maf ]]; then
	 				echo "SAVING *.maf AS *.maf.txt"
	 				echo $file
	 				mv $file "$file.txt"
	 			fi
	 		done
	 		cd ..
	 	fi
 	done
 	echo "---"
 	cd ..
done