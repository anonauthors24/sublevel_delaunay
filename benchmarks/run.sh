if [ $# -eq 0 ]; then
    echo "./run.sh results_dir"
    exit 1
fi

#TMPFILE=$(mktemp -d results.XXXXX)
#TMPFILE=results
TMPFILE=$1
mkdir -p "$TMPFILE"/examples
mkdir -p "$TMPFILE"/outputs

dir_example="$TMPFILE"/examples
dir_output="$TMPFILE"/outputs


echo "Computing all examples:"
directory=./"$dir_example"

# Check if the directory exists
if [ -d "$directory" ]; then
    # Declare an empty array
    files=()

    # Iterate over each file in the directory
    for file in "$directory"/*; do
        # Check if the current item is a file (not a directory)
        if [ -f "$file" ]; then
            # Add the file name to the array
            files+=("$(basename "$file")")
        fi
    done

    # Loop through the array of file names and compute the complex
    for filename in "${files[@]}"; do
        echo "Processing... $filename"
        ./main --no-output ./"$dir_example"/"$filename" > ./"$dir_output"/${filename}.txt
        ./main --no-output --no-delaunay-compare ./"$dir_example"/"$filename" > ./"$dir_output"/${filename}_nocompare.txt
        ./main --no-output --only-complex-size ./"$dir_example"/"$filename" > ./"$dir_output"/${filename}_complex.txt
        ./main --multi-chunk --no-output ./"$dir_example"/"$filename" > ./"$dir_output"/${filename}_multichunk.txt
        ./main --multi-chunk --minpres 1 --no-delaunay-compare  --no-output ./"$dir_example"/"$filename" > ./"$dir_output"/${filename}_minpres1.txt
        ./main --multi-chunk --minpres 0 --no-delaunay-compare --no-output ./"$dir_example"/"$filename" > ./"$dir_output"/${filename}_minpres0.txt
        ./main --multi-chunk --minpres 2 --no-delaunay-compare --no-output ./"$dir_example"/"$filename" > ./"$dir_output"/${filename}_minpres2.txt
    done
else
    echo "Directory not found: $directory"
fi

#python3 ./generate_summaries.py ./"$dir_output" ./"$TMPFILE"/summaries.csv

